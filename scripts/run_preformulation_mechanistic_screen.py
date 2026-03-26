#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from ase import Atoms
from ase.io import write
from ase.optimize import FIRE
from fairchem.core import FAIRChemCalculator
from run_formulation_descriptor_pilot import (
    DEFAULT_BASE_DIR,
    POLYMER_SPECS,
    build_rdkit_molecule,
    choose_device,
    rdkit_to_ase,
)


@dataclass(frozen=True)
class MoleculeSpec:
    name: str
    smiles: str
    charge: int = 0
    spin: int = 1


API_LIBRARY: dict[str, MoleculeSpec] = {
    "paracetamol": MoleculeSpec(
        name="Paracetamol",
        smiles="CC(=O)NC1=CC(O)=CC=C1O",
        charge=0,
        spin=1,
    ),
    "diclofenac_anion": MoleculeSpec(
        name="Diclofenac anion",
        smiles="O=C([O-])C1=CC=CC=C1Nc1c(Cl)cccc1Cl",
        charge=-1,
        spin=1,
    ),
}

SMALL_MOLECULES: dict[str, MoleculeSpec] = {
    "water": MoleculeSpec("Water", "O", 0, 1),
    "ethanol": MoleculeSpec("Ethanol", "CCO", 0, 1),
    "glycerol": MoleculeSpec("Glycerol", "OCC(O)CO", 0, 1),
    "sodium": MoleculeSpec("Sodium", "[Na+]", 1, 1),
    "chloride": MoleculeSpec("Chloride", "[Cl-]", -1, 1),
}

COSOLVENT_PRESETS: dict[str, dict[str, int]] = {
    "none": {},
    "water": {"water": 3},
    "ethanol_water": {"ethanol": 1, "water": 2},
    "glycerol_water": {"glycerol": 1, "water": 2},
}

RISK_THRESHOLDS = {
    "high_barrier_eV": 4.2,
    "medium_barrier_eV": 4.5,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run a checkpoint-based mechanistic preformulation screen on API/polymer/cosolvent candidate clusters."
    )
    parser.add_argument(
        "--candidate-matrix",
        type=Path,
        default=DEFAULT_BASE_DIR / "data" / "preformulation" / "paracetamol_mechanistic_screen_demo.csv",
    )
    parser.add_argument(
        "--checkpoint-path",
        type=Path,
        default=DEFAULT_BASE_DIR / "checkpoints" / "best_inference_ckpt.pt",
    )
    parser.add_argument("--task-name", default="omol")
    parser.add_argument("--device", choices=["auto", "cpu", "cuda"], default="auto")
    parser.add_argument("--relax-steps", type=int, default=20)
    parser.add_argument("--fmax", type=float, default=0.10)
    parser.add_argument("--n-replicates", type=int, default=1)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--ranking-group", choices=["all", "neutral", "ionic"], default="all")
    parser.add_argument("--max-candidates", type=int, default=None)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_BASE_DIR / "results" / "mechanistic_screen_run",
    )
    return parser.parse_args()


def build_atoms_from_spec(spec: MoleculeSpec) -> Atoms:
    mol, _ = build_rdkit_molecule(spec.smiles)
    return rdkit_to_ase(mol, charge=spec.charge, spin=spec.spin)


def translated(atoms: Atoms, shift: np.ndarray) -> Atoms:
    shifted = atoms.copy()
    shifted.set_positions(shifted.get_positions() + shift[None, :])
    shifted.info["charge"] = atoms.info.get("charge", 0)
    shifted.info["spin"] = atoms.info.get("spin", 1)
    return shifted


def component_cache_key(spec: MoleculeSpec) -> str:
    return json.dumps(
        {
            "smiles": spec.smiles,
            "charge": spec.charge,
            "spin": spec.spin,
        },
        sort_keys=True,
    )


def combine_components(components: list[dict[str, Any]]) -> tuple[Atoms, list[dict[str, Any]]]:
    merged_numbers: list[int] = []
    merged_positions: list[np.ndarray] = []
    merged_groups: list[dict[str, Any]] = []
    cursor = 0
    for component in components:
        atoms = component["atoms"]
        n_atoms = len(atoms)
        merged_numbers.extend(atoms.get_atomic_numbers().tolist())
        merged_positions.extend(atoms.get_positions())
        merged_groups.append(
            {
                "role": component["role"],
                "label": component["label"],
                "start": cursor,
                "stop": cursor + n_atoms,
            }
        )
        cursor += n_atoms
    merged = Atoms(numbers=merged_numbers, positions=np.asarray(merged_positions, dtype=float))
    merged.set_pbc(False)
    merged.center(vacuum=8.0)
    merged.info["charge"] = int(sum(int(component["atoms"].info.get("charge", 0)) for component in components))
    merged.info["spin"] = 1
    return merged, merged_groups


def candidate_ranking_group(row: pd.Series) -> str:
    return "ionic" if int(row.get("ion_count", 0) or 0) > 0 or int(row.get("api_charge", 0) or 0) != 0 else "neutral"


def build_component_templates(row: pd.Series) -> list[dict[str, Any]]:
    components: list[dict[str, Any]] = []

    polymer_key = row["polymer_key"]
    if polymer_key not in POLYMER_SPECS:
        raise KeyError(f"Unsupported polymer key: {polymer_key}")
    polymer_spec = POLYMER_SPECS[polymer_key]
    polymer_molecule = MoleculeSpec(polymer_spec.name, polymer_spec.smiles, polymer_spec.charge, polymer_spec.spin)

    api_smiles = str(row.get("api_smiles", "") or "").strip()
    api_name = str(row.get("api_name", "") or "").strip().lower()
    api_count = int(row.get("api_count", 0) or 0)
    if api_count > 0:
        api_spec = _resolve_api_spec(api_name, api_smiles, row)
        for idx in range(api_count):
            components.append(
                {
                    "role": "api",
                    "label": f"api_{idx+1}",
                    "cache_key": component_cache_key(api_spec),
                    "template_atoms": build_atoms_from_spec(api_spec),
                }
            )

    polymer_count = int(row.get("polymer_count", 1) or 1)
    for idx in range(polymer_count):
        components.append(
            {
                "role": "polymer",
                "label": f"polymer_{idx+1}",
                "cache_key": component_cache_key(polymer_molecule),
                "template_atoms": build_atoms_from_spec(polymer_molecule),
            }
        )

    cosolvent = str(row.get("cosolvent_system", "none") or "none").strip().lower()
    for mol_key, count in COSOLVENT_PRESETS.get(cosolvent, {}).items():
        spec = SMALL_MOLECULES[mol_key]
        for idx in range(count):
            components.append(
                {
                    "role": "solvent",
                    "label": f"{mol_key}_{idx+1}",
                    "cache_key": component_cache_key(spec),
                    "template_atoms": build_atoms_from_spec(spec),
                }
            )

    ion_smiles = str(row.get("ion_smiles", "") or "").strip()
    ion_count = int(row.get("ion_count", 0) or 0)
    if ion_count > 0:
        ion_spec = _resolve_ion_spec(row, ion_smiles)
        for idx in range(ion_count):
            components.append(
                {
                    "role": "ion",
                    "label": f"ion_{idx+1}",
                    "cache_key": component_cache_key(ion_spec),
                    "template_atoms": build_atoms_from_spec(ion_spec),
                }
            )

    return components


def random_rotation_matrix(rng: np.random.Generator) -> np.ndarray:
    u1, u2, u3 = rng.random(3)
    q1 = np.sqrt(1.0 - u1) * np.sin(2.0 * np.pi * u2)
    q2 = np.sqrt(1.0 - u1) * np.cos(2.0 * np.pi * u2)
    q3 = np.sqrt(u1) * np.sin(2.0 * np.pi * u3)
    q4 = np.sqrt(u1) * np.cos(2.0 * np.pi * u3)
    return np.array(
        [
            [1.0 - 2.0 * (q3**2 + q4**2), 2.0 * (q2 * q3 - q1 * q4), 2.0 * (q2 * q4 + q1 * q3)],
            [2.0 * (q2 * q3 + q1 * q4), 1.0 - 2.0 * (q2**2 + q4**2), 2.0 * (q3 * q4 - q1 * q2)],
            [2.0 * (q2 * q4 - q1 * q3), 2.0 * (q3 * q4 + q1 * q2), 1.0 - 2.0 * (q2**2 + q3**2)],
        ],
        dtype=float,
    )


def randomly_oriented(atoms: Atoms, rng: np.random.Generator) -> Atoms:
    rotated = atoms.copy()
    positions = rotated.get_positions()
    centered_positions = positions - positions.mean(axis=0, keepdims=True)
    rotation = random_rotation_matrix(rng)
    rotated.set_positions(centered_positions @ rotation.T)
    rotated.info["charge"] = atoms.info.get("charge", 0)
    rotated.info["spin"] = atoms.info.get("spin", 1)
    return rotated


def stochastic_placement_vectors(n_components: int, rng: np.random.Generator, jitter_A: float = 1.25) -> list[np.ndarray]:
    frame_rotation = random_rotation_matrix(rng)
    base_vectors = [vector @ frame_rotation.T for vector in placement_vectors(n_components)]
    vectors: list[np.ndarray] = []
    for idx in rng.permutation(len(base_vectors)):
        jitter = rng.normal(loc=0.0, scale=jitter_A, size=3)
        vectors.append(np.asarray(base_vectors[idx] + jitter, dtype=float))
    return vectors


def prepare_components(row: pd.Series, rng: np.random.Generator) -> list[dict[str, Any]]:
    templates = build_component_templates(row)

    positioned = []
    for component, shift in zip(templates, stochastic_placement_vectors(len(templates), rng)):
        oriented = randomly_oriented(component["template_atoms"], rng)
        positioned.append(
            {
                "role": component["role"],
                "label": component["label"],
                "cache_key": component["cache_key"],
                "isolated_atoms": component["template_atoms"],
                "atoms": translated(oriented, shift),
            }
        )
    return positioned


def _resolve_api_spec(api_name: str, api_smiles: str, row: pd.Series) -> MoleculeSpec:
    if api_name and api_name in API_LIBRARY:
        base = API_LIBRARY[api_name]
        charge = int(row.get("api_charge", base.charge) or base.charge)
        spin = int(row.get("api_spin", base.spin) or base.spin)
        return MoleculeSpec(base.name, base.smiles, charge, spin)
    if not api_smiles:
        raise ValueError("Each candidate with api_count > 0 needs api_smiles or a known api_name.")
    return MoleculeSpec(
        name=row.get("api_name", "API"),
        smiles=api_smiles,
        charge=int(row.get("api_charge", 0) or 0),
        spin=int(row.get("api_spin", 1) or 1),
    )


def _resolve_ion_spec(row: pd.Series, ion_smiles: str) -> MoleculeSpec:
    ion_name = str(row.get("ion_name", "") or "").strip().lower()
    if ion_name and ion_name in SMALL_MOLECULES:
        base = SMALL_MOLECULES[ion_name]
        return MoleculeSpec(base.name, base.smiles, base.charge, base.spin)
    if not ion_smiles:
        raise ValueError("Candidates with ion_count > 0 need ion_smiles or ion_name.")
    return MoleculeSpec(
        name=row.get("ion_name", "Ion"),
        smiles=ion_smiles,
        charge=int(row.get("ion_charge", 0) or 0),
        spin=int(row.get("ion_spin", 1) or 1),
    )


def placement_vectors(n_components: int) -> list[np.ndarray]:
    base = [
        np.array([0.0, 0.0, 0.0]),
        np.array([5.0, 0.0, 0.0]),
        np.array([-5.0, 0.0, 0.0]),
        np.array([0.0, 5.0, 0.0]),
        np.array([0.0, -5.0, 0.0]),
        np.array([5.0, 5.0, 0.0]),
        np.array([-5.0, 5.0, 0.0]),
        np.array([5.0, -5.0, 0.0]),
        np.array([-5.0, -5.0, 0.0]),
    ]
    if n_components <= len(base):
        return base[:n_components]
    extra = []
    radius = 5.0
    for idx in range(n_components - len(base)):
        angle = 2.0 * np.pi * idx / max(1, n_components - len(base))
        extra.append(np.array([radius * np.cos(angle), radius * np.sin(angle), 0.0]))
    return base + extra


def replicate_seed_for_candidate(base_seed: int, candidate_id: Any, replicate_index: int) -> int:
    digest = hashlib.sha256(f"{base_seed}:{candidate_id}:{replicate_index}".encode("utf-8")).digest()
    return int.from_bytes(digest[:8], "big") % (2**32)


def isolated_component_energy(calc: FAIRChemCalculator, atoms: Atoms) -> float:
    if len(atoms) == 1:
        # The checkpoint/inference path used here does not expose `omol` single-atom
        # references through the ASE calculator. For proxy ranking, keep monoatomic
        # ions on a neutral zero reference rather than failing the whole screen.
        return 0.0
    isolated = atoms.copy()
    isolated.calc = calc
    return float(isolated.get_potential_energy())


def sanitized_atoms(atoms: Atoms) -> Atoms:
    clean = Atoms(
        numbers=atoms.get_atomic_numbers(),
        positions=atoms.get_positions(),
        cell=atoms.get_cell(),
        pbc=atoms.get_pbc(),
    )
    clean.info.update(dict(atoms.info))
    return clean


def relax_atoms(calc: FAIRChemCalculator, atoms: Atoms, steps: int, fmax: float) -> tuple[Atoms, float]:
    relaxed = atoms.copy()
    relaxed.calc = calc
    if steps > 0:
        dyn = FIRE(relaxed, logfile=None)
        dyn.run(fmax=fmax, steps=steps)
    energy = float(relaxed.get_potential_energy())
    return relaxed, energy


def heavy_atom_indices(atoms: Atoms, start: int, stop: int) -> np.ndarray:
    numbers = atoms.get_atomic_numbers()[start:stop]
    local = np.where(numbers > 1)[0]
    if len(local) == 0:
        local = np.arange(stop - start)
    return local + start


def pair_min_distance(atoms: Atoms, groups: list[dict[str, Any]], role_a: str, role_b: str) -> float | None:
    idx_a_parts = [heavy_atom_indices(atoms, g["start"], g["stop"]) for g in groups if g["role"] == role_a]
    idx_b_parts = [heavy_atom_indices(atoms, g["start"], g["stop"]) for g in groups if g["role"] == role_b]
    if not idx_a_parts or not idx_b_parts:
        return None
    idx_a = np.concatenate(idx_a_parts)
    idx_b = np.concatenate(idx_b_parts)
    dist = atoms.get_all_distances(mic=False)[np.ix_(idx_a, idx_b)]
    return float(np.min(dist))


def pair_contact_count(atoms: Atoms, groups: list[dict[str, Any]], role_a: str, role_b: str, cutoff: float = 3.5) -> int | None:
    idx_a_parts = [heavy_atom_indices(atoms, g["start"], g["stop"]) for g in groups if g["role"] == role_a]
    idx_b_parts = [heavy_atom_indices(atoms, g["start"], g["stop"]) for g in groups if g["role"] == role_b]
    if not idx_a_parts or not idx_b_parts:
        return None
    idx_a = np.concatenate(idx_a_parts)
    idx_b = np.concatenate(idx_b_parts)
    dist = atoms.get_all_distances(mic=False)[np.ix_(idx_a, idx_b)]
    return int(np.sum(dist < cutoff))


def same_role_min_distance(atoms: Atoms, groups: list[dict[str, Any]], role: str) -> float | None:
    same = [g for g in groups if g["role"] == role]
    if len(same) < 2:
        return None
    minima = []
    distances = atoms.get_all_distances(mic=False)
    for i in range(len(same)):
        idx_i = heavy_atom_indices(atoms, same[i]["start"], same[i]["stop"])
        for j in range(i + 1, len(same)):
            idx_j = heavy_atom_indices(atoms, same[j]["start"], same[j]["stop"])
            minima.append(float(np.min(distances[np.ix_(idx_i, idx_j)])))
    return min(minima) if minima else None


def same_role_contact_count(atoms: Atoms, groups: list[dict[str, Any]], role: str, cutoff: float = 3.5) -> int | None:
    same = [g for g in groups if g["role"] == role]
    if len(same) < 2:
        return None
    distances = atoms.get_all_distances(mic=False)
    total = 0
    for i in range(len(same)):
        idx_i = heavy_atom_indices(atoms, same[i]["start"], same[i]["stop"])
        for j in range(i + 1, len(same)):
            idx_j = heavy_atom_indices(atoms, same[j]["start"], same[j]["stop"])
            total += int(np.sum(distances[np.ix_(idx_i, idx_j)] < cutoff))
    return total


def reactive_risk_from_descriptor(row: pd.Series) -> tuple[str, int]:
    barrier = row.get("descriptor_barrier_proxy_eV")
    if barrier is None or pd.isna(barrier):
        barrier = row.get("barrier_proxy_eV")
    if pd.isna(barrier):
        return "unknown", 0
    if barrier < RISK_THRESHOLDS["high_barrier_eV"]:
        return "high", 1
    if barrier < RISK_THRESHOLDS["medium_barrier_eV"]:
        return "medium", 0
    return "low", 0


def metric_or_zero(value: Any) -> float:
    if value is None or pd.isna(value):
        return 0.0
    return float(value)


def screening_score_proxy(metrics: dict[str, Any]) -> float:
    interaction = metrics.get("interaction_energy_eV")
    api_polymer_contacts = metric_or_zero(metrics.get("api_polymer_contact_count"))
    api_api_contacts = metric_or_zero(metrics.get("api_api_contact_count"))
    polymer_polymer_contacts = metric_or_zero(metrics.get("polymer_polymer_contact_count"))
    hotspot_penalty = 1.0 if metrics.get("reactive_hotspot_flag") else 0.0
    score = 0.0
    if interaction is not None and not pd.isna(interaction):
        score += -float(interaction)
    score += 0.10 * api_polymer_contacts
    score -= 0.10 * api_api_contacts
    score += 0.05 * polymer_polymer_contacts
    score -= 0.20 * hotspot_penalty
    return float(score)


def stage_metrics(atoms: Atoms, groups: list[dict[str, Any]], energy_eV: float, isolated_energy_eV: float) -> dict[str, Any]:
    metrics = {
        "energy_eV": energy_eV,
        "interaction_energy_eV": energy_eV - isolated_energy_eV,
        "api_polymer_min_distance_A": pair_min_distance(atoms, groups, "api", "polymer"),
        "api_polymer_contact_count": pair_contact_count(atoms, groups, "api", "polymer"),
        "api_api_min_distance_A": same_role_min_distance(atoms, groups, "api"),
        "api_api_contact_count": same_role_contact_count(atoms, groups, "api"),
        "polymer_polymer_min_distance_A": same_role_min_distance(atoms, groups, "polymer"),
        "polymer_polymer_contact_count": same_role_contact_count(atoms, groups, "polymer"),
        "solvent_polymer_min_distance_A": pair_min_distance(atoms, groups, "solvent", "polymer"),
        "solvent_api_min_distance_A": pair_min_distance(atoms, groups, "solvent", "api"),
        "ion_polymer_min_distance_A": pair_min_distance(atoms, groups, "ion", "polymer"),
        "ion_api_min_distance_A": pair_min_distance(atoms, groups, "ion", "api"),
    }
    return metrics


def score_components(metrics: dict[str, Any]) -> dict[str, float]:
    interaction = metrics.get("interaction_energy_eV")
    api_polymer_contacts = metric_or_zero(metrics.get("api_polymer_contact_count"))
    api_api_contacts = metric_or_zero(metrics.get("api_api_contact_count"))
    polymer_polymer_contacts = metric_or_zero(metrics.get("polymer_polymer_contact_count"))
    hotspot_penalty = 1.0 if metrics.get("reactive_hotspot_flag") else 0.0
    return {
        "score_interaction_term": float(-metric_or_zero(interaction)),
        "score_api_polymer_term": float(0.10 * api_polymer_contacts),
        "score_api_api_term": float(-0.10 * api_api_contacts),
        "score_polymer_polymer_term": float(0.05 * polymer_polymer_contacts),
        "score_hotspot_term": float(-0.20 * hotspot_penalty),
    }


def candidate_metrics(
    row: pd.Series,
    calc: FAIRChemCalculator,
    descriptor_library: dict[str, dict[str, Any]],
    component_cache: dict[str, float],
    relax_steps: int,
    fmax: float,
    output_dir: Path,
    replicate_index: int,
    replicate_seed: int,
) -> dict[str, Any]:
    rng = np.random.default_rng(replicate_seed)
    components = prepare_components(row, rng)
    complex_initial, groups = combine_components(components)
    initial_with_calc = complex_initial.copy()
    initial_with_calc.calc = calc
    pre_energy = float(initial_with_calc.get_potential_energy())
    relaxed, post_energy = relax_atoms(calc, complex_initial, steps=relax_steps, fmax=fmax)

    isolated_energy = 0.0
    for component in components:
        key = component["cache_key"]
        if key not in component_cache:
            component_cache[key] = isolated_component_energy(calc, component["isolated_atoms"])
        isolated_energy += component_cache[key]

    descriptor = descriptor_library.get(row["polymer_key"], {})
    hotspot_risk, hotspot_flag = reactive_risk_from_descriptor(pd.Series(descriptor))

    candidate_dir = output_dir / "structures" / str(row["candidate_id"]) / f"replicate_{replicate_index:03d}"
    candidate_dir.mkdir(parents=True, exist_ok=True)
    write(str(candidate_dir / "initial.xyz"), sanitized_atoms(complex_initial))
    write(str(candidate_dir / "relaxed.xyz"), sanitized_atoms(relaxed))

    pre_stage = stage_metrics(complex_initial, groups, pre_energy, isolated_energy)
    post_stage = stage_metrics(relaxed, groups, post_energy, isolated_energy)

    metrics = {
        "candidate_id": row["candidate_id"],
        "replicate_index": replicate_index,
        "replicate_seed": replicate_seed,
        "api_name": row.get("api_name"),
        "api_charge": row.get("api_charge"),
        "polymer_key": row.get("polymer_key"),
        "cosolvent_system": row.get("cosolvent_system"),
        "ion_count": row.get("ion_count"),
        "priority_tier": row.get("priority_tier"),
        "isolated_component_energy_eV": isolated_energy,
        "reactive_hotspot_risk": hotspot_risk,
        "reactive_hotspot_flag": hotspot_flag,
        "descriptor_barrier_proxy_eV": descriptor.get("barrier_proxy_eV"),
        "descriptor_stretch_tolerance_A": descriptor.get("stretch_tolerance_A"),
        "descriptor_peak_force_eV_per_A": descriptor.get("peak_force_eV_per_A"),
        "descriptor_hotspot_bond_type": descriptor.get("hotspot_bond_type"),
        "structure_initial_xyz": str(candidate_dir / "initial.xyz"),
        "structure_relaxed_xyz": str(candidate_dir / "relaxed.xyz"),
    }
    metrics.update({f"pre_{key}": value for key, value in pre_stage.items()})
    metrics.update({f"post_{key}": value for key, value in post_stage.items()})
    metrics["relaxed_energy_eV"] = post_stage["energy_eV"]
    metrics["interaction_energy_eV"] = post_stage["interaction_energy_eV"]
    metrics["api_polymer_min_distance_A"] = post_stage["api_polymer_min_distance_A"]
    metrics["api_polymer_contact_count"] = post_stage["api_polymer_contact_count"]
    metrics["api_api_min_distance_A"] = post_stage["api_api_min_distance_A"]
    metrics["api_api_contact_count"] = post_stage["api_api_contact_count"]
    metrics["polymer_polymer_min_distance_A"] = post_stage["polymer_polymer_min_distance_A"]
    metrics["polymer_polymer_contact_count"] = post_stage["polymer_polymer_contact_count"]
    metrics["solvent_polymer_min_distance_A"] = post_stage["solvent_polymer_min_distance_A"]
    metrics["solvent_api_min_distance_A"] = post_stage["solvent_api_min_distance_A"]
    metrics["ion_polymer_min_distance_A"] = post_stage["ion_polymer_min_distance_A"]
    metrics["ion_api_min_distance_A"] = post_stage["ion_api_min_distance_A"]
    pre_score_input = {
        **pre_stage,
        "reactive_hotspot_flag": metrics["reactive_hotspot_flag"],
    }
    metrics["pre_screening_score_proxy"] = screening_score_proxy(pre_score_input)
    post_score_terms = score_components({**post_stage, **metrics})
    metrics.update(post_score_terms)
    metrics["screening_score_proxy"] = float(sum(post_score_terms.values()))
    metrics["score_delta_relaxation"] = metrics["screening_score_proxy"] - metrics["pre_screening_score_proxy"]
    metrics["ranking_group"] = candidate_ranking_group(row)
    return metrics


def load_descriptor_library(base_dir: Path) -> dict[str, dict[str, Any]]:
    df = pd.read_csv(base_dir / "data" / "preformulation" / "polymer_descriptor_library.csv")
    return df.set_index("polymer_key").to_dict(orient="index")


def filter_candidates_by_ranking_group(candidates: pd.DataFrame, ranking_group: str) -> pd.DataFrame:
    annotated = candidates.copy()
    annotated["ranking_group"] = annotated.apply(candidate_ranking_group, axis=1)
    if ranking_group != "all":
        annotated = annotated[annotated["ranking_group"] == ranking_group].copy()
    return annotated.drop(columns=["ranking_group"])


def add_replicate_ranks(per_replicate_df: pd.DataFrame) -> pd.DataFrame:
    if per_replicate_df.empty:
        return per_replicate_df.copy()
    ranked = per_replicate_df.copy()
    ranked["replicate_rank"] = (
        ranked.groupby(["ranking_group", "replicate_index"])["screening_score_proxy"]
        .rank(method="first", ascending=False)
        .astype(int)
    )
    return ranked


def aggregate_replicate_results(per_replicate_df: pd.DataFrame) -> pd.DataFrame:
    if per_replicate_df.empty:
        return pd.DataFrame()

    grouped = per_replicate_df.groupby("candidate_id", sort=False)
    metadata_cols = [
        "api_name",
        "api_charge",
        "polymer_key",
        "cosolvent_system",
        "ion_count",
        "priority_tier",
        "reactive_hotspot_risk",
        "reactive_hotspot_flag",
        "descriptor_barrier_proxy_eV",
        "descriptor_stretch_tolerance_A",
        "descriptor_peak_force_eV_per_A",
        "descriptor_hotspot_bond_type",
        "ranking_group",
    ]
    metadata_df = grouped[metadata_cols].first().reset_index()

    exclude_numeric = {"replicate_index", "replicate_seed"} | set(metadata_cols)
    numeric_cols = [
        column
        for column in per_replicate_df.columns
        if pd.api.types.is_numeric_dtype(per_replicate_df[column]) and column not in exclude_numeric
    ]
    numeric_stats = grouped[numeric_cols].agg(["mean", "std"])
    numeric_stats.columns = [
        f"{column}_{'sd' if statistic == 'std' else statistic}" for column, statistic in numeric_stats.columns
    ]
    numeric_stats = numeric_stats.reset_index().fillna(0.0)

    counts_df = grouped.size().rename("n_replicates_completed").reset_index()
    rank_counts = grouped["replicate_rank"].agg(
        top1_count=lambda values: int((values == 1).sum()),
        top2_count=lambda values: int((values <= 2).sum()),
    )
    rank_counts = rank_counts.reset_index().merge(counts_df, on="candidate_id")
    rank_counts["top1_frequency"] = rank_counts["top1_count"] / rank_counts["n_replicates_completed"]
    rank_counts["top2_frequency"] = rank_counts["top2_count"] / rank_counts["n_replicates_completed"]

    best_rows = per_replicate_df.loc[grouped["screening_score_proxy"].idxmax()].copy()
    best_rows = best_rows[
        [
            "candidate_id",
            "replicate_index",
            "replicate_seed",
            "screening_score_proxy",
            "structure_initial_xyz",
            "structure_relaxed_xyz",
        ]
    ].rename(
        columns={
            "replicate_index": "best_replicate_index",
            "replicate_seed": "best_replicate_seed",
            "screening_score_proxy": "best_replicate_screening_score_proxy",
            "structure_initial_xyz": "best_structure_initial_xyz",
            "structure_relaxed_xyz": "best_structure_relaxed_xyz",
        }
    )

    aggregated = metadata_df.merge(counts_df, on="candidate_id")
    aggregated = aggregated.merge(numeric_stats, on="candidate_id")
    aggregated = aggregated.merge(rank_counts.drop(columns=["n_replicates_completed"]), on="candidate_id")
    aggregated = aggregated.merge(best_rows, on="candidate_id")
    if "replicate_rank_mean" in aggregated.columns:
        aggregated = aggregated.rename(columns={"replicate_rank_mean": "mean_rank", "replicate_rank_sd": "rank_sd"})
    return aggregated


def markdown_or_empty(df: pd.DataFrame, columns: list[str], empty_message: str) -> str:
    if df.empty:
        return empty_message
    present_columns = [column for column in columns if column in df.columns]
    return df[present_columns].to_markdown(index=False)


def build_summary_markdown(
    aggregated_df: pd.DataFrame,
    candidate_path: Path,
    checkpoint_path: Path,
    n_replicates: int,
    seed: int,
    ranking_group: str,
) -> str:
    top = aggregated_df.sort_values("screening_score_proxy_mean", ascending=False)
    neutral = top[top["ranking_group"] == "neutral"]
    ionic = top[top["ranking_group"] == "ionic"]
    relaxed = bool((aggregated_df["pre_energy_eV_mean"] != aggregated_df["post_energy_eV_mean"]).any())
    confidence_tier = "ionic proxy branch only" if ranking_group == "ionic" else "neutral candidates only"
    lines = [
        "# Mechanistic Screen v1",
        "",
        f"- candidate matrix: `{candidate_path}`",
        f"- checkpoint: `{checkpoint_path}`",
        f"- ranking group filter: `{ranking_group}`",
        f"- replicates per candidate: `{n_replicates}`",
        f"- base seed: `{seed}`",
        f"- current status: `{'relaxed local compatibility screen' if relaxed else 'runnable static proxy screen'}`",
        f"- confidence tier: `{confidence_tier}`",
        "- current use: `shortlist generation, not decision-making`",
        f"- next gate: `{'DFT spot-checks' if relaxed else 'GPU relaxation + DFT spot-checks'}`",
        "- do not over-interpret: `monoatomic ion cases use a zero-reference proxy for isolated ions`",
        "",
        "## Neutral Ranking",
        "",
        markdown_or_empty(
            neutral,
            [
                "candidate_id",
                "api_name",
                "polymer_key",
                "cosolvent_system",
                "n_replicates_completed",
                "screening_score_proxy_mean",
                "screening_score_proxy_sd",
                "post_interaction_energy_eV_mean",
                "post_interaction_energy_eV_sd",
                "api_polymer_contact_count_mean",
                "api_polymer_contact_count_sd",
                "mean_rank",
                "rank_sd",
                "top1_frequency",
                "top2_frequency",
                "reactive_hotspot_risk",
            ],
            "_No neutral candidates in this run._",
        ),
        "",
        "## Ionic Ranking",
        "",
        markdown_or_empty(
            ionic,
            [
                "candidate_id",
                "api_name",
                "polymer_key",
                "cosolvent_system",
                "n_replicates_completed",
                "screening_score_proxy_mean",
                "screening_score_proxy_sd",
                "post_interaction_energy_eV_mean",
                "post_interaction_energy_eV_sd",
                "post_api_polymer_min_distance_A_mean",
                "post_ion_polymer_min_distance_A_mean",
                "mean_rank",
                "rank_sd",
                "top1_frequency",
                "top2_frequency",
                "reactive_hotspot_risk",
            ],
            "_No ionic candidates in this run._",
        ),
        "",
        "## Score Components",
        "",
        markdown_or_empty(
            top,
            [
                "candidate_id",
                "ranking_group",
                "score_interaction_term_mean",
                "score_api_polymer_term_mean",
                "score_api_api_term_mean",
                "score_polymer_polymer_term_mean",
                "score_hotspot_term_mean",
                "screening_score_proxy_mean",
                "screening_score_proxy_sd",
            ],
            "_No candidates in this run._",
        ),
        "",
        "## Notes",
        "",
        "- This is a local-cluster proxy screen, not a rigorous free-energy calculation.",
        "- `screening_score_proxy` is now decomposed into explicit interaction, aggregation, cohesion, and hotspot terms and reported as replicate mean plus standard deviation.",
        "- Use the aggregated neutral ranking, not a single replicate, to choose the DFT spot-check set.",
        "- Neutral and ionic candidates are separated because the ionic branch is not yet cross-comparable to the neutral branch.",
        "- HPMC and Carbopol-like proxies are not yet wired into this first runnable screen; current supported polymer keys come from the existing descriptor workflow.",
    ]
    return "\n".join(lines) + "\n"


def main() -> int:
    args = parse_args()
    if args.n_replicates < 1:
        raise ValueError("--n-replicates must be at least 1.")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    calc = FAIRChemCalculator.from_model_checkpoint(
        str(args.checkpoint_path),
        task_name=args.task_name,
        device=choose_device(args.device),
    )

    descriptor_library = load_descriptor_library(DEFAULT_BASE_DIR)
    candidates = filter_candidates_by_ranking_group(pd.read_csv(args.candidate_matrix), args.ranking_group)
    if args.max_candidates is not None:
        candidates = candidates.head(args.max_candidates).copy()
    if candidates.empty:
        raise ValueError(f"No candidates remain after applying ranking-group filter `{args.ranking_group}`.")

    component_cache: dict[str, float] = {}
    rows = []
    for replicate_index in range(1, args.n_replicates + 1):
        for _, row in candidates.iterrows():
            rows.append(
                candidate_metrics(
                    row,
                    calc=calc,
                    descriptor_library=descriptor_library,
                    component_cache=component_cache,
                    relax_steps=args.relax_steps,
                    fmax=args.fmax,
                    output_dir=args.output_dir,
                    replicate_index=replicate_index,
                    replicate_seed=replicate_seed_for_candidate(args.seed, row["candidate_id"], replicate_index),
                )
            )

    per_replicate_df = add_replicate_ranks(pd.DataFrame(rows))
    per_replicate_df = per_replicate_df.sort_values(
        ["replicate_index", "ranking_group", "replicate_rank", "candidate_id"],
        ascending=[True, True, True, True],
    )
    aggregated_df = aggregate_replicate_results(per_replicate_df).sort_values(
        "screening_score_proxy_mean", ascending=False
    )
    aggregated_df["shortlist_rank"] = np.arange(1, len(aggregated_df) + 1)
    neutral_df = aggregated_df[aggregated_df["ranking_group"] == "neutral"].copy()
    ionic_df = aggregated_df[aggregated_df["ranking_group"] == "ionic"].copy()
    neutral_df["shortlist_rank_neutral"] = np.arange(1, len(neutral_df) + 1)
    ionic_df["shortlist_rank_ionic"] = np.arange(1, len(ionic_df) + 1)

    per_replicate_csv = args.output_dir / "per_replicate.csv"
    aggregated_csv = args.output_dir / "aggregated.csv"
    aggregated_neutral_csv = args.output_dir / "aggregated_neutral.csv"
    aggregated_ionic_csv = args.output_dir / "aggregated_ionic.csv"
    results_csv = args.output_dir / "mechanistic_screen_results.csv"
    results_neutral_csv = args.output_dir / "mechanistic_screen_results_neutral.csv"
    results_ionic_csv = args.output_dir / "mechanistic_screen_results_ionic.csv"
    summary_md = args.output_dir / "summary.md"
    config_json = args.output_dir / "run_config.json"

    per_replicate_df.to_csv(per_replicate_csv, index=False)
    aggregated_df.to_csv(aggregated_csv, index=False)
    neutral_df.to_csv(aggregated_neutral_csv, index=False)
    ionic_df.to_csv(aggregated_ionic_csv, index=False)
    aggregated_df.to_csv(results_csv, index=False)
    neutral_df.to_csv(results_neutral_csv, index=False)
    ionic_df.to_csv(results_ionic_csv, index=False)
    summary_md.write_text(
        build_summary_markdown(
            aggregated_df,
            args.candidate_matrix,
            args.checkpoint_path,
            n_replicates=args.n_replicates,
            seed=args.seed,
            ranking_group=args.ranking_group,
        )
    )
    config_json.write_text(
        json.dumps(
            {
                "candidate_matrix": str(args.candidate_matrix),
                "checkpoint_path": str(args.checkpoint_path),
                "task_name": args.task_name,
                "device": choose_device(args.device),
                "relax_steps": args.relax_steps,
                "fmax": args.fmax,
                "n_replicates": args.n_replicates,
                "seed": args.seed,
                "ranking_group": args.ranking_group,
            },
            indent=2,
        )
        + "\n"
    )

    print(
        json.dumps(
            {
                "per_replicate_csv": str(per_replicate_csv),
                "aggregated_csv": str(aggregated_csv),
                "aggregated_neutral_csv": str(aggregated_neutral_csv),
                "aggregated_ionic_csv": str(aggregated_ionic_csv),
                "results_csv": str(results_csv),
                "results_neutral_csv": str(results_neutral_csv),
                "results_ionic_csv": str(results_ionic_csv),
                "summary_md": str(summary_md),
                "config_json": str(config_json),
                "n_candidates": int(len(aggregated_df)),
                "n_replicates": int(args.n_replicates),
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
