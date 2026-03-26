#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from collections import deque
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Optional

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from ase import Atoms
from ase.io import write
from fairchem.core import FAIRChemCalculator
from rdkit import Chem
from rdkit.Chem import AllChem


DEFAULT_BASE_DIR = Path(__file__).resolve().parents[1]


@dataclass(frozen=True)
class PolymerSpec:
    name: str
    family: str
    smiles: str
    charge: int
    spin: int
    note: str
    fragment_id: str
    repeat_units: int
    end_capping_strategy: str
    probe_rationale: str
    probe_atom_indices: Optional[tuple[int, int]] = None


POLYMER_SPECS: dict[str, PolymerSpec] = {
    "peg": PolymerSpec(
        name="PEG/PEO trimer-like fragment",
        family="polyether",
        smiles="COCC[O:1][CH2:2]COCCOC",
        charge=0,
        spin=1,
        note="Representative capped PEG-like oligomer with a central C-O backbone probe.",
        fragment_id="peg_trimer_dimethyl_capped",
        repeat_units=3,
        end_capping_strategy="Dimethyl-capped oligomer ends",
        probe_rationale="Central ether backbone O-C bond chosen as a representative polyether scission motif.",
    ),
    "pvp": PolymerSpec(
        name="PVP trimer-like fragment",
        family="vinyl lactam polymer",
        smiles="CC(N1CCCC1=O)[CH2:1][CH:2](N1CCCC1=O)CC(N1CCCC1=O)C",
        charge=0,
        spin=1,
        note="Representative capped PVP-like oligomer with a central backbone C-C probe.",
        fragment_id="pvp_trimer_methyl_capped",
        repeat_units=3,
        end_capping_strategy="Methyl-capped oligomer ends",
        probe_rationale="Central backbone C-C bond chosen to emulate chain scission in the vinyl polymer backbone.",
    ),
    "eudragit_l100_like": PolymerSpec(
        name="Eudragit L100-like fragment",
        family="methacrylate copolymer",
        smiles="C[C](C(=O)O)([CH2:1])[C:2](C)(C(=O)OC)CC(C)(C(=O)O)C",
        charge=0,
        spin=1,
        note="Simplified methacrylic-acid / methyl-methacrylate-like fragment with a central backbone C-C probe.",
        fragment_id="eudragit_l100_like_trimer_simplified",
        repeat_units=3,
        end_capping_strategy="Small alkyl-capped simplified copolymer fragment",
        probe_rationale="Use the central backbone C-C linkage between the left and middle methacrylate-like units.",
        probe_atom_indices=(1, 6),
    ),
    "gelatin": PolymerSpec(
        name="Gelatin-like peptide fragment",
        family="polypeptide",
        smiles="CC(=O)NC[C:1](=O)[N:2]CC(=O)NC",
        charge=0,
        spin=1,
        note="Simplified capped peptide proxy for a gelatin-family backbone scission probe.",
        fragment_id="gelatin_like_ac_glygly_nme_proxy",
        repeat_units=2,
        end_capping_strategy="N-acetyl / N-methylamide caps",
        probe_rationale="Central peptide C-N bond chosen as a representative backbone scission motif for a gelatin-like proteinaceous polymer proxy.",
    ),
    "pectin": PolymerSpec(
        name="Pectin-like digalacturonic acid dimer proxy",
        family="polygalacturonate polysaccharide",
        smiles="O=C(O)[C@H]1O[C@H]([O:1][C@@H:2]2[C@H](O)[C@@H](O)C(O)O[C@@H]2C(=O)O)[C@H](O)[C@@H](O)[C@H]1O",
        charge=0,
        spin=1,
        note="Simplified digalacturonic-acid-like dimer proxy for a pectin-family backbone scission probe.",
        fragment_id="pectin_like_digalacturonic_dimer_proxy",
        repeat_units=2,
        end_capping_strategy="Native dimer termini of the digalacturonic-acid proxy; no extra alkyl caps",
        probe_rationale="Central glycosidic O-C bond chosen as a first pectin-like backbone scission motif under the standardized probe.",
    ),
    "hpc": PolymerSpec(
        name="HPC-like hydroxypropyl cellobiose proxy",
        family="cellulose ether",
        smiles="CC(O)CO[C@H]1[C@H]([O:1][C@@H:2]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@H](CO)O1",
        charge=0,
        spin=1,
        note="Simplified hydroxypropyl-cellulose-like dimer proxy with one representative hydroxypropyl substitution.",
        fragment_id="hpc_like_hydroxypropyl_cellobiose_proxy",
        repeat_units=2,
        end_capping_strategy="Cellobiose-like dimer with one fixed hydroxypropyl substitution pattern",
        probe_rationale="Central glycosidic O-C bond chosen to keep the first cellulose-ether probe aligned with the pectin and xanthan scans.",
    ),
    "xanthan": PolymerSpec(
        name="Xanthan-like branched glucan proxy",
        family="branched polysaccharide",
        smiles="OC[C@H]1O[C@@H]([O:1][C@H:2]2[C@H](O)[C@@H](O[C@@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@H]3O)C(O)O[C@@H]2CO)[C@H](O)[C@@H](O)[C@H]1O",
        charge=0,
        spin=1,
        note="Simplified branched glucan proxy intended to capture xanthan-like backbone-plus-branch topology without claiming a full exact repeat unit.",
        fragment_id="xanthan_like_branched_glucan_proxy",
        repeat_units=3,
        end_capping_strategy="Minimal branched saccharide proxy; no extra alkyl caps",
        probe_rationale="Backbone glycosidic O-C bond chosen so the first xanthan-like descriptor stays comparable to the other polysaccharide families.",
    ),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build representative formulation-polymer fragments and run standardized mapper bond-stretch scans."
    )
    parser.add_argument("--base-dir", type=Path, default=DEFAULT_BASE_DIR)
    parser.add_argument(
        "--checkpoint-path",
        type=Path,
        default=DEFAULT_BASE_DIR / "checkpoints" / "best_inference_ckpt.pt",
    )
    parser.add_argument(
        "--polymers",
        nargs="+",
        default=["peg", "pvp", "eudragit_l100_like"],
        choices=sorted(POLYMER_SPECS),
    )
    parser.add_argument("--task-name", default="omol")
    parser.add_argument("--device", choices=["auto", "cpu", "cuda"], default="auto")
    parser.add_argument("--delta-step", type=float, default=0.1, help="Bond stretch increment in A.")
    parser.add_argument("--max-delta", type=float, default=2.0, help="Maximum bond stretch in A.")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Defaults to results/formulation_descriptor_pilot.",
    )
    return parser.parse_args()


def choose_device(requested: str) -> str:
    if requested in {"cpu", "cuda"}:
        return requested
    import torch

    return "cuda" if torch.cuda.is_available() else "cpu"


def rdkit_warning_banner() -> str:
    return (
        "RDKit is being used in the project venv. If you later want cleaner imports, "
        "pin NumPy<2 in a dedicated fragment-building environment."
    )


def build_rdkit_molecule(smiles: str) -> tuple[Chem.Mol, dict[int, int]]:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    mol = Chem.AddHs(mol)
    status = AllChem.EmbedMolecule(mol, randomSeed=20260325)
    if status != 0:
        raise RuntimeError(f"RDKit embedding failed with status {status}")
    AllChem.MMFFOptimizeMolecule(mol)
    mapped_atoms = {
        atom.GetAtomMapNum(): atom.GetIdx()
        for atom in mol.GetAtoms()
        if atom.GetAtomMapNum() > 0
    }
    if 1 not in mapped_atoms or 2 not in mapped_atoms:
        # Some fragments use explicit heavy-atom indices instead of atom-map labels.
        return mol, mapped_atoms
    return mol, mapped_atoms


def rdkit_to_ase(mol: Chem.Mol, charge: int, spin: int) -> Atoms:
    conf = mol.GetConformer()
    positions = []
    numbers = []
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        positions.append([pos.x, pos.y, pos.z])
        numbers.append(atom.GetAtomicNum())
    atoms = Atoms(numbers=numbers, positions=np.asarray(positions, dtype=float))
    atoms.set_pbc(False)
    atoms.center(vacuum=8.0)
    atoms.info["charge"] = charge
    atoms.info["spin"] = spin
    return atoms


def bond_graph(mol: Chem.Mol) -> dict[int, set[int]]:
    graph: dict[int, set[int]] = {atom.GetIdx(): set() for atom in mol.GetAtoms()}
    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        graph[i].add(j)
        graph[j].add(i)
    return graph


def component_after_cut(graph: dict[int, set[int]], atom_i: int, atom_j: int) -> list[int]:
    trimmed = {node: set(neigh) for node, neigh in graph.items()}
    trimmed[atom_i].discard(atom_j)
    trimmed[atom_j].discard(atom_i)
    seen = {atom_j}
    queue: deque[int] = deque([atom_j])
    while queue:
        node = queue.popleft()
        for neigh in trimmed[node]:
            if neigh not in seen:
                seen.add(neigh)
                queue.append(neigh)
    if atom_i in seen:
        raise RuntimeError(
            "Target bond does not split the fragment into two components. Choose a non-cyclic backbone bond."
        )
    return sorted(seen)


def resolve_probe_atoms(
    spec: PolymerSpec,
    graph: dict[int, set[int]],
    mapped_atoms: dict[int, int],
) -> tuple[int, int]:
    if spec.probe_atom_indices is not None:
        atom_i, atom_j = spec.probe_atom_indices
    else:
        if 1 not in mapped_atoms or 2 not in mapped_atoms:
            raise RuntimeError(
                f"{spec.name} does not define a usable probe bond. "
                "Provide adjacent atom-map labels 1/2 or explicit probe_atom_indices."
            )
        atom_i, atom_j = mapped_atoms[1], mapped_atoms[2]

    if atom_j not in graph[atom_i]:
        raise RuntimeError(
            f"{spec.name} probe atoms ({atom_i}, {atom_j}) are not directly bonded. "
            "Choose a bonded pair for the rigid bond-stretch scan."
        )
    return atom_i, atom_j


def stretch_frames(
    atoms: Atoms,
    moving_atoms: list[int],
    atom_i: int,
    atom_j: int,
    delta_values: list[float],
) -> list[tuple[float, Atoms]]:
    base_positions = atoms.get_positions()
    axis = base_positions[atom_j] - base_positions[atom_i]
    axis_norm = np.linalg.norm(axis)
    if axis_norm == 0:
        raise RuntimeError("Target bond has zero length in the optimized fragment.")
    axis = axis / axis_norm

    frames: list[tuple[float, Atoms]] = []
    for delta in delta_values:
        stretched = atoms.copy()
        new_positions = base_positions.copy()
        new_positions[moving_atoms] += axis * float(delta)
        stretched.set_positions(new_positions)
        stretched.info["charge"] = atoms.info["charge"]
        stretched.info["spin"] = atoms.info["spin"]
        frames.append((float(delta), stretched))
    return frames


def evaluate_scan(
    calc: FAIRChemCalculator,
    frames: list[tuple[float, Atoms]],
    atom_i: int,
    atom_j: int,
) -> list[dict]:
    rows = []
    pred_energies = []
    for delta, atoms in frames:
        test_atoms = atoms.copy()
        test_atoms.calc = calc
        pred_energy = float(test_atoms.get_potential_energy())
        pred_forces = np.asarray(test_atoms.get_forces(), dtype=float)
        bond_length = float(test_atoms.get_distance(atom_i, atom_j, mic=False))
        atom_i_force = float(np.linalg.norm(pred_forces[atom_i]))
        atom_j_force = float(np.linalg.norm(pred_forces[atom_j]))
        mean_force = 0.5 * (atom_i_force + atom_j_force)
        pred_energies.append(pred_energy)
        rows.append(
            {
                "delta_A": delta,
                "bond_length_A": bond_length,
                "pred_energy_eV": pred_energy,
                "atom_i_pred_force_norm_eV_per_A": atom_i_force,
                "atom_j_pred_force_norm_eV_per_A": atom_j_force,
                "mean_target_force_norm_eV_per_A": mean_force,
            }
        )

    ref_energy = pred_energies[0]
    for row in rows:
        row["relative_energy_eV"] = float(row["pred_energy_eV"] - ref_energy)
    return rows


def write_scan_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def plot_scan(rows: list[dict], title: str, output_prefix: Path) -> None:
    bond_length = [float(row["bond_length_A"]) for row in rows]
    rel_energy = [float(row["relative_energy_eV"]) for row in rows]
    mean_force = [float(row["mean_target_force_norm_eV_per_A"]) for row in rows]

    fig, axes = plt.subplots(2, 1, figsize=(7.5, 6.2), dpi=220, sharex=True)
    axes[0].plot(bond_length, rel_energy, marker="o", linewidth=2.0, color="#9b2226")
    axes[0].set_ylabel("Relative energy (eV)")
    axes[0].set_title(title)
    axes[0].grid(True, alpha=0.25)

    axes[1].plot(bond_length, mean_force, marker="s", linewidth=2.0, color="#005f73")
    axes[1].set_xlabel("Bond length (A)")
    axes[1].set_ylabel("Mean target force (eV/A)")
    axes[1].grid(True, alpha=0.25)

    fig.tight_layout()
    fig.savefig(output_prefix.with_suffix(".png"), bbox_inches="tight")
    fig.savefig(output_prefix.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def extract_descriptors(
    rows: list[dict],
    atom_i: int,
    atom_j: int,
    atoms: Atoms,
) -> dict:
    best_energy_row = max(rows, key=lambda row: float(row["relative_energy_eV"]))
    best_force_row = max(rows, key=lambda row: float(row["mean_target_force_norm_eV_per_A"]))
    symbols = atoms.get_chemical_symbols()
    return {
        "barrier_proxy_eV": float(best_energy_row["relative_energy_eV"]),
        "hotspot_bond_type": f"{symbols[atom_i]}-{symbols[atom_j]}",
        "stretch_tolerance_A": float(best_force_row["bond_length_A"]),
        "peak_force_eV_per_A": float(best_force_row["mean_target_force_norm_eV_per_A"]),
        "target_atom_i": int(atom_i),
        "target_atom_j": int(atom_j),
    }


def write_summary_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    base_dir = args.base_dir
    output_dir = args.output_dir or (base_dir / "results" / "formulation_descriptor_pilot")
    output_dir.mkdir(parents=True, exist_ok=True)

    calc = FAIRChemCalculator.from_model_checkpoint(
        str(args.checkpoint_path),
        task_name=args.task_name,
        device=choose_device(args.device),
    )

    delta_values = list(np.arange(0.0, args.max_delta + 1e-8, args.delta_step))
    summary_rows = []

    for polymer_key in args.polymers:
        spec = POLYMER_SPECS[polymer_key]
        polymer_dir = output_dir / polymer_key
        polymer_dir.mkdir(parents=True, exist_ok=True)

        mol, mapped_atoms = build_rdkit_molecule(spec.smiles)
        atoms = rdkit_to_ase(mol, charge=spec.charge, spin=spec.spin)
        graph = bond_graph(mol)
        atom_i, atom_j = resolve_probe_atoms(spec, graph, mapped_atoms)
        moving_atoms = component_after_cut(graph, atom_i, atom_j)
        frames = stretch_frames(atoms, moving_atoms, atom_i, atom_j, delta_values)
        scan_rows = evaluate_scan(calc, frames, atom_i, atom_j)

        optimized_xyz = polymer_dir / "optimized_fragment.xyz"
        scan_xyz = polymer_dir / "bond_scan_frames.xyz"
        scan_csv = polymer_dir / "bond_scan_metrics.csv"
        descriptors_json = polymer_dir / "descriptors.json"
        plot_prefix = polymer_dir / "bond_scan_plot"

        write(str(optimized_xyz), atoms)
        write(str(scan_xyz), [frame for _, frame in frames])
        write_scan_csv(scan_csv, scan_rows)
        plot_scan(scan_rows, f"{spec.name}: standardized bond-stretch probe", plot_prefix)

        descriptors = extract_descriptors(scan_rows, atom_i, atom_j, atoms)
        payload = {
            "polymer_key": polymer_key,
            "polymer_name": spec.name,
            "family": spec.family,
            "note": spec.note,
            "warning": rdkit_warning_banner(),
            "spec": asdict(spec),
            "protocol": {
                "geometry_handling": "Rigid bond-stretch scan. One graph component is translated along the target bond axis; no constrained relaxation is performed.",
                "scan_settings": {
                    "starting_bond_length_A": float(scan_rows[0]["bond_length_A"]),
                    "final_bond_length_A": float(scan_rows[-1]["bond_length_A"]),
                    "step_size_A": float(args.delta_step),
                    "num_frames": len(scan_rows),
                },
                "probe_bond": {
                    "atom_i": int(atom_i),
                    "atom_j": int(atom_j),
                    "bond_type": descriptors["hotspot_bond_type"],
                    "selection_rationale": spec.probe_rationale,
                },
                "inference": {
                    "checkpoint_path": str(args.checkpoint_path),
                    "task_name": args.task_name,
                    "device": choose_device(args.device),
                    "script": str(Path(__file__).resolve()),
                },
                "descriptor_definitions": {
                    "barrier_proxy_eV": "Maximum predicted relative energy along the bond-stretch scan.",
                    "stretch_tolerance_A": "Bond length at the maximum mean target force along the scan.",
                    "peak_force_eV_per_A": "Maximum mean force norm across the two probe atoms.",
                    "hotspot_bond_type": "Element-pair label of the explicitly probed bond.",
                },
            },
            "descriptors": descriptors,
            "outputs": {
                "optimized_xyz": str(optimized_xyz),
                "scan_xyz": str(scan_xyz),
                "scan_csv": str(scan_csv),
                "plot_png": str(plot_prefix.with_suffix(".png")),
                "plot_pdf": str(plot_prefix.with_suffix(".pdf")),
            },
        }
        descriptors_json.write_text(json.dumps(payload, indent=2))

        summary_rows.append(
            {
                "polymer_key": polymer_key,
                "polymer_name": spec.name,
                "family": spec.family,
                "barrier_proxy_eV": descriptors["barrier_proxy_eV"],
                "hotspot_bond_type": descriptors["hotspot_bond_type"],
                "stretch_tolerance_A": descriptors["stretch_tolerance_A"],
                "peak_force_eV_per_A": descriptors["peak_force_eV_per_A"],
                "target_atom_i": descriptors["target_atom_i"],
                "target_atom_j": descriptors["target_atom_j"],
                "status": "done",
                "scan_csv": str(scan_csv),
                "descriptors_json": str(descriptors_json),
            }
        )

    summary_csv = output_dir / "descriptor_summary.csv"
    summary_json = output_dir / "descriptor_summary.json"
    write_summary_csv(summary_csv, summary_rows)
    summary_json.write_text(json.dumps(summary_rows, indent=2))
    print(json.dumps({"output_dir": str(output_dir), "polymers": args.polymers}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
