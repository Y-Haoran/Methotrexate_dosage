#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors

from chemistry_registry import API_LIBRARY, SMALL_MOLECULES


DEFAULT_BASE_DIR = Path(__file__).resolve().parents[1]


ACIDIC_PATTERNS = [
    Chem.MolFromSmarts("[CX3](=O)[OX2H1]"),
    Chem.MolFromSmarts("[$([SX4](=O)(=O)[OX2H1])]"),
    Chem.MolFromSmarts("[nH]1cccc1"),
]

BASIC_PATTERNS = [
    Chem.MolFromSmarts("[NX3;H2,H1,H0;!$(NC=O)]"),
    Chem.MolFromSmarts("[nX2;$([nH0;+0])][c,n]"),
]


def clip(value: float, low: float = 0.0, high: float = 1.0) -> float:
    return max(low, min(high, float(value)))


def slugify(text: str) -> str:
    return "".join(char.lower() if char.isalnum() else "_" for char in text).strip("_")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate API descriptors, assign formulation-relevant classes, "
            "recommend formulation families, and write a downstream mechanistic screen matrix."
        )
    )
    parser.add_argument(
        "--input-json",
        type=Path,
        default=DEFAULT_BASE_DIR / "data" / "preformulation" / "api_family_recommendation_input_template.json",
    )
    parser.add_argument(
        "--prior-table",
        type=Path,
        default=DEFAULT_BASE_DIR / "data" / "preformulation" / "family_recommendation_priors.csv",
    )
    parser.add_argument(
        "--polymer-library",
        type=Path,
        default=DEFAULT_BASE_DIR / "data" / "preformulation" / "polymer_descriptor_library.csv",
    )
    parser.add_argument("--max-families-per-state", type=int, default=None)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_BASE_DIR / "results" / "api_family_recommendation_run",
    )
    return parser.parse_args()


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def resolve_parent_smiles(payload: dict[str, Any]) -> str:
    explicit_smiles = str(payload.get("api_smiles", "") or "").strip()
    if explicit_smiles:
        if Chem.MolFromSmiles(explicit_smiles) is None:
            raise ValueError(f"Invalid api_smiles in {payload.get('api_name', 'input payload')}.")
        return explicit_smiles

    api_name = str(payload.get("api_name", "") or "").strip().lower()
    if api_name in API_LIBRARY:
        return API_LIBRARY[api_name].smiles
    raise ValueError("Input needs api_smiles or a known api_name from chemistry_registry.API_LIBRARY.")


def normalize_state_options(payload: dict[str, Any], parent_smiles: str) -> list[dict[str, Any]]:
    state_options = payload.get("state_options") or []
    if not state_options:
        mol = Chem.MolFromSmiles(parent_smiles)
        default_charge = Chem.GetFormalCharge(mol)
        state_options = [
            {
                "state_id": "neutral" if default_charge == 0 else f"charge_{default_charge:+d}".replace("+", "plus_"),
                "label": "neutral" if default_charge == 0 else f"charge {default_charge:+d}",
                "formal_charge": default_charge,
                "api_spin": 1,
                "state_smiles": None,
                "include_in_screen": True,
                "screen_priority": "primary",
            }
        ]

    normalized: list[dict[str, Any]] = []
    for index, state in enumerate(state_options, start=1):
        state_id = str(state.get("state_id", f"state_{index}") or f"state_{index}").strip()
        state_smiles = str(state.get("state_smiles", "") or "").strip() or parent_smiles
        if Chem.MolFromSmiles(state_smiles) is None:
            raise ValueError(f"Invalid state_smiles for state '{state_id}'.")
        normalized.append(
            {
                "state_id": state_id,
                "label": str(state.get("label", state_id) or state_id),
                "formal_charge": int(state.get("formal_charge", 0) or 0),
                "api_spin": int(state.get("api_spin", 1) or 1),
                "state_smiles": state_smiles,
                "include_in_screen": bool(state.get("include_in_screen", True)),
                "screen_priority": str(state.get("screen_priority", "primary") or "primary"),
            }
        )
    return normalized


def normalize_salt_options(payload: dict[str, Any]) -> list[dict[str, Any]]:
    salts = []
    for index, salt in enumerate(payload.get("salt_options") or [], start=1):
        ion_name = str(salt.get("ion_name", "") or "").strip().lower()
        ion_smiles = str(salt.get("ion_smiles", "") or "").strip()
        if ion_name in SMALL_MOLECULES and not ion_smiles:
            ion_spec = SMALL_MOLECULES[ion_name]
            ion_smiles = ion_spec.smiles
        if ion_smiles and Chem.MolFromSmiles(ion_smiles) is None:
            raise ValueError(f"Invalid ion_smiles for salt option {salt.get('salt_id', index)}.")
        salts.append(
            {
                "salt_id": str(salt.get("salt_id", f"salt_{index}") or f"salt_{index}"),
                "state_id": str(salt.get("state_id", "") or "").strip(),
                "ion_name": ion_name,
                "ion_smiles": ion_smiles,
                "ion_charge": int(salt.get("ion_charge", 0) or 0),
                "ion_spin": int(salt.get("ion_spin", 1) or 1),
                "ion_count": int(salt.get("ion_count", 0) or 0),
                "include_in_screen": bool(salt.get("include_in_screen", True)),
            }
        )
    return salts


def count_matches(mol: Chem.Mol, patterns: list[Chem.Mol]) -> int:
    count = 0
    for pattern in patterns:
        if pattern is None:
            continue
        count += len(mol.GetSubstructMatches(pattern))
    return count


def descriptor_smiles_note(parent_smiles: str, state_smiles: str, formal_charge: int) -> str:
    if state_smiles != parent_smiles:
        return "state-specific smiles"
    if formal_charge == 0:
        return "parent neutral smiles"
    return "parent smiles with metadata-only charge override"


def compute_api_descriptors(parent_smiles: str, state: dict[str, Any]) -> dict[str, Any]:
    mol = Chem.MolFromSmiles(state["state_smiles"])
    if mol is None:
        raise ValueError(f"Invalid state_smiles for state {state['state_id']}.")

    molecular_weight = float(Descriptors.MolWt(mol))
    heavy_atom_count = int(mol.GetNumHeavyAtoms())
    hetero_atom_count = int(sum(atom.GetAtomicNum() not in (1, 6) for atom in mol.GetAtoms()))
    clogp = float(Crippen.MolLogP(mol))
    tpsa = float(rdMolDescriptors.CalcTPSA(mol))
    hbond_donors = int(Lipinski.NumHDonors(mol))
    hbond_acceptors = int(Lipinski.NumHAcceptors(mol))
    rotatable_bonds = int(Lipinski.NumRotatableBonds(mol))
    aromatic_ring_count = int(rdMolDescriptors.CalcNumAromaticRings(mol))
    acidic_group_count = int(count_matches(mol, ACIDIC_PATTERNS))
    basic_group_count = int(count_matches(mol, BASIC_PATTERNS))
    formal_charge = int(state["formal_charge"])

    flexibility_proxy = clip((rotatable_bonds / max(heavy_atom_count, 1)) * 3.5)
    rigidity_proxy = clip(1.0 - flexibility_proxy)
    polarity_balance_proxy = clip(
        0.55 * min(tpsa / 180.0, 1.0)
        + 0.25 * min((hbond_donors + hbond_acceptors) / 12.0, 1.0)
        + 0.20 * clip((2.5 - max(clogp, -1.0)) / 4.5),
    )
    water_affinity_proxy = clip(
        0.45 * min(tpsa / 180.0, 1.0)
        + 0.25 * min((hbond_donors + hbond_acceptors) / 12.0, 1.0)
        + 0.20 * (1.0 / (1.0 + math.exp(clogp - 1.5)))
        + 0.10 * min(abs(formal_charge), 2),
    )
    self_aggregation_risk_proxy = clip(
        0.45 * min(aromatic_ring_count / 4.0, 1.0)
        + 0.25 * min(max(clogp, 0.0) / 4.0, 1.0)
        + 0.20 * rigidity_proxy
        + 0.10 * (1.0 - water_affinity_proxy),
    )
    crystallization_risk_proxy = clip(
        0.40 * self_aggregation_risk_proxy
        + 0.30 * rigidity_proxy
        + 0.20 * min(heavy_atom_count / 45.0, 1.0)
        + 0.10 * (1.0 - polarity_balance_proxy),
    )

    return {
        "state_id": state["state_id"],
        "state_label": state["label"],
        "formal_charge": formal_charge,
        "api_spin": state["api_spin"],
        "descriptor_source_smiles": state["state_smiles"],
        "descriptor_source_note": descriptor_smiles_note(parent_smiles, state["state_smiles"], formal_charge),
        "molecular_weight": molecular_weight,
        "heavy_atom_count": heavy_atom_count,
        "hetero_atom_fraction": hetero_atom_count / max(heavy_atom_count, 1),
        "cLogP": clogp,
        "TPSA_A2": tpsa,
        "hbond_donors": hbond_donors,
        "hbond_acceptors": hbond_acceptors,
        "rotatable_bonds": rotatable_bonds,
        "aromatic_ring_count": aromatic_ring_count,
        "acidic_group_count": acidic_group_count,
        "basic_group_count": basic_group_count,
        "ionizable_group_count": acidic_group_count + basic_group_count,
        "flexibility_proxy": flexibility_proxy,
        "rigidity_proxy": rigidity_proxy,
        "polarity_balance_proxy": polarity_balance_proxy,
        "water_affinity_proxy": water_affinity_proxy,
        "self_aggregation_risk_proxy": self_aggregation_risk_proxy,
        "crystallization_risk_proxy": crystallization_risk_proxy,
    }


def assign_api_classes(descriptors: dict[str, Any]) -> tuple[str, list[str]]:
    classes: set[str] = set()
    formal_charge = int(descriptors["formal_charge"])
    acidic_group_count = int(descriptors["acidic_group_count"])
    basic_group_count = int(descriptors["basic_group_count"])

    if formal_charge < 0 or acidic_group_count > 0:
        classes.add("acidic_or_anionic")
    if formal_charge > 0 or (basic_group_count > 0 and acidic_group_count == 0):
        classes.add("basic_or_cationic")
    if formal_charge == 0:
        if descriptors["water_affinity_proxy"] >= 0.58 and descriptors["TPSA_A2"] >= 75.0:
            classes.add("hydrophilic_neutral")
        else:
            classes.add("poorly_soluble_neutral")
    if acidic_group_count > 0 and basic_group_count > 0:
        classes.add("amphiphilic")
    elif 0.30 <= descriptors["water_affinity_proxy"] <= 0.80 and (
        descriptors["aromatic_ring_count"] >= 1 or abs(formal_charge) > 0
    ):
        classes.add("amphiphilic")
    if descriptors["rotatable_bonds"] >= 7 and descriptors["TPSA_A2"] >= 80.0:
        classes.add("highly_flexible_polar")
    if descriptors["aromatic_ring_count"] >= 2 and descriptors["self_aggregation_risk_proxy"] >= 0.42:
        classes.add("aromatic_aggregation_prone")
    if not classes:
        classes.add("hydrophilic_neutral" if formal_charge == 0 else "amphiphilic")

    if formal_charge < 0:
        primary_class = "acidic_or_anionic"
    elif formal_charge > 0:
        primary_class = "basic_or_cationic"
    elif "poorly_soluble_neutral" in classes and descriptors["self_aggregation_risk_proxy"] >= 0.45:
        primary_class = "poorly_soluble_neutral"
    elif "hydrophilic_neutral" in classes:
        primary_class = "hydrophilic_neutral"
    elif "highly_flexible_polar" in classes:
        primary_class = "highly_flexible_polar"
    else:
        primary_class = sorted(classes)[0]

    return primary_class, sorted(classes)


def load_prior_table(path: Path) -> pd.DataFrame:
    priors = pd.read_csv(path)
    required = {
        "api_class",
        "formulation_context",
        "route_of_administration",
        "polymer_key",
        "cosolvent_system",
        "charge_scope",
        "base_prior_score",
        "rationale",
    }
    missing = required - set(priors.columns)
    if missing:
        raise ValueError(f"Prior table is missing columns: {sorted(missing)}")
    return priors


def load_polymer_library(path: Path) -> pd.DataFrame:
    library = pd.read_csv(path)
    required = {"polymer_key", "polymer_name", "family"}
    missing = required - set(library.columns)
    if missing:
        raise ValueError(f"Polymer library is missing columns: {sorted(missing)}")
    return library


def matches_scope(prior_value: str, request_value: str) -> bool:
    prior_value = str(prior_value or "").strip().lower()
    request_value = str(request_value or "").strip().lower()
    return prior_value == "any" or prior_value == request_value


def charge_scope_for_state(formal_charge: int) -> str:
    return "neutral" if int(formal_charge) == 0 else "ionic"


def score_prior_match(
    prior: pd.Series,
    state: dict[str, Any],
    descriptors: dict[str, Any],
    primary_class: str,
    payload: dict[str, Any],
) -> float:
    score = float(prior["base_prior_score"])
    if str(prior["api_class"]) == primary_class:
        score += 0.35
    else:
        score += 0.15

    if matches_scope(prior["formulation_context"], payload.get("formulation_context", "")):
        if str(prior["formulation_context"]).strip().lower() != "any":
            score += 0.20
    if matches_scope(prior["route_of_administration"], payload.get("route_of_administration", "")):
        if str(prior["route_of_administration"]).strip().lower() != "any":
            score += 0.10

    processing_context = str(payload.get("processing_context", "") or "").lower()
    polymer_key = str(prior["polymer_key"])
    cosolvent_system = str(prior["cosolvent_system"])

    if "dry" in processing_context and cosolvent_system == "ethanol_water":
        score += 0.10
    if "gel" in processing_context and cosolvent_system == "glycerol_water":
        score += 0.08
    if descriptors["self_aggregation_risk_proxy"] >= 0.55 and polymer_key == "pvp":
        score += 0.20
    if descriptors["water_affinity_proxy"] >= 0.65 and cosolvent_system == "water":
        score += 0.10
    if int(state["formal_charge"]) < 0 and polymer_key == "hpc":
        score += 0.15
    return score


def recommend_families(
    payload: dict[str, Any],
    state: dict[str, Any],
    descriptors: dict[str, Any],
    primary_class: str,
    api_classes: list[str],
    priors: pd.DataFrame,
    polymer_library: pd.DataFrame,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    state_charge_scope = charge_scope_for_state(state["formal_charge"])

    for _, prior in priors.iterrows():
        if str(prior["api_class"]) not in api_classes:
            continue
        if not matches_scope(prior["formulation_context"], payload.get("formulation_context", "")):
            continue
        if not matches_scope(prior["route_of_administration"], payload.get("route_of_administration", "")):
            continue
        prior_charge_scope = str(prior["charge_scope"]).strip().lower()
        if prior_charge_scope not in {"all", state_charge_scope}:
            continue

        rows.append(
            {
                "state_id": state["state_id"],
                "state_label": state["label"],
                "formal_charge": state["formal_charge"],
                "primary_api_class": primary_class,
                "matched_api_class": str(prior["api_class"]),
                "polymer_key": str(prior["polymer_key"]),
                "cosolvent_system": str(prior["cosolvent_system"]),
                "charge_scope": prior_charge_scope,
                "base_prior_score": float(prior["base_prior_score"]),
                "match_score": score_prior_match(prior, state, descriptors, primary_class, payload),
                "rationale": str(prior["rationale"]),
            }
        )

    if not rows:
        return pd.DataFrame()

    frame = pd.DataFrame(rows)
    grouped_rows: list[dict[str, Any]] = []
    for (state_id, polymer_key, cosolvent_system), group in frame.groupby(
        ["state_id", "polymer_key", "cosolvent_system"], sort=False
    ):
        matched_api_classes = sorted(set(group["matched_api_class"].tolist()))
        best_row = group.sort_values("match_score", ascending=False).iloc[0].to_dict()
        best_row["matched_api_classes"] = ",".join(matched_api_classes)
        best_row["api_class_match_count"] = len(matched_api_classes)
        best_row["recommendation_score"] = float(best_row["match_score"]) + 0.08 * (len(matched_api_classes) - 1)
        grouped_rows.append(best_row)

    recommendations = pd.DataFrame(grouped_rows)
    recommendations = recommendations.merge(
        polymer_library[["polymer_key", "polymer_name", "family"]],
        on="polymer_key",
        how="left",
    )
    recommendations = recommendations.sort_values(
        ["state_id", "recommendation_score", "polymer_key", "cosolvent_system"],
        ascending=[True, False, True, True],
    ).reset_index(drop=True)
    recommendations["state_rank"] = recommendations.groupby("state_id").cumcount() + 1
    return recommendations


def choose_salt_for_state(state_id: str, salts: list[dict[str, Any]]) -> dict[str, Any] | None:
    matches = [salt for salt in salts if salt["state_id"] == state_id and salt["include_in_screen"]]
    if not matches:
        return None
    return matches[0]


def build_candidate_matrix(
    payload: dict[str, Any],
    parent_smiles: str,
    state_options: list[dict[str, Any]],
    salts: list[dict[str, Any]],
    recommendations: pd.DataFrame,
    max_families_per_state: int,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    api_name = str(payload.get("api_name", "api") or "api")
    api_slug = slugify(api_name)

    for state in state_options:
        if not state["include_in_screen"]:
            continue
        state_rows = recommendations[recommendations["state_id"] == state["state_id"]].head(max_families_per_state)
        salt = choose_salt_for_state(state["state_id"], salts)
        for index, (_, row) in enumerate(state_rows.iterrows(), start=1):
            candidate_id = f"{api_slug}_{slugify(state['state_id'])}_{row['polymer_key']}_{row['cosolvent_system']}_{index:02d}"
            ionic = int(state["formal_charge"]) != 0
            priority_tier = "primary" if int(row["state_rank"]) <= 2 else "screen"
            notes = (
                f"family_recommendation_score={row['recommendation_score']:.2f}; "
                f"matched_api_classes={row['matched_api_classes']}; "
                f"rationale={row['rationale']}"
            )
            if ionic and state["state_smiles"] == parent_smiles:
                notes += "; warning=replace parent-smiles charge override with an explicit ionized state_smiles before production mechanistic screening"
            rows.append(
                {
                    "candidate_id": candidate_id,
                    "api_name": api_name,
                    "api_smiles": state["state_smiles"],
                    "api_charge": int(state["formal_charge"]),
                    "api_spin": int(state["api_spin"]),
                    "api_count": 1,
                    "polymer_key": row["polymer_key"],
                    "polymer_count": 1,
                    "cosolvent_system": row["cosolvent_system"],
                    "ion_name": salt["ion_name"] if ionic and salt else "",
                    "ion_smiles": salt["ion_smiles"] if ionic and salt else "",
                    "ion_charge": int(salt["ion_charge"]) if ionic and salt else 0,
                    "ion_spin": int(salt["ion_spin"]) if ionic and salt else 1,
                    "ion_count": int(salt["ion_count"]) if ionic and salt else 0,
                    "priority_tier": priority_tier,
                    "notes": notes,
                }
            )
    return pd.DataFrame(rows)


def write_summary(
    path: Path,
    payload: dict[str, Any],
    state_rows: pd.DataFrame,
    recommendations: pd.DataFrame,
    candidate_matrix: pd.DataFrame,
    max_families_per_state: int,
) -> None:
    lines = [
        f"# Family Recommendation Summary: {payload.get('api_name', 'API')}",
        "",
        "## Input Context",
        "",
        f"- formulation_context: `{payload.get('formulation_context', 'unknown')}`",
        f"- route_of_administration: `{payload.get('route_of_administration', 'unknown')}`",
        f"- processing_context: `{payload.get('processing_context', 'unknown')}`",
        f"- pH_context: `{payload.get('pH_context', 'unknown')}`",
        f"- top_k_families_per_state: `{max_families_per_state}`",
        "",
        "## API State Classification",
        "",
    ]

    for _, row in state_rows.iterrows():
        lines.extend(
            [
                f"### {row['state_id']}",
                "",
                f"- label: `{row['state_label']}`",
                f"- formal_charge: `{int(row['formal_charge'])}`",
                f"- primary_api_class: `{row['primary_api_class']}`",
                f"- assigned_api_classes: `{row['assigned_api_classes']}`",
                f"- water_affinity_proxy: `{row['water_affinity_proxy']:.2f}`",
                f"- self_aggregation_risk_proxy: `{row['self_aggregation_risk_proxy']:.2f}`",
                f"- crystallization_risk_proxy: `{row['crystallization_risk_proxy']:.2f}`",
                f"- descriptor_source_note: `{row['descriptor_source_note']}`",
                "",
            ]
        )

    lines.extend(["## Top Family Recommendations", ""])
    for state_id, group in recommendations.groupby("state_id", sort=False):
        lines.append(f"### {state_id}")
        lines.append("")
        for _, row in group.head(max_families_per_state).iterrows():
            lines.append(
                (
                    f"- rank {int(row['state_rank'])}: `{row['polymer_key']} + {row['cosolvent_system']}` "
                    f"(score `{row['recommendation_score']:.2f}`, classes `{row['matched_api_classes']}`)"
                )
            )
            lines.append(f"  rationale: {row['rationale']}")
        lines.append("")

    lines.extend(
        [
            "## Downstream Mechanistic Screen Matrix",
            "",
            f"- generated rows: `{len(candidate_matrix)}`",
            "- neutral and ionic states are both supported in the recommendation layer, but their mechanistic rankings should still be interpreted separately.",
            "- descriptor-level charge-state handling is metadata-based unless an explicit `state_smiles` is supplied for the charged form.",
            "",
        ]
    )
    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    payload = load_json(args.input_json)
    parent_smiles = resolve_parent_smiles(payload)
    state_options = normalize_state_options(payload, parent_smiles)
    salts = normalize_salt_options(payload)
    priors = load_prior_table(args.prior_table)
    polymer_library = load_polymer_library(args.polymer_library)
    max_families_per_state = int(args.max_families_per_state or payload.get("top_k_families", 4) or 4)

    state_descriptor_rows: list[dict[str, Any]] = []
    recommendation_frames: list[pd.DataFrame] = []

    for state in state_options:
        descriptor_row = compute_api_descriptors(parent_smiles, state)
        primary_class, api_classes = assign_api_classes(descriptor_row)
        descriptor_row["primary_api_class"] = primary_class
        descriptor_row["assigned_api_classes"] = ",".join(api_classes)
        state_descriptor_rows.append(descriptor_row)

        recommendation_frame = recommend_families(
            payload=payload,
            state=state,
            descriptors=descriptor_row,
            primary_class=primary_class,
            api_classes=api_classes,
            priors=priors,
            polymer_library=polymer_library,
        )
        if not recommendation_frame.empty:
            recommendation_frames.append(recommendation_frame)

    state_rows = pd.DataFrame(state_descriptor_rows)
    recommendations = (
        pd.concat(recommendation_frames, ignore_index=True)
        if recommendation_frames
        else pd.DataFrame(
            columns=[
                "state_id",
                "polymer_key",
                "cosolvent_system",
                "recommendation_score",
                "state_rank",
            ]
        )
    )
    candidate_matrix = build_candidate_matrix(
        payload=payload,
        parent_smiles=parent_smiles,
        state_options=state_options,
        salts=salts,
        recommendations=recommendations,
        max_families_per_state=max_families_per_state,
    )

    descriptor_json_path = args.output_dir / "descriptor_summary.json"
    state_csv_path = args.output_dir / "api_state_descriptors.csv"
    recommendation_csv_path = args.output_dir / "family_recommendations.csv"
    candidate_csv_path = args.output_dir / "mechanistic_screen_candidates.csv"
    summary_md_path = args.output_dir / "summary.md"

    descriptor_json_path.write_text(
        json.dumps(
            {
                "input_json": str(args.input_json),
                "api_name": payload.get("api_name"),
                "api_smiles": parent_smiles,
                "state_options": state_options,
                "salt_options": salts,
                "state_descriptors": state_rows.to_dict(orient="records"),
            },
            indent=2,
        )
        + "\n"
    )
    state_rows.to_csv(state_csv_path, index=False)
    recommendations.to_csv(recommendation_csv_path, index=False)
    candidate_matrix.to_csv(candidate_csv_path, index=False)
    write_summary(summary_md_path, payload, state_rows, recommendations, candidate_matrix, max_families_per_state)

    print(
        json.dumps(
            {
                "output_dir": str(args.output_dir),
                "state_csv": str(state_csv_path),
                "recommendation_csv": str(recommendation_csv_path),
                "candidate_csv": str(candidate_csv_path),
                "summary_md": str(summary_md_path),
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
