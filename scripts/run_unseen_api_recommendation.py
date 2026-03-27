#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path
from typing import Any

import pandas as pd

from classify_api import classify_state_descriptor_rows
from explain_results import write_user_facing_outputs
from extract_api_descriptors import build_default_salt_options, build_state_descriptor_rows, load_json
from generate_candidate_matrix import generate_family_recommendations


DEFAULT_BASE_DIR = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "User-facing entry point for unseen-API formulation-family recommendation. "
            "Provide an API, and the platform will generate candidate systems, rank them, "
            "and optionally launch mechanistic screening if a checkpoint is available."
        )
    )
    parser.add_argument("--input-json", type=Path, default=None)
    parser.add_argument("--api-name", default="")
    parser.add_argument("--api-smiles", default="")
    parser.add_argument("--context", default="semi_solid_printing")
    parser.add_argument("--route", default="topical")
    parser.add_argument("--processing-context", default="drying_after_extrusion")
    parser.add_argument("--ph-context", default="unknown")
    parser.add_argument("--charge-states", default="auto", help="Comma-separated, for example `0,-1,-2`, or `auto`.")
    parser.add_argument("--top-k", type=int, default=4)
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
    parser.add_argument(
        "--checkpoint-path",
        type=Path,
        default=DEFAULT_BASE_DIR / "checkpoints" / "best_inference_ckpt.pt",
    )
    parser.add_argument("--device", choices=["auto", "cpu", "cuda"], default="auto")
    parser.add_argument("--relax-steps", type=int, default=20)
    parser.add_argument("--n-replicates", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--skip-mechanistic-screen", action="store_true")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_BASE_DIR / "results" / "unseen_api_recommendation_run",
    )
    return parser.parse_args()


def state_entry_for_charge(formal_charge: int, include_in_screen: bool) -> dict[str, Any]:
    if formal_charge == 0:
        state_id, label = "neutral", "neutral"
    elif formal_charge == -1:
        state_id, label = "monoanion", "monoanionic"
    elif formal_charge == -2:
        state_id, label = "dianion", "dianionic"
    elif formal_charge == 1:
        state_id, label = "monocation", "monocationic"
    else:
        state_id = f"charge_{formal_charge:+d}".replace("+", "plus_").replace("-", "minus_")
        label = f"charge {formal_charge:+d}"
    return {
        "state_id": state_id,
        "label": label,
        "formal_charge": formal_charge,
        "api_spin": 1,
        "state_smiles": None,
        "include_in_screen": include_in_screen,
        "screen_priority": "primary" if include_in_screen else "secondary",
    }


def parse_charge_states(charge_states: str) -> list[int] | None:
    if not charge_states or charge_states.strip().lower() == "auto":
        return None
    return [int(token.strip()) for token in charge_states.split(",") if token.strip()]


def build_payload_from_args(args: argparse.Namespace) -> dict[str, Any]:
    if args.input_json is not None:
        return load_json(args.input_json)

    if not args.api_name:
        raise ValueError("--api-name is required when --input-json is not provided.")

    payload: dict[str, Any] = {
        "api_name": args.api_name,
        "api_smiles": args.api_smiles,
        "formulation_context": args.context,
        "route_of_administration": args.route,
        "processing_context": args.processing_context,
        "pH_context": args.ph_context,
        "top_k_families": args.top_k,
    }

    parsed_charge_states = parse_charge_states(args.charge_states)
    if parsed_charge_states is not None:
        state_options = []
        for index, charge in enumerate(parsed_charge_states):
            include_in_screen = index < max(args.top_k, 1)
            if charge == -2:
                include_in_screen = False
            state_options.append(state_entry_for_charge(charge, include_in_screen))
        payload["state_options"] = state_options
        payload["salt_options"] = build_default_salt_options(state_options)
    return payload


def write_internal_outputs(
    output_dir: Path,
    payload: dict[str, Any],
    state_rows: pd.DataFrame,
    recommendations: pd.DataFrame,
    candidate_matrix: pd.DataFrame,
    parent_smiles: str,
    state_options: list[dict[str, Any]],
    salts: list[dict[str, Any]],
) -> None:
    internal_dir = output_dir / "internal"
    internal_dir.mkdir(parents=True, exist_ok=True)
    (internal_dir / "input_payload.json").write_text(json.dumps(payload, indent=2) + "\n")
    (internal_dir / "descriptor_summary.json").write_text(
        json.dumps(
            {
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
    state_rows.to_csv(internal_dir / "api_state_descriptors.csv", index=False)
    recommendations.to_csv(internal_dir / "family_recommendations_full.csv", index=False)
    candidate_matrix.to_csv(internal_dir / "mechanistic_screen_candidates_internal.csv", index=False)


def maybe_run_mechanistic_screen(
    args: argparse.Namespace,
    candidate_matrix_path: Path,
    output_dir: Path,
) -> str:
    if args.skip_mechanistic_screen:
        return "skipped (user requested --skip-mechanistic-screen)"
    if not args.checkpoint_path.exists():
        return f"skipped (checkpoint not found at {args.checkpoint_path})"

    mechanistic_output_dir = output_dir / "mechanistic_screen"
    cmd = [
        sys.executable,
        str(DEFAULT_BASE_DIR / "scripts" / "run_preformulation_mechanistic_screen.py"),
        "--candidate-matrix",
        str(candidate_matrix_path),
        "--checkpoint-path",
        str(args.checkpoint_path),
        "--device",
        args.device,
        "--relax-steps",
        str(args.relax_steps),
        "--n-replicates",
        str(args.n_replicates),
        "--seed",
        str(args.seed),
        "--ranking-group",
        "all",
        "--output-dir",
        str(mechanistic_output_dir),
    ]
    subprocess.run(cmd, check=True)
    return f"completed ({mechanistic_output_dir})"


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    payload = build_payload_from_args(args)
    parent_smiles, state_options, salts, descriptor_rows = build_state_descriptor_rows(payload)
    annotated_descriptor_rows = classify_state_descriptor_rows(descriptor_rows)
    state_rows = pd.DataFrame(annotated_descriptor_rows)

    recommendations, candidate_matrix, _ = generate_family_recommendations(
        payload=payload,
        parent_smiles=parent_smiles,
        state_options=state_options,
        salts=salts,
        annotated_descriptor_rows=annotated_descriptor_rows,
        prior_table_path=args.prior_table,
        polymer_library_path=args.polymer_library,
        max_families_per_state=args.top_k,
    )

    write_internal_outputs(
        output_dir=args.output_dir,
        payload=payload,
        state_rows=state_rows,
        recommendations=recommendations,
        candidate_matrix=candidate_matrix,
        parent_smiles=parent_smiles,
        state_options=state_options,
        salts=salts,
    )

    candidate_matrix_path = args.output_dir / "candidate_matrix_generated.csv"
    candidate_matrix.to_csv(candidate_matrix_path, index=False)
    mechanistic_status = maybe_run_mechanistic_screen(args, candidate_matrix_path, args.output_dir)
    write_user_facing_outputs(
        output_dir=args.output_dir,
        payload=payload,
        recommendations=recommendations,
        candidate_matrix=candidate_matrix,
        mechanistic_status=mechanistic_status,
    )

    print(
        json.dumps(
            {
                "output_dir": str(args.output_dir),
                "summary": str(args.output_dir / "summary.md"),
                "ranking": str(args.output_dir / "ranking.csv"),
                "candidate_matrix": str(candidate_matrix_path),
                "buy_list": str(args.output_dir / "buy_list.txt"),
                "lab_test_plan": str(args.output_dir / "lab_test_plan.txt"),
                "mechanistic_status": mechanistic_status,
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
