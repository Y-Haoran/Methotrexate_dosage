#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd


def build_user_ranking_frame(api_name: str, recommendations: pd.DataFrame) -> pd.DataFrame:
    if recommendations.empty:
        return pd.DataFrame(
            columns=[
                "api_name",
                "state_id",
                "rank",
                "polymer_key",
                "polymer_name",
                "cosolvent_system",
                "recommendation_score",
                "matched_api_classes",
                "rationale",
            ]
        )

    ranking = recommendations[
        [
            "state_id",
            "state_rank",
            "polymer_key",
            "polymer_name",
            "cosolvent_system",
            "recommendation_score",
            "matched_api_classes",
            "rationale",
        ]
    ].copy()
    ranking.insert(0, "api_name", api_name)
    ranking = ranking.rename(columns={"state_rank": "rank"})
    return ranking


def write_buy_list(path: Path, ranking: pd.DataFrame, candidate_matrix: pd.DataFrame) -> None:
    polymer_entries = sorted(
        {
            f"{row['polymer_key']}: {row['polymer_name']}"
            for _, row in ranking[["polymer_key", "polymer_name"]].drop_duplicates().iterrows()
        }
    )
    cosolvent_entries = sorted(set(candidate_matrix["cosolvent_system"].fillna("").astype(str)) - {""})
    ion_entries = sorted(set(candidate_matrix["ion_name"].fillna("").astype(str)) - {""})

    lines = ["Buy List", "", "Polymer Families To Source First:"]
    lines.extend(polymer_entries or ["none"])
    lines.extend(["", "Solvent / Co-Solvent Systems To Prepare:"])
    lines.extend(cosolvent_entries or ["none"])
    lines.extend(["", "Counterions Or Salts To Prepare:"])
    lines.extend(ion_entries or ["none"])
    path.write_text("\n".join(lines) + "\n")


def write_lab_test_plan(path: Path, api_name: str, ranking: pd.DataFrame) -> None:
    neutral_rows = ranking[ranking["state_id"] == "neutral"].head(2)
    ionic_rows = ranking[ranking["state_id"] != "neutral"].head(2)

    lines = [
        f"Lab Test Plan: {api_name}",
        "",
        "First Pass Candidates:",
    ]
    for _, row in neutral_rows.iterrows():
        lines.append(f"- neutral: {row['polymer_key']} + {row['cosolvent_system']}")
    for _, row in ionic_rows.iterrows():
        lines.append(f"- ionic: {row['polymer_key']} + {row['cosolvent_system']}")

    lines.extend(
        [
            "",
            "Recommended First Assays:",
            "- DSC / PXRD / FTIR on the top neutral candidates",
            "- rheology and print-focused handling on the top cellulose / gel-network routes",
            "- drying observation for ethanol-water candidates",
            "- DFT spot checks on one top, one middle, and one weak candidate before stronger claims",
        ]
    )
    path.write_text("\n".join(lines) + "\n")


def write_platform_summary(
    path: Path,
    payload: dict[str, Any],
    ranking: pd.DataFrame,
    candidate_matrix: pd.DataFrame,
    mechanistic_status: str,
) -> None:
    top_lines: list[str] = []
    for state_id, group in ranking.groupby("state_id", sort=False):
        top = group.head(3)
        if top.empty:
            continue
        top_lines.append(f"## {state_id}")
        top_lines.append("")
        for _, row in top.iterrows():
            top_lines.append(
                f"- rank {int(row['rank'])}: `{row['polymer_key']} + {row['cosolvent_system']}` "
                f"(score `{row['recommendation_score']:.2f}`)"
            )
            top_lines.append(f"  explanation: {row['rationale']}")
        top_lines.append("")

    lines = [
        f"# Unseen API Recommendation Summary: {payload.get('api_name', 'API')}",
        "",
        "Give us an API, and this platform recommends the most promising formulation families to test first.",
        "",
        "## Input",
        "",
        f"- api_name: `{payload.get('api_name', 'unknown')}`",
        f"- formulation_context: `{payload.get('formulation_context', 'unknown')}`",
        f"- route_of_administration: `{payload.get('route_of_administration', 'unknown')}`",
        f"- processing_context: `{payload.get('processing_context', 'unknown')}`",
        "",
        "## What The Platform Returned",
        "",
        f"- ranked family rows: `{len(ranking)}`",
        f"- generated candidate systems: `{len(candidate_matrix)}`",
        f"- mechanistic_screen_status: `{mechanistic_status}`",
        "",
        "## Recommended Formulation Families",
        "",
    ]
    lines.extend(top_lines or ["No recommendation rows were produced.", ""])
    lines.extend(
        [
            "## Output Files",
            "",
            "- `summary.md`",
            "- `ranking.csv`",
            "- `candidate_matrix_generated.csv`",
            "- `buy_list.txt`",
            "- `lab_test_plan.txt`",
            "",
            "## Limits",
            "",
            "- this platform recommends and ranks candidate formulation families",
            "- it does not predict final formulation ratios",
            "- it does not predict exact printability, dissolution, or long-term stability",
        ]
    )
    path.write_text("\n".join(lines) + "\n")


def write_user_facing_outputs(
    output_dir: Path,
    payload: dict[str, Any],
    recommendations: pd.DataFrame,
    candidate_matrix: pd.DataFrame,
    mechanistic_status: str,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    ranking = build_user_ranking_frame(str(payload.get("api_name", "api")), recommendations)
    ranking.to_csv(output_dir / "ranking.csv", index=False)
    candidate_matrix.to_csv(output_dir / "candidate_matrix_generated.csv", index=False)
    write_buy_list(output_dir / "buy_list.txt", ranking, candidate_matrix)
    write_lab_test_plan(output_dir / "lab_test_plan.txt", str(payload.get("api_name", "api")), ranking)
    write_platform_summary(output_dir / "summary.md", payload, ranking, candidate_matrix, mechanistic_status)


__all__ = [
    "build_user_ranking_frame",
    "write_buy_list",
    "write_lab_test_plan",
    "write_platform_summary",
    "write_user_facing_outputs",
]
