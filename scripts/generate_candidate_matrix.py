#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from run_api_family_recommendation import (
    build_candidate_matrix,
    load_polymer_library,
    load_prior_table,
    recommend_families,
)


def generate_family_recommendations(
    payload: dict[str, Any],
    parent_smiles: str,
    state_options: list[dict[str, Any]],
    salts: list[dict[str, Any]],
    annotated_descriptor_rows: list[dict[str, Any]],
    prior_table_path: Path,
    polymer_library_path: Path,
    max_families_per_state: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    priors = load_prior_table(prior_table_path)
    polymer_library = load_polymer_library(polymer_library_path)

    recommendation_frames: list[pd.DataFrame] = []
    state_lookup = {state["state_id"]: state for state in state_options}

    for descriptor_row in annotated_descriptor_rows:
        state = state_lookup[descriptor_row["state_id"]]
        recommendation_frame = recommend_families(
            payload=payload,
            state=state,
            descriptors=descriptor_row,
            primary_class=descriptor_row["primary_api_class"],
            api_classes=str(descriptor_row["assigned_api_classes"]).split(","),
            priors=priors,
            polymer_library=polymer_library,
        )
        if not recommendation_frame.empty:
            recommendation_frames.append(recommendation_frame)

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
    return recommendations, candidate_matrix, polymer_library


__all__ = [
    "build_candidate_matrix",
    "generate_family_recommendations",
    "load_polymer_library",
    "load_prior_table",
    "recommend_families",
]
