#!/usr/bin/env python3
from __future__ import annotations

from typing import Any

from rdkit import Chem

from run_api_family_recommendation import (
    ACIDIC_PATTERNS,
    BASIC_PATTERNS,
    compute_api_descriptors,
    count_matches,
    load_json,
    normalize_salt_options,
    normalize_state_options,
    resolve_parent_smiles,
)


def infer_state_options_from_parent_smiles(parent_smiles: str) -> list[dict[str, Any]]:
    mol = Chem.MolFromSmiles(parent_smiles)
    if mol is None:
        raise ValueError("Invalid parent_smiles passed into infer_state_options_from_parent_smiles.")

    acidic_group_count = count_matches(mol, ACIDIC_PATTERNS)
    basic_group_count = count_matches(mol, BASIC_PATTERNS)

    state_options: list[dict[str, Any]] = [
        {
            "state_id": "neutral",
            "label": "neutral",
            "formal_charge": 0,
            "api_spin": 1,
            "state_smiles": None,
            "include_in_screen": True,
            "screen_priority": "primary",
        }
    ]

    if acidic_group_count > 0:
        state_options.append(
            {
                "state_id": "monoanion",
                "label": "monoanionic",
                "formal_charge": -1,
                "api_spin": 1,
                "state_smiles": None,
                "include_in_screen": True,
                "screen_priority": "primary",
            }
        )
        if acidic_group_count > 1:
            state_options.append(
                {
                    "state_id": "dianion",
                    "label": "dianionic",
                    "formal_charge": -2,
                    "api_spin": 1,
                    "state_smiles": None,
                    "include_in_screen": False,
                    "screen_priority": "secondary",
                }
            )
    elif basic_group_count > 0:
        state_options.append(
            {
                "state_id": "monocation",
                "label": "monocationic",
                "formal_charge": 1,
                "api_spin": 1,
                "state_smiles": None,
                "include_in_screen": True,
                "screen_priority": "primary",
            }
        )

    return state_options


def build_default_salt_options(state_options: list[dict[str, Any]]) -> list[dict[str, Any]]:
    salts: list[dict[str, Any]] = []
    for state in state_options:
        formal_charge = int(state["formal_charge"])
        if formal_charge < 0:
            salts.append(
                {
                    "salt_id": f"sodium_{state['state_id']}",
                    "state_id": state["state_id"],
                    "ion_name": "sodium",
                    "ion_smiles": "[Na+]",
                    "ion_charge": 1,
                    "ion_spin": 1,
                    "ion_count": abs(formal_charge),
                    "include_in_screen": bool(state.get("include_in_screen", True)),
                }
            )
    return salts


def prepare_payload_for_descriptor_extraction(payload: dict[str, Any]) -> tuple[str, list[dict[str, Any]], list[dict[str, Any]]]:
    parent_smiles = resolve_parent_smiles(payload)
    normalized_payload = dict(payload)
    if not normalized_payload.get("state_options"):
        normalized_payload["state_options"] = infer_state_options_from_parent_smiles(parent_smiles)
    if normalized_payload.get("salt_options") is None:
        normalized_payload["salt_options"] = build_default_salt_options(normalized_payload["state_options"])

    state_options = normalize_state_options(normalized_payload, parent_smiles)
    salts = normalize_salt_options(normalized_payload)
    return parent_smiles, state_options, salts


def build_state_descriptor_rows(payload: dict[str, Any]) -> tuple[str, list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    parent_smiles, state_options, salts = prepare_payload_for_descriptor_extraction(payload)
    descriptor_rows = [compute_api_descriptors(parent_smiles, state) for state in state_options]
    return parent_smiles, state_options, salts, descriptor_rows


__all__ = [
    "build_default_salt_options",
    "build_state_descriptor_rows",
    "compute_api_descriptors",
    "infer_state_options_from_parent_smiles",
    "load_json",
    "normalize_salt_options",
    "normalize_state_options",
    "prepare_payload_for_descriptor_extraction",
    "resolve_parent_smiles",
]
