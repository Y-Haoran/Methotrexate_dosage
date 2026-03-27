#!/usr/bin/env python3
from __future__ import annotations

from typing import Any

from run_api_family_recommendation import assign_api_classes


def classify_state_descriptor_rows(descriptor_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    annotated_rows: list[dict[str, Any]] = []
    for descriptor_row in descriptor_rows:
        primary_class, api_classes = assign_api_classes(descriptor_row)
        annotated = dict(descriptor_row)
        annotated["primary_api_class"] = primary_class
        annotated["assigned_api_classes"] = ",".join(api_classes)
        annotated_rows.append(annotated)
    return annotated_rows


__all__ = ["assign_api_classes", "classify_state_descriptor_rows"]
