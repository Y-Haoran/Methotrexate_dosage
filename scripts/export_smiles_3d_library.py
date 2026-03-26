#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

from ase.io import write
from chemistry_registry import API_LIBRARY, COSOLVENT_PRESETS, SMALL_MOLECULES, MoleculeSpec
from run_formulation_descriptor_pilot import DEFAULT_BASE_DIR, POLYMER_SPECS, build_rdkit_molecule, rdkit_to_ase


OUTPUT_DIR = DEFAULT_BASE_DIR / "SMILES_3D"


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def build_atoms(spec: MoleculeSpec):
    mol, _ = build_rdkit_molecule(spec.smiles)
    return rdkit_to_ase(mol, charge=spec.charge, spin=spec.spin)


def write_smiles_and_3d(path: Path, key: str, spec: MoleculeSpec, metadata: dict[str, object]) -> dict[str, object]:
    ensure_dir(path)
    smiles_path = path / f"{key}.smi"
    xyz_path = path / f"{key}.xyz"
    json_path = path / f"{key}.json"

    atoms = build_atoms(spec)
    smiles_path.write_text(spec.smiles + "\n")
    write(str(xyz_path), atoms)
    json_path.write_text(json.dumps(metadata, indent=2) + "\n")

    return {
        "key": key,
        "name": spec.name,
        "smiles_file": str(smiles_path.relative_to(DEFAULT_BASE_DIR)),
        "xyz_file": str(xyz_path.relative_to(DEFAULT_BASE_DIR)),
        "json_file": str(json_path.relative_to(DEFAULT_BASE_DIR)),
    }


def export_api_library() -> list[dict[str, object]]:
    entries: list[dict[str, object]] = []
    for key, spec in sorted(API_LIBRARY.items()):
        entries.append(
            write_smiles_and_3d(
                OUTPUT_DIR / "api",
                key,
                spec,
                {
                    "key": key,
                    "category": "api",
                    "name": spec.name,
                    "smiles": spec.smiles,
                    "charge": spec.charge,
                    "spin": spec.spin,
                },
            )
        )
    return entries


def export_polymer_library() -> list[dict[str, object]]:
    entries: list[dict[str, object]] = []
    for key, spec in sorted(POLYMER_SPECS.items()):
        entries.append(
            write_smiles_and_3d(
                OUTPUT_DIR / "polymer",
                key,
                MoleculeSpec(spec.name, spec.smiles, spec.charge, spec.spin),
                {
                    "key": key,
                    "category": "polymer",
                    "name": spec.name,
                    "family": spec.family,
                    "smiles": spec.smiles,
                    "charge": spec.charge,
                    "spin": spec.spin,
                    "fragment_id": spec.fragment_id,
                    "repeat_units": spec.repeat_units,
                    "end_capping_strategy": spec.end_capping_strategy,
                    "probe_rationale": spec.probe_rationale,
                    "note": spec.note,
                },
            )
        )
    return entries


def export_solvent_library() -> list[dict[str, object]]:
    entries: list[dict[str, object]] = []
    solvent_keys = ["water", "ethanol", "glycerol"]
    for key in solvent_keys:
        spec = SMALL_MOLECULES[key]
        entries.append(
            write_smiles_and_3d(
                OUTPUT_DIR / "solvent",
                key,
                spec,
                {
                    "key": key,
                    "category": "solvent",
                    "name": spec.name,
                    "smiles": spec.smiles,
                    "charge": spec.charge,
                    "spin": spec.spin,
                },
            )
        )
    return entries


def export_cosolvent_presets() -> list[dict[str, object]]:
    ensure_dir(OUTPUT_DIR / "co_solvent")
    entries: list[dict[str, object]] = []
    for key, components in sorted(COSOLVENT_PRESETS.items()):
        payload = {
            "key": key,
            "category": "co_solvent",
            "kind": "mixture_preset",
            "components": [
                {
                    "molecule_key": molecule_key,
                    "name": SMALL_MOLECULES[molecule_key].name,
                    "smiles": SMALL_MOLECULES[molecule_key].smiles,
                    "count": count,
                }
                for molecule_key, count in sorted(components.items())
            ],
            "note": "Co-solvent presets are stored as component manifests rather than one combined XYZ structure.",
        }
        json_path = OUTPUT_DIR / "co_solvent" / f"{key}.json"
        json_path.write_text(json.dumps(payload, indent=2) + "\n")
        entries.append(
            {
                "key": key,
                "json_file": str(json_path.relative_to(DEFAULT_BASE_DIR)),
            }
        )
    return entries


def write_manifest(payload: dict[str, object]) -> None:
    manifest_path = OUTPUT_DIR / "library_manifest.json"
    manifest_path.write_text(json.dumps(payload, indent=2) + "\n")


def main() -> int:
    ensure_dir(OUTPUT_DIR)
    payload = {
        "output_dir": str(OUTPUT_DIR.relative_to(DEFAULT_BASE_DIR)),
        "api": export_api_library(),
        "polymer": export_polymer_library(),
        "solvent": export_solvent_library(),
        "co_solvent": export_cosolvent_presets(),
    }
    write_manifest(payload)
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
