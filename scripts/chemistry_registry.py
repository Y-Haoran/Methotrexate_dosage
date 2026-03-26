from __future__ import annotations

from dataclasses import dataclass


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
