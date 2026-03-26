# SMILES_3D Library

This folder stores the small chemical building blocks used in the mechanistic preformulation screen.

## Structure

- `api/`: API entries with one `.smi`, one `.xyz`, and one `.json` metadata file per molecule
- `polymer/`: polymer proxy fragments with one `.smi`, one `.xyz`, and one `.json` metadata file per fragment
- `solvent/`: single-molecule solvent entries with one `.smi`, one `.xyz`, and one `.json` metadata file per solvent
- `co_solvent/`: mixture presets stored as `.json` manifests because a co-solvent system is not one single molecule
- `library_manifest.json`: machine-readable index of all exported entries

## Notes

- The polymer entries are proxy fragments, not full commercial polymer chains.
- The `.xyz` files are RDKit-built 3D starting structures.
- The co-solvent presets point back to component solvents and counts rather than pretending the mixture is one molecule.
- Regenerate this library with `python scripts/export_smiles_3d_library.py`.
