# Polymer Proxy Guide

These are the polymer fragments currently built into the screening workflow.

## Why These Are Proxies

The screen does not use full bulk polymer chains. It uses compact oligomer or motif-level fragments that preserve the local chemistry most relevant for:

- API-polymer contacts
- local solvation competition
- local aggregation suppression
- reactive hotspot screening

## Current Polymer Entries

### `peg`

- name: PEG/PEO trimer-like fragment
- family: polyether
- repeat units: 3
- end capping: dimethyl-capped oligomer ends
- reason: compact PEG-like ether backbone probe for local compatibility screening

### `pvp`

- name: PVP trimer-like fragment
- family: vinyl lactam polymer
- repeat units: 3
- end capping: methyl-capped oligomer ends
- reason: local lactam-rich environment for neutral API stabilization tests

### `eudragit_l100_like`

- name: Eudragit L100-like fragment
- family: methacrylate copolymer
- repeat units: 3
- end capping: simplified alkyl-capped copolymer fragment
- reason: acidic methacrylate-style proxy for pH-sensitive and ionic formulation behavior

### `gelatin`

- name: gelatin-like peptide fragment
- family: polypeptide
- repeat units: 2
- end capping: N-acetyl / N-methylamide caps
- reason: simple proteinaceous backbone proxy for hydrogen-bond-rich local environments

### `pectin`

- name: pectin-like digalacturonic acid dimer proxy
- family: polygalacturonate polysaccharide
- repeat units: 2
- end capping: native dimer termini
- reason: carbohydrate-rich, acid-containing proxy for gel-like polysaccharide chemistry

### `hpc`

- name: HPC-like hydroxypropyl cellobiose proxy
- family: cellulose ether
- repeat units: 2
- end capping: fixed hydroxypropyl-substituted cellobiose-like dimer
- reason: cellulose-ether proxy for semi-solid excipient screening

### `xanthan`

- name: xanthan-like branched glucan proxy
- family: branched polysaccharide
- repeat units: 3
- end capping: minimal branched saccharide proxy
- reason: branched polysaccharide environment for network-forming excipient comparisons

## Where The Files Are

Each polymer has:

- `SMILES_3D/polymer/<key>.smi`
- `SMILES_3D/polymer/<key>.xyz`
- `SMILES_3D/polymer/<key>.json`

The `.json` file includes the fragment id, family, repeat-unit count, end-capping strategy, and the original rationale.
