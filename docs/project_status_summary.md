# Project Status Summary

## Current Project State

This repository is now in a usable `pre-formulation recommendation + mechanistic screening` state.

Implemented layers:

- organized `SMILES_3D` library for APIs, polymer proxies, solvents, and co-solvent presets
- rule-based API family recommendation from `API SMILES`
- mechanistic cluster screening with FAIRChem-style checkpoints
- neutral versus ionic branch separation in mechanistic scoring
- replicate-ready screening workflow for stability-aware ranking
- archived benchmark outputs for both a relaxed GPU demo and a methotrexate family recommendation run

The project is no longer only a concept document. It now has runnable scripts, documented inputs, tracked chemistry assets, and archived outputs.

## What The System Can Do Now

### Layer 1: API family recommendation

From an input JSON containing API SMILES and charge-state options, the repo can now:

- compute formulation-relevant molecular descriptors
- assign formulation-relevant API classes
- apply a curated family prior table
- rank formulation families
- generate a downstream mechanistic screening candidate matrix

Main script:

- `scripts/run_api_family_recommendation.py`

### Layer 2: mechanistic screening

For shortlisted candidates, the repo can:

- build local `API + polymer + solvent/co-solvent` clusters
- relax them with a checkpoint-based FAIRChem workflow
- score local interaction and contact patterns
- separate neutral and ionic branches
- run multi-seed replicate screens for more stable ranking

Main script:

- `scripts/run_preformulation_mechanistic_screen.py`

## Current Archived Findings

### 1. Relaxed GPU demo

Archived here:

- `results/mechanistic_screen_relaxed_gpu_v3/`

Main finding:

- for the paracetamol demo, `PVP + water` remained the top neutral candidate after relaxation

Interpretation:

- the workflow is directionally useful for shortlist generation
- relaxation matters strongly
- static single-geometry rankings are not enough

### 2. Methotrexate family recommendation v1

Archived here:

- `results/methotrexate_family_recommendation_v1/`

Main descriptor-level finding:

- methotrexate currently looks like a highly polar, flexible, acid-bearing API rather than a classic hydrophobic aggregation-dominant case

Current descriptor picture:

- molecular weight: about `454.45`
- `cLogP`: about `0.27`
- `TPSA`: about `210.54`
- rotatable bonds: `9`
- water-affinity proxy: `0.85` to `1.00`
- self-aggregation risk proxy: `0.36` to `0.38`

Current family shortlist:

Neutral branch:

1. `PVP + water`
2. `HPC + water`
3. `HPC + ethanol/water`
4. `Gelatin + water`

Monoanionic branch:

1. `HPC + ethanol/water`
2. `PVP + water`
3. `Xanthan + water`
4. `HPC + water`

Interpretation:

- `PVP + water` is the best current neutral binder-first route
- `HPC + ethanol/water` is the best current ionic or process-oriented route
- `Xanthan + water` is a useful gel-network comparator, not the leading overall choice

## What We Can Honestly Claim Now

- the repository supports a real unseen-API family recommendation workflow
- methotrexate can now be fed through the workflow end-to-end
- the project can produce a rational shortlist for mechanistic screening and later DFT selection
- the current methotrexate shortlist is chemically plausible enough to guide the next computational step

## What We Should Not Claim Yet

- that the current system predicts final formulation success
- that the current methotrexate ranking is already experimentally validated
- that neutral and ionic candidates are directly cross-comparable in mechanistic scoring
- that metadata-only charge overrides are good enough for production ionic screening

## Current Scientific Limits

The main current limits are:

- ionic recommendation still uses metadata-level charge overrides unless explicit ionized `state_smiles` are supplied
- family recommendation is rule-based and prior-based, not trained on experimental outcome labels
- mechanistic ranking for methotrexate has not yet been run with replicate relaxation on the methotrexate shortlist
- DFT spot checks have not yet been done for methotrexate-specific clusters

## Best Next Steps

The next technically serious path is:

1. build explicit ionized methotrexate state SMILES instead of metadata-only charge overrides
2. generate a replicate mechanistic screen for the neutral methotrexate top three:
   - `PVP + water`
   - `HPC + water`
   - `HPC + ethanol/water`
3. then screen the top ionic branch:
   - `HPC + ethanol/water`
   - `PVP + water`
4. choose top, middle, and weak candidates for DFT spot checks

## Bottom Line

The project has moved from:

- a polymer-oriented idea plus a demo screen

to:

- a documented unseen-API recommendation workflow with runnable code, methotrexate chemistry assets, and a first methotrexate shortlist

That is a real project milestone, but it is still a shortlist-generation system, not a final formulation predictor.
