# DFT Spot-Check Plan v1

Source ranking:
- `results/mechanistic_screen_relaxed_gpu_v3/mechanistic_screen_results_neutral.csv`

Scope:
- validate the relaxed neutral screen before using it for stronger formulation claims
- compare MLIP vs DFT on relative ordering, retained contact motifs, and relaxation trends

Recommended neutral candidates:

1. Top candidate
- `screen_001`
- `paracetamol + pvp + water`
- best relaxed neutral score in the current screen

2. Second candidate
- `screen_003`
- `paracetamol + hpc + ethanol_water`
- next-best relaxed neutral score and chemically distinct from PVP

3. Mid candidate
- `screen_002`
- `paracetamol + peg + water`
- intermediate relaxed score and a simple polyether environment

4. Weak candidate
- `screen_004`
- `paracetamol + pectin + glycerol_water`
- weakest relaxed API-containing neutral candidate in the current shortlist

Deferred ionic check:
- `screen_007`
- `diclofenac_anion + pvp + water + sodium`
- keep out of the first DFT batch until the ionic reference treatment is improved

For each selected candidate:
- single-point DFT on the initial cluster
- geometry relaxation at the same charge/spin state used in the MLIP screen
- compare relative energy ordering against MLIP
- inspect whether API stays associated with the polymer after relaxation
- record the dominant contact motifs before and after relaxation

Acceptance criterion:
- the MLIP top/mid/weak ordering is broadly preserved
- the main API-polymer contact motif remains chemically consistent after DFT relaxation
- no candidate that looked favorable in MLIP collapses into an obviously non-interacting geometry under DFT
