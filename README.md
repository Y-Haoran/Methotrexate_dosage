# Methotrexate Dosage

Mechanistic preformulation screening workspace for semi-solid and polymer-assisted drug formulation design.

This repository now contains a clean standalone slice of the atomistic screening workflow built around FAIRChem-style molecular potentials, polymer fragment proxies, candidate formulation matrices, and archived benchmark outputs. The immediate purpose is to support fast ranking of API-polymer-cosolvent candidates before DFT spot checks and wet-lab validation.

## What This Repo Does

- screens local `API + polymer + solvent/cosolvent` clusters
- compares API-polymer association against API self-aggregation
- separates neutral and ionic branches so unstable ionic references are not over-interpreted
- supports multi-seed replicate screening with `mean ± sd`, rank stability, and top-1/top-2 frequencies
- keeps a compact archive of the first relaxed GPU demo and the first DFT shortlist

## Current Scope

The screening layer is useful for:

- local compatibility
- local solvation competition
- early excipient ranking
- shortlist generation for DFT

It is not a direct predictor of:

- bulk rheology
- final printability windows
- dissolution curves
- long-term storage success

## Repository Layout

```text
.
├── data/
│   └── preformulation/
│       ├── paracetamol_mechanistic_screen_demo.csv
│       ├── candidate_matrix_template.csv
│       ├── candidate_matrix_augmented_template.csv
│       ├── polymer_descriptor_library.csv
│       ├── polymer_family_aliases.csv
│       └── screening_results_template.csv
├── docs/
│   └── dft_spotcheck_plan_v1.md
├── jobs/
│   └── run_preformulation_mechanistic_screen_gpu_normal.lsf.sh
├── results/
│   └── mechanistic_screen_relaxed_gpu_v3/
├── scripts/
│   ├── run_formulation_descriptor_pilot.py
│   ├── run_preformulation_mechanistic_screen.py
│   └── utils/sql_queries/
└── checkpoints/
```

`scripts/utils/sql_queries/` contains the original SQL assets already present in the upstream repository and has been left intact.

## Main Workflow

1. Put a compatible FAIRChem checkpoint into `checkpoints/best_inference_ckpt.pt`.
2. Start from `data/preformulation/candidate_matrix_template.csv` or replace the demo matrix with a methotrexate-specific matrix.
3. Run the screening script locally or through the LSF wrapper.
4. Use the aggregated neutral ranking to choose DFT spot-check candidates.
5. Validate the final shortlist experimentally.

## Quick Start

### Local run

```bash
python scripts/run_preformulation_mechanistic_screen.py \
  --candidate-matrix data/preformulation/paracetamol_mechanistic_screen_demo.csv \
  --checkpoint-path checkpoints/best_inference_ckpt.pt \
  --device cpu \
  --relax-steps 20 \
  --n-replicates 5 \
  --seed 42 \
  --ranking-group neutral \
  --output-dir results/mechanistic_screen_run
```

### GPU batch run

```bash
bsub < jobs/run_preformulation_mechanistic_screen_gpu_normal.lsf.sh
```

Environment variables accepted by the batch wrapper:

- `CANDIDATE_MATRIX`
- `CHECKPOINT_PATH`
- `OUTPUT_DIR`
- `RELAX_STEPS`
- `N_REPLICATES`
- `SEED`
- `RANKING_GROUP`
- `VENV`

## Key Files

- `scripts/run_preformulation_mechanistic_screen.py`: main mechanistic screening runner
- `scripts/run_formulation_descriptor_pilot.py`: polymer fragment/descriptor helper definitions used by the screen
- `docs/dft_spotcheck_plan_v1.md`: first DFT shortlist plan
- `results/mechanistic_screen_relaxed_gpu_v3/summary.md`: archived relaxed GPU demo summary

## Current Archived Result

The repo includes the first relaxed GPU demo archive in `results/mechanistic_screen_relaxed_gpu_v3/`.

High-level takeaway from that archive:

- top neutral demo candidate: `paracetamol + pvp + water`
- neutral and ionic rankings are separated
- the ranking is useful for shortlist generation, not final decision-making
- the next strong scientific step is replicate-based neutral screening followed by DFT checks

## Adapting This Repo To Methotrexate

The included demo matrix is still a paracetamol placeholder. For methotrexate-focused work, the first changes should be:

- replace the demo matrix with methotrexate plus relevant protonation or salt states
- add the excipient families that match the real formulation program
- run neutral replicate screening first
- use the aggregated ranking to define the DFT batch

## Notes

- `checkpoints/` is intentionally empty in git; place model weights there locally.
- `logs/` and new run outputs are ignored by default.
- If you use a temporary GitHub token to push this repo, revoke and rotate it afterward.
