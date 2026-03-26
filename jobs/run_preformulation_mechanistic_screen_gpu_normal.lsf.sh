#!/usr/bin/env bash
#BSUB -J mtx_preform_screen
#BSUB -q gpu-normal
#BSUB -gpu "num=1:mode=shared"
#BSUB -n 1
#BSUB -M 8000
#BSUB -R "select[mem>8000] span[hosts=1] rusage[mem=8000]"
#BSUB -W 04:00
#BSUB -oo logs/mtx_preform_screen.%J.out
#BSUB -eo logs/mtx_preform_screen.%J.err

set -euo pipefail

ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
VENV="${VENV:-$ROOT_DIR/.venv}"
CANDIDATE_MATRIX="${CANDIDATE_MATRIX:-$ROOT_DIR/data/preformulation/paracetamol_mechanistic_screen_demo.csv}"
OUTPUT_DIR="${OUTPUT_DIR:-$ROOT_DIR/results/mechanistic_screen_gpu_run}"
CHECKPOINT_PATH="${CHECKPOINT_PATH:-$ROOT_DIR/checkpoints/best_inference_ckpt.pt}"
RELAX_STEPS="${RELAX_STEPS:-20}"
N_REPLICATES="${N_REPLICATES:-5}"
SEED="${SEED:-42}"
RANKING_GROUP="${RANKING_GROUP:-neutral}"

mkdir -p "$ROOT_DIR/logs" "$OUTPUT_DIR"
source "$VENV/bin/activate"

python "$ROOT_DIR/scripts/run_preformulation_mechanistic_screen.py" \
  --candidate-matrix "$CANDIDATE_MATRIX" \
  --checkpoint-path "$CHECKPOINT_PATH" \
  --task-name omol \
  --device cuda \
  --relax-steps "$RELAX_STEPS" \
  --n-replicates "$N_REPLICATES" \
  --seed "$SEED" \
  --ranking-group "$RANKING_GROUP" \
  --output-dir "$OUTPUT_DIR"
