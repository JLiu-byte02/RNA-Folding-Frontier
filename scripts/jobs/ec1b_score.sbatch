#!/usr/bin/env bash
#SBATCH --job-name=rff_ec1b
#SBATCH --partition=g-3090-2
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/home/junxuanliu/projects/RNA-Folding-Frontier/logs/ec1b_%j.out
#
# E-C1b decomposition and one-hot specificity control.
set -euo pipefail

ROOT="/home/junxuanliu/projects/RNA-Folding-Frontier"
PY="/home/junxuanliu/.conda/envs/rff/bin/python"
cd "${ROOT}"

"${PY}" scripts/ec1b_score.py \
  --ec1 experiments/002_e1prime_runs/ec1_p1_p2.json \
  --seeds 17 23 41 \
  --out experiments/002_e1prime_runs/ec1b_specificity.json \
  --md-out docs/ec1b_2026-11-13.md
