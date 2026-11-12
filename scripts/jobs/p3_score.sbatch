#!/usr/bin/env bash
#SBATCH --job-name=rff_p3score
#SBATCH --partition=g-3090-2
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/home/junxuanliu/projects/RNA-Folding-Frontier/logs/p3score_%j.out
#
# CPU-only D1-core table and verdict on withheld.
set -euo pipefail

ROOT="/home/junxuanliu/projects/RNA-Folding-Frontier"
PY="/home/junxuanliu/.conda/envs/rff/bin/python"
cd "${ROOT}"

"${PY}" scripts/score_p3_d1core.py \
  --split test_withheld \
  --nested-report \
  --out experiments/002_e1prime_runs/d1core_test_withheld.json \
  --md-out docs/d1core_2026-11-12.md
