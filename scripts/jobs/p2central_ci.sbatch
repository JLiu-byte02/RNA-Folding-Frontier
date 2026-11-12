#!/usr/bin/env bash
#SBATCH --job-name=rff_p2ci
#SBATCH --partition=g-3090-2
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/home/junxuanliu/projects/RNA-Folding-Frontier/logs/p2ci_%j.out
#
# CPU-only bootstrap CIs for the corrected P2-central retention ratios.
set -euo pipefail

ROOT="/home/junxuanliu/projects/RNA-Folding-Frontier"
PY="/home/junxuanliu/.conda/envs/rff/bin/python"
cd "${ROOT}"

"${PY}" scripts/bootstrap_p2central.py \
  --out experiments/002_e1prime_runs/p2central_ratios_ci.json \
  --md-out docs/p2central_ratios_ci_2026-11-12.md
