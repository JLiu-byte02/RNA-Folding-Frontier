#!/usr/bin/env bash
#SBATCH --job-name=rff_ec1
#SBATCH --partition=g-3090-2
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/home/junxuanliu/projects/RNA-Folding-Frontier/logs/ec1_score_%j.out
#
# CPU-only E-C1 outcome and partial correlations.
set -euo pipefail

ROOT="/home/junxuanliu/projects/RNA-Folding-Frontier"
PY="/home/junxuanliu/.conda/envs/rff/bin/python"
cd "${ROOT}"

"${PY}" scripts/ec1_score.py \
  --families experiments/002_e1prime_runs/ec1_families.json \
  --split splits/s1prime_filtered/test_withheld.jsonl \
  --fm-dump experiments/002_e1prime_runs/dump_p2fix_seed17/a_test_withheld.jsonl \
  --rfam-family data/raw/rfam/family.txt \
  --p2 experiments/002_e1prime_runs/ec1_p2_pppl.json \
  --out experiments/002_e1prime_runs/ec1_p1_p2.json \
  --md-out docs/ec1_2026-11-13.md
