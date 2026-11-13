#!/usr/bin/env bash
#SBATCH --job-name=rff_ec1p2
#SBATCH --partition=g-3090-2
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/home/junxuanliu/projects/RNA-Folding-Frontier/logs/ec1_p2_%j.out
#
# RiNALMo masked pseudo-perplexity for E-C1.
set -euo pipefail

ROOT="/home/junxuanliu/projects/RNA-Folding-Frontier"
PY="/home/junxuanliu/.conda/envs/rff/bin/python"
cd "${ROOT}"

"${PY}" scripts/ec1_p2_pppl.py \
  --families experiments/002_e1prime_runs/ec1_families.json \
  --split splits/s1prime_filtered/test_withheld.jsonl \
  --checkpoint data/models/rinalmo-giga \
  --masks 4 --mask-frac 0.15 --seeds 17 23 41 53 \
  --out experiments/002_e1prime_runs/ec1_p2_pppl.json
