#!/usr/bin/env bash
#SBATCH --job-name=rff_ec1fam
#SBATCH --partition=g-3090-2
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/home/junxuanliu/projects/RNA-Folding-Frontier/logs/ec1_families_%j.out
#
# Build the frozen E-C1 family list from the P3 dumps.
set -euo pipefail

ROOT="/home/junxuanliu/projects/RNA-Folding-Frontier"
PY="/home/junxuanliu/.conda/envs/rff/bin/python"
cd "${ROOT}"

"${PY}" scripts/build_ec1_family_list.py \
  --dump-root experiments/002_e1prime_runs \
  --split test_withheld --seeds 17 23 41 --min-true-pairs 4 \
  --out experiments/002_e1prime_runs/ec1_families.json
