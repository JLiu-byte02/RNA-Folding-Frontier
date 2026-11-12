#!/usr/bin/env bash
#SBATCH --job-name=rff_sffast
#SBATCH --partition=c-80-255
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/home/junxuanliu/projects/RNA-Folding-Frontier/logs/seqfree_fast_%j.out
#
# Fixed-0.5 nested scoring for all-N / shuffle dumps.
set -euo pipefail

ROOT="/home/junxuanliu/projects/RNA-Folding-Frontier"
PY="/home/junxuanliu/.conda/envs/rff/bin/python"
cd "${ROOT}"

"${PY}" scripts/score_template_controls.py \
  --root experiments/002_e1prime_runs \
  --threshold 0.5 --min-true-pairs 4 \
  --out experiments/002_e1prime_runs/template_controls.json
