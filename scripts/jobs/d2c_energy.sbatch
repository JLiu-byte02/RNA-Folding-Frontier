#!/usr/bin/env bash
#SBATCH --job-name=rff_d2c
#SBATCH --partition=c-80-255
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/home/junxuanliu/projects/RNA-Folding-Frontier/logs/d2c_%j.out
#
# Exploratory H2.4: FP-helix stacking energy and RNAfold-MFE overlap.
set -euo pipefail

ROOT="/home/junxuanliu/projects/RNA-Folding-Frontier"
PY="/home/junxuanliu/.conda/envs/rff/bin/python"
cd "${ROOT}"

for SPLIT in id_test test_withheld; do
  "${PY}" scripts/d2c_helix_energy.py \
    --dump-dir experiments/002_e1prime_runs/dump_seed17 \
    --split "splits/s1prime_filtered/${SPLIT}.jsonl" \
    --name "a_${SPLIT}" \
    --out "experiments/002_e1prime_runs/d2c_seed17_${SPLIT}_v3primary.json"
done
