#!/usr/bin/env bash
#SBATCH --job-name=rff_d2
#SBATCH --partition=c-80-255
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/home/junxuanliu/projects/RNA-Folding-Frontier/logs/p2_d2_%j.out
#
# H2.1/H2.2 on the native (a) seed17 dumps, under v1-all and v3-primary.
set -euo pipefail

ROOT="/home/junxuanliu/projects/RNA-Folding-Frontier"
PY="/home/junxuanliu/.conda/envs/rff/bin/python"
cd "${ROOT}"

DUMP="experiments/002_e1prime_runs/dump_seed17"
OUT="experiments/002_e1prime_runs"

for SPLIT in id_test test_withheld; do
  "${PY}" scripts/d2_structure_errors.py \
    --dump-dir "${DUMP}" --names "a_${SPLIT}" \
    --split "splits/s1prime_filtered/${SPLIT}.jsonl" \
    --label-version v1 --min-true-pairs 0 \
    --out "${OUT}/d2_seed17_${SPLIT}_v1all.json"
  "${PY}" scripts/d2_structure_errors.py \
    --dump-dir "${DUMP}" --names "a_${SPLIT}" \
    --split "splits/s1prime_filtered/${SPLIT}.jsonl" \
    --label-version v3 --min-true-pairs 4 \
    --out "${OUT}/d2_seed17_${SPLIT}_v3primary.json"
done
