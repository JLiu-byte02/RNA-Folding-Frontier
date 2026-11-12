#!/usr/bin/env bash
#SBATCH --job-name=rff_d1n
#SBATCH --partition=c-80-255
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/home/junxuanliu/projects/RNA-Folding-Frontier/logs/p2_d1n_%j.out
#
# Primary-style D1: truth is reduced to its maximum nested subset before
# scoring, as required by Amendment 6. Covers the label contrast too.
set -euo pipefail

ROOT="/home/junxuanliu/projects/RNA-Folding-Frontier"
PY="/home/junxuanliu/.conda/envs/rff/bin/python"
cd "${ROOT}"

for SEED in 17 23 41; do
  for SPLIT in id_test test_withheld; do
    "${PY}" scripts/d1_scoring.py \
      --dump-dir "experiments/002_e1prime_runs/dump_seed${SEED}" \
      --names "a_${SPLIT}" "e_${SPLIT}" \
      --label-version v3 --min-true-pairs 4 --truth-mode nested \
      --out "experiments/002_e1prime_runs/d1_seed${SEED}_${SPLIT}_v3primary_nested.json"
    "${PY}" scripts/d1_scoring.py \
      --dump-dir "experiments/002_e1prime_runs/dump_seed${SEED}" \
      --names "a_${SPLIT}" "e_${SPLIT}" \
      --label-version v1 --min-true-pairs 0 --truth-mode nested \
      --out "experiments/002_e1prime_runs/d1_seed${SEED}_${SPLIT}_v1all_nested.json"
  done
done
