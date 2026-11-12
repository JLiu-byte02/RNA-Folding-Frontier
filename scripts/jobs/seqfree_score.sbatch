#!/usr/bin/env bash
#SBATCH --job-name=rff_sfscore
#SBATCH --partition=c-80-255
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/home/junxuanliu/projects/RNA-Folding-Frontier/logs/seqfree_score_%j.out
#
# D1-style scoring of all-N and shuffle dumps; v3 primary label policy.
set -euo pipefail

ROOT="/home/junxuanliu/projects/RNA-Folding-Frontier"
PY="/home/junxuanliu/.conda/envs/rff/bin/python"
cd "${ROOT}"

for SEED in 17 23 41; do
  for NAME in all_n_id_test all_n_test_withheld shuffle_id_test shuffle_test_withheld; do
    "${PY}" scripts/d1_scoring.py \
      --dump-dir "experiments/002_e1prime_runs/dump_seqfree_seed${SEED}" \
      --names "${NAME}" \
      --label-version v3 --min-true-pairs 4 --truth-mode nested \
      --out "experiments/002_e1prime_runs/d1_seqfree_seed${SEED}_${NAME}.json"
  done
done
