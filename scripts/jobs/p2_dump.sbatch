#!/usr/bin/env bash
#SBATCH --job-name=rff_dump
#SBATCH --partition=g-3090-2
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
# Absolute: Slurm starts batch jobs in $HOME.
#SBATCH --output=/home/junxuanliu/projects/RNA-Folding-Frontier/logs/p2_dump_%j.out
#
# Dump native (a) and masked (e) pair probabilities for D1/K1.
#   sbatch scripts/jobs/p2_dump.sbatch <seed>
#
# The dump carries both v1 and v3 truth. Scoring chooses the label policy, so
# retraining is not needed to re-score the same probabilities.
set -euo pipefail

ROOT="/home/junxuanliu/projects/RNA-Folding-Frontier"
SEED="${1:?seed}"
PY="/home/junxuanliu/.conda/envs/rff/bin/python"

cd "${ROOT}"
OUT="experiments/002_e1prime_runs/dump_seed${SEED}"

"${PY}" scripts/e1prime_dump.py \
  --checkpoint "runs/p2_seed${SEED}/best.pt" \
  --out-dir "${OUT}" \
  --eval "a_id_test=splits/s1prime_filtered/id_test.jsonl:data/cache/giga_a_id_test" \
  --eval "e_id_test=splits/s1prime_filtered/id_test.jsonl:data/cache/giga_e_id_test" \
  --eval "a_test_withheld=splits/s1prime_filtered/test_withheld.jsonl:data/cache/giga_a_test_withheld" \
  --eval "e_test_withheld=splits/s1prime_filtered/test_withheld.jsonl:data/cache/giga_e_test_withheld"
