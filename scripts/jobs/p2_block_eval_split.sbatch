#!/usr/bin/env bash
#SBATCH --job-name=rff_p2beval
#SBATCH --partition=g-3090-2
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/home/junxuanliu/projects/RNA-Folding-Frontier/logs/p2beval_%j.out
#
# Dump native-head probabilities for the P2-block c condition and run H2.3
# against the existing full-context a dump.
#   sbatch scripts/jobs/p2_block_eval_split.sbatch <id_test|test_withheld>
set -euo pipefail

ROOT="/home/junxuanliu/projects/RNA-Folding-Frontier"
PY="/home/junxuanliu/.conda/envs/rff/bin/python"
SPLIT="${1:?id_test or test_withheld}"
OUT="experiments/002_e1prime_runs/dump_p2block_seed17/${SPLIT}"
DUMP_A="experiments/002_e1prime_runs/dump_seed17/a_${SPLIT}.jsonl"

cd "${ROOT}"
mkdir -p "${OUT}"

"${PY}" scripts/e1prime_dump.py \
  --checkpoint runs/p2_seed17/best.pt \
  --out-dir "${OUT}" \
  --eval "c_${SPLIT}=splits/s1prime_filtered/${SPLIT}.jsonl:data/cache/giga_block_c_${SPLIT}"

"${PY}" scripts/d2b_block_boundary.py \
  --a-dump "${DUMP_A}" \
  --c-dump "${OUT}/c_${SPLIT}.jsonl" \
  --split "splits/s1prime_filtered/${SPLIT}.jsonl" \
  --name "p2block_${SPLIT}_a_vs_c" \
  --block 16 --flank 32 \
  --label-version v3 --min-true-pairs 4 \
  --out "experiments/002_e1prime_runs/h23_${SPLIT}.json"
