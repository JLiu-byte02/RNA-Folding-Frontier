#!/usr/bin/env bash
#SBATCH --job-name=rff_block_dump
#SBATCH --partition=g-3090-2
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/home/junxuanliu/projects/RNA-Folding-Frontier/logs/block_dump_%j.out
#
# Dump native-head probabilities for all five P2-block conditions.
#   sbatch scripts/jobs/p2_block_dump.sbatch <17|23|41> <id_test|test_withheld>
set -euo pipefail

ROOT="/home/junxuanliu/projects/RNA-Folding-Frontier"
PY="/home/junxuanliu/.conda/envs/rff/bin/python"
SEED="${1:?seed}"
SPLIT="${2:?id_test or test_withheld}"
DATA="splits/s1prime_filtered/${SPLIT}.jsonl"
OUT="experiments/002_e1prime_runs/dump_p2block_seed${SEED}"

cd "${ROOT}"
mkdir -p "${OUT}"

ARGS=()
for COND in a b c d e; do
  if [ "${COND}" = "a" ]; then
    CACHE="data/cache/giga_a_${SPLIT}"
  else
    CACHE="data/cache/giga_block_${COND}_${SPLIT}"
  fi
  test -f "${CACHE}/manifest.json"
  ARGS+=(--eval "${COND}_${SPLIT}=${DATA}:${CACHE}")
done

"${PY}" scripts/e1prime_dump.py \
  --checkpoint "runs/p2fix_seed${SEED}/best.pt" \
  --out-dir "${OUT}" \
  "${ARGS[@]}"

