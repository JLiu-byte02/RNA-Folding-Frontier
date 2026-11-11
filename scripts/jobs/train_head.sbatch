#!/usr/bin/env bash
#SBATCH --job-name=rff_head
#SBATCH --partition=g-3090-2
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=12:00:00
# Absolute: Slurm starts batch jobs in $HOME, not the submit directory.
#SBATCH --output=/home/junxuanliu/projects/RNA-Folding-Frontier/logs/head_%j.out
#
# Train the decoder-B pair head on S1' condition (a) and evaluate the P2 table.
#   sbatch scripts/jobs/train_head.sbatch <seed> <out-dir> [extra args...]
set -euo pipefail

ROOT="/home/junxuanliu/projects/RNA-Folding-Frontier"
PY="/home/junxuanliu/.conda/envs/rff/bin/python"
SEED="${1:?seed}"
OUT="${2:?out dir}"
shift 2

cd "${ROOT}"
"${PY}" scripts/train_head_s1.py \
  --train splits/s1prime_filtered/train.jsonl \
  --train-cache data/cache/giga_a_train \
  --seed "${SEED}" \
  --out "${OUT}" \
  "$@"
