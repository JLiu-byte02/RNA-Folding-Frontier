#!/usr/bin/env bash
#SBATCH --job-name=rff_p3
#SBATCH --partition=g-3090-2
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/home/junxuanliu/projects/RNA-Folding-Frontier/logs/p3_%j.out
#
# Train one P3 head on a block cache. Evaluation is dumped separately so that
# the v3 label swap and the D1-core bootstrap are computed without a GPU.
#   sbatch scripts/jobs/p3_train_block.sbatch <b|e> <17|23|41>
set -euo pipefail

ROOT="/home/junxuanliu/projects/RNA-Folding-Frontier"
PY="/home/junxuanliu/.conda/envs/rff/bin/python"
COND="${1:?b or e}"
SEED="${2:?seed}"
TRAIN_CACHE="data/cache/giga_block_${COND}_train"
OUT="runs/p3_${COND}_seed${SEED}"

cd "${ROOT}"
test -f "${TRAIN_CACHE}/manifest.json"

"${PY}" scripts/train_head_s1.py \
  --train-cache "${TRAIN_CACHE}" \
  --out "${OUT}" \
  --seed "${SEED}"

