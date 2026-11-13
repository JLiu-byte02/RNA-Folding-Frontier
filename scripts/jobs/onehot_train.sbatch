#!/usr/bin/env bash
#SBATCH --job-name=rff_onehothead
#SBATCH --partition=g-3090-2
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/home/junxuanliu/projects/RNA-Folding-Frontier/logs/onehothead_%j.out
#
# Train one no-pretraining one-hot PairMapHead for E-C1b.
#   sbatch scripts/jobs/onehot_train.sbatch <17|23|41>
set -euo pipefail

ROOT="/home/junxuanliu/projects/RNA-Folding-Frontier"
PY="/home/junxuanliu/.conda/envs/rff/bin/python"
SEED="${1:?seed}"
cd "${ROOT}"
test -f data/cache/onehot_train/manifest.json

"${PY}" scripts/train_head_s1.py \
  --train-cache data/cache/onehot_train \
  --out "runs/onehot_seed${SEED}" --seed "${SEED}"
