#!/usr/bin/env bash
#SBATCH --job-name=rff_p3a
#SBATCH --partition=g-3090-2
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/home/junxuanliu/projects/RNA-Folding-Frontier/logs/p3a_%j.out
#
# Retrain the P3 condition-a arm on the common split (Amendment 11).
#   sbatch scripts/jobs/p3_train_a.sbatch <17|23|41>
set -euo pipefail

ROOT="/home/junxuanliu/projects/RNA-Folding-Frontier"
PY="/home/junxuanliu/.conda/envs/rff/bin/python"
SEED="${1:?seed}"

cd "${ROOT}"
test -f data/cache/giga_a_train/manifest.json

"${PY}" scripts/train_head_s1.py \
  --train-cache data/cache/giga_a_train \
  --out "runs/p3_a_seed${SEED}" \
  --seed "${SEED}"
