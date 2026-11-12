#!/usr/bin/env bash
#SBATCH --job-name=rff_traindump
#SBATCH --partition=g-3090-2
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/home/junxuanliu/projects/RNA-Folding-Frontier/logs/train_dump_%j.out
#
# Native-head probability dump on the frozen training split for the (s,b) fit.
set -euo pipefail

ROOT="/home/junxuanliu/projects/RNA-Folding-Frontier"
PY="/home/junxuanliu/.conda/envs/rff/bin/python"
cd "${ROOT}"

"${PY}" scripts/e1prime_dump.py \
  --checkpoint runs/p2_seed17/best.pt \
  --out-dir experiments/002_e1prime_runs/dump_seed17 \
  --eval "train=splits/s1prime_filtered/train.jsonl:data/cache/giga_full_train"
