#!/usr/bin/env bash
#SBATCH --job-name=rff_mx_smoke
#SBATCH --partition=g-3090-2
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --gres=gpu:1
#SBATCH --exclude=gpu0004
#SBATCH --output=/home/junxuanliu/projects/RNA-Folding-Frontier/logs/mx_smoke_%j.out
#
# Smoke test: does `mxfold2 train` run end-to-end on our exported BPSEQ data?
set -euo pipefail

ROOT="/home/junxuanliu/projects/RNA-Folding-Frontier"
PY="/home/junxuanliu/.conda/envs/rff/bin/python"
MXP="/home/junxuanliu/.conda/envs/mxfold2/bin/python"
export PYTHONPATH="/home/junxuanliu/projects/vendor/MXfold2"
cd "${ROOT}"

OUT="experiments/013_mxfold2/smoke"
mkdir -p "${OUT}" logs

nvidia-smi --query-gpu=name,memory.total,memory.used --format=csv,noheader

"${PY}" scripts/export_mxfold2_bpseq.py \
  --split splits/s1prime_filtered/train.jsonl \
  --out-dir "${OUT}/data" --max-len 120 --limit 300 \
  --manifest "${OUT}/export_manifest.json"

time "${MXP}" -m mxfold2 train \
  --model MixC --epochs 1 --gpu 0 --disable-progress-bar \
  --param "${OUT}/model.pth" --save-config "${OUT}/model.conf" \
  "${OUT}/data/train.lst"

ls -l "${OUT}"