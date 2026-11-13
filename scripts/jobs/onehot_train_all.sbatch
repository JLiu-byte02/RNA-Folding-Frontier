#!/usr/bin/env bash
#SBATCH --job-name=rff_onehotall
#SBATCH --partition=g-3090-2
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --exclude=gpu0004
#SBATCH --output=/home/junxuanliu/projects/RNA-Folding-Frontier/logs/onehotall_%j.out
#
# Train all three one-hot PairMapHead seeds serially in one allocation.
# Rationale: three parallel submissions on one node collided with another
# user's processes (cuBLAS handle alloc failures surfaced as
# "GET was unable to find an engine").  Serial + retry is robust and cheap.
set -uo pipefail

ROOT="/home/junxuanliu/projects/RNA-Folding-Frontier"
PY="/home/junxuanliu/.conda/envs/rff/bin/python"
cd "${ROOT}"
test -f data/cache/onehot_train/manifest.json
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

for SEED in 17 23 41; do
  if [ -f "runs/onehot_seed${SEED}/best.pt" ]; then
    echo "[skip] seed ${SEED}: best.pt exists"
    continue
  fi
  ok=0
  for ATTEMPT in 1 2 3; do
    echo "[train] seed ${SEED} attempt ${ATTEMPT} $(date -Is)"
    if "${PY}" scripts/train_head_s1.py \
         --train-cache data/cache/onehot_train \
         --out "runs/onehot_seed${SEED}" --seed "${SEED}"; then
      ok=1; break
    fi
    echo "[retry] seed ${SEED} attempt ${ATTEMPT} failed; sleeping 120s"
    sleep 120
  done
  if [ "${ok}" -ne 1 ]; then echo "[fatal] seed ${SEED} failed 3x"; exit 1; fi
done
echo "[done] $(date -Is)"
