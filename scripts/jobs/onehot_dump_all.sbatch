#!/usr/bin/env bash
#SBATCH --job-name=rff_onehotdumpall
#SBATCH --partition=g-3090-2
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --exclude=gpu0004
#SBATCH --output=/home/junxuanliu/projects/RNA-Folding-Frontier/logs/onehotdumpall_%j.out
#
# Dump one-hot head probabilities for all three seeds, serially.
set -uo pipefail

ROOT="/home/junxuanliu/projects/RNA-Folding-Frontier"
PY="/home/junxuanliu/.conda/envs/rff/bin/python"
cd "${ROOT}"
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

for SEED in 17 23 41; do
  CKPT="runs/onehot_seed${SEED}/best.pt"
  OUT="experiments/002_e1prime_runs/dump_onehot_seed${SEED}"
  test -f "${CKPT}"
  if [ -f "${OUT}/test_withheld.jsonl" ]; then
    echo "[skip] seed ${SEED}: dump exists"
    continue
  fi
  mkdir -p "${OUT}"
  ok=0
  for ATTEMPT in 1 2 3; do
    echo "[dump] seed ${SEED} attempt ${ATTEMPT} $(date -Is)"
    if "${PY}" scripts/e1prime_dump.py \
         --checkpoint "${CKPT}" --out-dir "${OUT}" \
         --eval "test_withheld=splits/s1prime_filtered/test_withheld.jsonl:data/cache/onehot_test_withheld" \
         --eval "id_test=splits/s1prime_filtered/id_test.jsonl:data/cache/onehot_id_test"; then
      ok=1; break
    fi
    echo "[retry] seed ${SEED} attempt ${ATTEMPT} failed; sleeping 120s"
    sleep 120
  done
  if [ "${ok}" -ne 1 ]; then echo "[fatal] dump seed ${SEED} failed 3x"; exit 1; fi
done
echo "[done] $(date -Is)"
