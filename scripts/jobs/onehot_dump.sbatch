#!/usr/bin/env bash
#SBATCH --job-name=rff_onehotdump
#SBATCH --partition=g-3090-2
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/home/junxuanliu/projects/RNA-Folding-Frontier/logs/onehotdump_%j.out
#
# Dump one-hot head probabilities for E-C1b.
#   sbatch scripts/jobs/onehot_dump.sbatch <17|23|41>
set -euo pipefail

ROOT="/home/junxuanliu/projects/RNA-Folding-Frontier"
PY="/home/junxuanliu/.conda/envs/rff/bin/python"
SEED="${1:?seed}"
CKPT="runs/onehot_seed${SEED}/best.pt"
OUT="experiments/002_e1prime_runs/dump_onehot_seed${SEED}"
cd "${ROOT}"
test -f "${CKPT}"
mkdir -p "${OUT}"

"${PY}" scripts/e1prime_dump.py \
  --checkpoint "${CKPT}" --out-dir "${OUT}" \
  --eval "test_withheld=splits/s1prime_filtered/test_withheld.jsonl:data/cache/onehot_test_withheld" \
  --eval "id_test=splits/s1prime_filtered/id_test.jsonl:data/cache/onehot_id_test"
