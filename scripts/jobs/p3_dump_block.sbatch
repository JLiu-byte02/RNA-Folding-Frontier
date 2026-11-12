#!/usr/bin/env bash
#SBATCH --job-name=rff_p3dump
#SBATCH --partition=g-3090-2
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/home/junxuanliu/projects/RNA-Folding-Frontier/logs/p3dump_%j.out
#
# Dump one P3 block-trained head on id_test and withheld.
#   sbatch scripts/jobs/p3_dump_block.sbatch <b|e> <17|23|41>
set -euo pipefail

ROOT="/home/junxuanliu/projects/RNA-Folding-Frontier"
PY="/home/junxuanliu/.conda/envs/rff/bin/python"
COND="${1:?b or e}"
SEED="${2:?seed}"
CKPT="runs/p3_${COND}_seed${SEED}/best.pt"
OUT="experiments/002_e1prime_runs/dump_p3_${COND}_seed${SEED}"

cd "${ROOT}"
test -f "${CKPT}"
mkdir -p "${OUT}"

"${PY}" scripts/e1prime_dump.py \
  --checkpoint "${CKPT}" \
  --out-dir "${OUT}" \
  --eval "id_test=splits/s1prime_filtered/id_test.jsonl:data/cache/giga_block_${COND}_id_test" \
  --eval "test_withheld=splits/s1prime_filtered/test_withheld.jsonl:data/cache/giga_block_${COND}_test_withheld"
