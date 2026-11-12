#!/usr/bin/env bash
#SBATCH --job-name=rff_block_build
#SBATCH --partition=g-3090-2
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/home/junxuanliu/projects/RNA-Folding-Frontier/logs/block_build_%j.out
#
# Build one P2-block cache for one split and one condition.
#   sbatch scripts/jobs/p2_block_build_one.sbatch <id_test|test_withheld|train> <b|c|d|e>
#
# Condition a is the existing full-context cache and is never built here.
# The output directory is a marker: an existing manifest means the cache is
# complete and the job exits without touching it.
set -euo pipefail

ROOT="/home/junxuanliu/projects/RNA-Folding-Frontier"
PY="/home/junxuanliu/.conda/envs/rff/bin/python"
SPLIT="${1:?id_test, test_withheld or train}"
COND="${2:?b, c, d or e}"
DATA="splits/s1prime_filtered/${SPLIT}.jsonl"
DONOR_DIR="data/processed/block_donors/${SPLIT}"
POOL="splits/s1prime_filtered/train.jsonl"
OUT="data/cache/giga_block_${COND}_${SPLIT}"

cd "${ROOT}"
if [ -f "${OUT}/manifest.json" ]; then
  echo "cache already complete: ${OUT}"
  cat "${OUT}/manifest.json"
  exit 0
fi

mkdir -p "${DONOR_DIR}"
if [ "${COND}" = "c" ]; then
  HOMOLOG="data/processed/homologs/${SPLIT}.json"
  test -f "${HOMOLOG}"
  # The homolog map is built per split and points to another record in the
  # same split; the donor pool must therefore contain those records. Using the
  # train pool here silently produced an empty donor map and an empty cache.
  C_POOL="${DATA}"
  "${PY}" scripts/build_block_donors.py \
    --split "${DATA}" --out-dir "${DONOR_DIR}" \
    --conditions c --seed 17 --donor-pool "${C_POOL}" \
    --homolog-map "${HOMOLOG}"
else
  "${PY}" scripts/build_block_donors.py \
    --split "${DATA}" --out-dir "${DONOR_DIR}" \
    --conditions "${COND}" --seed 17 --donor-pool "${POOL}"
fi

"${PY}" scripts/build_block_cache.py \
  --condition "${COND}" --data "${DATA}" \
  --donors "${DONOR_DIR}/${COND}.json" \
  --out "${OUT}" \
  --backbone rinalmo-giga --model-path data/models/rinalmo-giga \
  --block 16 --flank 32 --max-len 1024

