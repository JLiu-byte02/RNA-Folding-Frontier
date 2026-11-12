#!/usr/bin/env bash
#SBATCH --job-name=rff_p2block
#SBATCH --partition=g-3090-2
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=/home/junxuanliu/projects/RNA-Folding-Frontier/logs/p2block_%j.out
#
# Build P2-block caches for one split and conditions b/c/d/e.
#   sbatch scripts/jobs/p2_block_split.sbatch <id_test|test_withheld>
#
# Condition a is the existing full-context cache and needs no extra pass.
set -euo pipefail

ROOT="/home/junxuanliu/projects/RNA-Folding-Frontier"
PY="/home/junxuanliu/.conda/envs/rff/bin/python"
SPLIT="${1:?id_test or test_withheld}"
DATA="splits/s1prime_filtered/${SPLIT}.jsonl"
DONOR_DIR="data/processed/block_donors/${SPLIT}"
POOL="splits/s1prime_filtered/train.jsonl"

cd "${ROOT}"
mkdir -p "${DONOR_DIR}"

"${PY}" scripts/build_block_donors.py \
  --split "${DATA}" --out-dir "${DONOR_DIR}" \
  --conditions b,d,e --seed 17 --donor-pool "${POOL}"

"${PY}" scripts/build_block_donors.py \
  --split "${DATA}" --out-dir "${DONOR_DIR}" \
  --conditions c --seed 17 --donor-pool "${POOL}" \
  --homolog-map "data/processed/homologs/${SPLIT}.json"

for COND in b c d e; do
  "${PY}" scripts/build_block_cache.py \
    --condition "${COND}" --data "${DATA}" \
    --donors "${DONOR_DIR}/${COND}.json" \
    --out "data/cache/giga_block_${COND}_${SPLIT}" \
    --backbone rinalmo-giga --model-path data/models/rinalmo-giga \
    --block 16 --flank 32 --max-len 1024
done
