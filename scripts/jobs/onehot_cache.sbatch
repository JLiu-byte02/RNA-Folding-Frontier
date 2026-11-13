#!/usr/bin/env bash
#SBATCH --job-name=rff_onehot
#SBATCH --partition=g-3090-2
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/home/junxuanliu/projects/RNA-Folding-Frontier/logs/onehot_%j.out
#
# Build one-hot caches for train/id_test/withheld.
set -euo pipefail

ROOT="/home/junxuanliu/projects/RNA-Folding-Frontier"
PY="/home/junxuanliu/.conda/envs/rff/bin/python"
cd "${ROOT}"

for SPLIT in train id_test test_withheld; do
  "${PY}" scripts/build_onehot_cache.py \
    --split "splits/s1prime_filtered/${SPLIT}.jsonl" \
    --out "data/cache/onehot_${SPLIT}" --max-len 1024
done
