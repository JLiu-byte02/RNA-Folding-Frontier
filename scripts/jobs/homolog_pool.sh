#!/usr/bin/env bash
# Build condition-(c) homolog maps for every S1' split.
#
# The CPU partition is fully allocated by other users, and this is a small,
# low-priority job, so it runs on the login node at nice 19.  Launch it as
#   nohup bash scripts/jobs/homolog_pool.sh > logs/homolog_pool.log 2>&1 &
# from the project root.
set -euo pipefail

ROOT="/home/junxuanliu/projects/RNA-Folding-Frontier"
PY="/home/junxuanliu/.conda/envs/rff/bin/python"
cd "${ROOT}"

for split in train id_test test_withheld; do
  echo "=== ${split} $(date -Is)"
  "${PY}" scripts/build_homolog_pool.py \
    --input "splits/s1prime_filtered/${split}.jsonl" \
    --out "data/processed/homologs/${split}.json" \
    --max-verify 12
  echo "=== done ${split} $(date -Is)"
done
