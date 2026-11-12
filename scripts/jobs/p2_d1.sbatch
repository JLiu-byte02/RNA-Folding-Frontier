#!/usr/bin/env bash
#SBATCH --job-name=rff_d1
#SBATCH --partition=c-80-255
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/home/junxuanliu/projects/RNA-Folding-Frontier/logs/p2_d1_%j.out
#
# Score the P2 probability dumps under four label policies.
#   sbatch scripts/jobs/p2_d1.sbatch <seed>
#
# v1all       v1 truth, empty-truth rows included
# v1nonempty  v1 truth, empty-truth rows excluded
# v3all       v3 truth, all rows with a parsed label
# v3primary   v3 truth, Amendment 8 primary rule (<4 true pairs excluded)
set -euo pipefail

ROOT="/home/junxuanliu/projects/RNA-Folding-Frontier"
SEED="${1:?seed}"
PY="/home/junxuanliu/.conda/envs/rff/bin/python"

cd "${ROOT}"
DUMP="experiments/002_e1prime_runs/dump_seed${SEED}"
OUTDIR="experiments/002_e1prime_runs"

for SPLIT in id_test test_withheld; do
  "${PY}" scripts/d1_scoring.py \
    --dump-dir "${DUMP}" --names "a_${SPLIT}" "e_${SPLIT}" \
    --label-version v1 --min-true-pairs 0 \
    --out "${OUTDIR}/d1_seed${SEED}_${SPLIT}_v1all.json"
  "${PY}" scripts/d1_scoring.py \
    --dump-dir "${DUMP}" --names "a_${SPLIT}" "e_${SPLIT}" \
    --label-version v1 --min-true-pairs 1 \
    --out "${OUTDIR}/d1_seed${SEED}_${SPLIT}_v1nonempty.json"
  "${PY}" scripts/d1_scoring.py \
    --dump-dir "${DUMP}" --names "a_${SPLIT}" "e_${SPLIT}" \
    --label-version v3 --min-true-pairs 0 \
    --out "${OUTDIR}/d1_seed${SEED}_${SPLIT}_v3all.json"
  "${PY}" scripts/d1_scoring.py \
    --dump-dir "${DUMP}" --names "a_${SPLIT}" "e_${SPLIT}" \
    --label-version v3 --min-true-pairs 4 \
    --out "${OUTDIR}/d1_seed${SEED}_${SPLIT}_v3primary.json"
done
