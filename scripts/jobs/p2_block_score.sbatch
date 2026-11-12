#!/usr/bin/env bash
#SBATCH --job-name=rff_blockscore
#SBATCH --partition=g-3090-2
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/home/junxuanliu/projects/RNA-Folding-Frontier/logs/blockscore_%j.out
#
# CPU-only H3 table for one split; no GPU is requested.
#   sbatch scripts/jobs/p2_block_score.sbatch <id_test|test_withheld>
set -euo pipefail

ROOT="/home/junxuanliu/projects/RNA-Folding-Frontier"
PY="/home/junxuanliu/.conda/envs/rff/bin/python"
SPLIT="${1:?id_test or test_withheld}"
cd "${ROOT}"

"${PY}" scripts/score_p2_block.py \
  --split "${SPLIT}" \
  --out "experiments/002_e1prime_runs/p2block_h3_${SPLIT}.json" \
  --md-out "docs/p2block_h3_${SPLIT}_2026-11-12.md"

