#!/usr/bin/env bash
#SBATCH --job-name=rff_mx_fold
#SBATCH --partition=g-3090-2
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/home/junxuanliu/projects/RNA-Folding-Frontier/logs/mx_fold_%j.out
#
# Fold both evaluation splits with the retrained MXfold2 (Amendment 13 A13.1).
# CPU only: the MixC forward pass is small and the O(L^3) folding DP dominates,
# so 16 CPU workers beat one GPU.  Resumable via the per-split .partial.jsonl
# cache; a killed job restarts from the last completed chunk.
set -euo pipefail
ROOT="/home/junxuanliu/projects/RNA-Folding-Frontier"
PY="/home/junxuanliu/.conda/envs/rff/bin/python"
OUT="experiments/013_mxfold2"
cd "${ROOT}"

for split in test_withheld id_test; do
  "${PY}" scripts/predict_mxfold2.py \
    --split "splits/s1prime_filtered/${split}.jsonl" \
    --out "${OUT}/mxfold2_${split}.jsonl" \
    --param "${OUT}/run/model.pth" \
    --conf "${OUT}/run/model.conf" \
    --model MixC --max-length 600 \
    --workers 16 --chunk-size 64
done
echo "[mxfold2 fold done] $(date -Iseconds)"