#!/usr/bin/env bash
set -euo pipefail
cd /home/junxuanliu/projects/RNA-Folding-Frontier
PY="/home/junxuanliu/.conda/envs/rff/bin/python"
R="experiments/002_e1prime_runs"
MX="experiments/013_mxfold2"
for s in test_withheld id_test; do
  echo "=== ${s}"
  "${PY}" scripts/cohort_band_interaction.py \
    --split "splits/s1prime_filtered/${s}.jsonl" \
    --cohort "${R}/cohort/${s}.jsonl" \
    --method "RNAfold=${R}/rnafold_${s}.jsonl" \
    --method "EternaFold=${R}/eternafold_${s}.jsonl" \
    --method "RNAstructure=${R}/rnastructure_${s}.jsonl" \
    --method "MXfold2=${MX}/mxfold2_${s}.jsonl" \
    --method "head17=${R}/headpreds/seed17_${s}.jsonl" \
    --method "head23=${R}/headpreds/seed23_${s}.jsonl" \
    --method "head41=${R}/headpreds/seed41_${s}.jsonl" \
    --mean-of head_mean=head17,head23,head41 \
    --compare head_mean:RNAfold \
    --compare head_mean:EternaFold \
    --compare head_mean:RNAstructure \
    --compare head_mean:MXfold2 \
    --out "${R}/cohort/band_interaction_${s}.json"
done
echo "[band interaction done] $(date -Iseconds)"
