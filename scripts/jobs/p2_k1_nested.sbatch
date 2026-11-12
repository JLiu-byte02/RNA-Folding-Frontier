#!/usr/bin/env bash
#SBATCH --job-name=rff_k1n
#SBATCH --partition=c-80-255
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/home/junxuanliu/projects/RNA-Folding-Frontier/logs/p2_k1n_%j.out
#
# K1 on the primary-style (maximum-nested-subset truth) D1 outputs.
set -euo pipefail

ROOT="/home/junxuanliu/projects/RNA-Folding-Frontier"
PY="/home/junxuanliu/.conda/envs/rff/bin/python"
cd "${ROOT}"

for SPLIT in id_test test_withheld; do
  for POLICY in v1all_nested v3primary_nested; do
    "${PY}" scripts/k1_check.py \
      --a \
        "experiments/002_e1prime_runs/d1_seed17_${SPLIT}_${POLICY}.json.per_sequence/a_${SPLIT}.jsonl" \
        "experiments/002_e1prime_runs/d1_seed23_${SPLIT}_${POLICY}.json.per_sequence/a_${SPLIT}.jsonl" \
        "experiments/002_e1prime_runs/d1_seed41_${SPLIT}_${POLICY}.json.per_sequence/a_${SPLIT}.jsonl" \
      --e \
        "experiments/002_e1prime_runs/d1_seed17_${SPLIT}_${POLICY}.json.per_sequence/e_${SPLIT}.jsonl" \
        "experiments/002_e1prime_runs/d1_seed23_${SPLIT}_${POLICY}.json.per_sequence/e_${SPLIT}.jsonl" \
        "experiments/002_e1prime_runs/d1_seed41_${SPLIT}_${POLICY}.json.per_sequence/e_${SPLIT}.jsonl" \
      --metric-label "nested-truth primary F1, oracle threshold, 80-299 nt, ${SPLIT}, ${POLICY}" \
      --out "experiments/002_e1prime_runs/k1_${SPLIT}_${POLICY}.json"
  done
done
