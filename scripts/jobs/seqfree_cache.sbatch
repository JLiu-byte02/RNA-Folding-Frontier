#!/usr/bin/env bash
#SBATCH --job-name=rff_seqfree
#SBATCH --partition=g-3090-2
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --output=/home/junxuanliu/projects/RNA-Folding-Frontier/logs/seqfree_%j.out
#
# Build all-N or dinucleotide-shuffled caches for id_test and withheld.
#   sbatch scripts/jobs/seqfree_cache.sbatch <all_n|shuffle>
set -euo pipefail

ROOT="/home/junxuanliu/projects/RNA-Folding-Frontier"
MODE="${1:?all_n or shuffle}"
PY="/home/junxuanliu/.conda/envs/rff/bin/python"
MODEL="data/models/rinalmo-giga"

cd "${ROOT}"
for SPLIT in id_test test_withheld; do
  DATA="data/processed/controls/${MODE}_${SPLIT}.jsonl"
  "${PY}" scripts/build_sequence_free_controls.py \
    --input "splits/s1prime_filtered/${SPLIT}.jsonl" \
    --out "${DATA}" --mode "${MODE}" --seed 17
  "${PY}" scripts/build_cache.py \
    --backbone rinalmo-giga --model-path "${MODEL}" \
    --k full --data "${DATA}" --out "data/cache/giga_${MODE}_${SPLIT}" \
    --max-len 1024
done
