#!/usr/bin/env bash
#SBATCH --job-name=rff_sfdump
#SBATCH --partition=g-3090-2
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --output=/home/junxuanliu/projects/RNA-Folding-Frontier/logs/seqfree_dump_%j.out
#
# Native-head probability dumps on all-N and dinucleotide-shuffle inputs.
#   sbatch scripts/jobs/seqfree_dump.sbatch <seed>
set -euo pipefail

ROOT="/home/junxuanliu/projects/RNA-Folding-Frontier"
SEED="${1:?seed}"
PY="/home/junxuanliu/.conda/envs/rff/bin/python"
OUT="experiments/002_e1prime_runs/dump_seqfree_seed${SEED}"

cd "${ROOT}"
"${PY}" scripts/e1prime_dump.py \
  --checkpoint "runs/p2_seed${SEED}/best.pt" \
  --out-dir "${OUT}" \
  --eval "all_n_id_test=splits/s1prime_filtered/id_test.jsonl:data/cache/giga_all_n_id_test" \
  --eval "all_n_test_withheld=splits/s1prime_filtered/test_withheld.jsonl:data/cache/giga_all_n_test_withheld" \
  --eval "shuffle_id_test=splits/s1prime_filtered/id_test.jsonl:data/cache/giga_shuffle_id_test" \
  --eval "shuffle_test_withheld=splits/s1prime_filtered/test_withheld.jsonl:data/cache/giga_shuffle_test_withheld"
