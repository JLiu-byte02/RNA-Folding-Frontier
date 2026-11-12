#!/usr/bin/env bash
#SBATCH --job-name=rff_p2fix
#SBATCH --partition=g-3090-2
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/home/junxuanliu/projects/RNA-Folding-Frontier/logs/p2fix_%j.out
#
# Dump corrected native-head probabilities for P2-central a-e.
#   sbatch scripts/jobs/p2central_dump_fix.sbatch <17|23|41>
set -euo pipefail

ROOT="/home/junxuanliu/projects/RNA-Folding-Frontier"
SEED="${1:?seed}"
PY="/home/junxuanliu/.conda/envs/rff/bin/python"
OUT="experiments/002_e1prime_runs/dump_p2fix_seed${SEED}"
ID="splits/s1prime_filtered/id_test.jsonl"
WH="splits/s1prime_filtered/test_withheld.jsonl"

cd "${ROOT}"
"${PY}" scripts/e1prime_dump.py \
  --checkpoint "runs/p2fix_seed${SEED}/best.pt" \
  --out-dir "${OUT}" \
  --eval "a_id_test=${ID}:data/cache/giga_a_id_test" \
  --eval "b_id_test=${ID}:data/cache/giga_b_id_test" \
  --eval "c_id_test=${ID}:data/cache/giga_c_id_test" \
  --eval "d_id_test=${ID}:data/cache/giga_d_id_test" \
  --eval "e_id_test=${ID}:data/cache/giga_e_id_test" \
  --eval "a_test_withheld=${WH}:data/cache/giga_a_test_withheld" \
  --eval "b_test_withheld=${WH}:data/cache/giga_b_test_withheld" \
  --eval "c_test_withheld=data/processed/conditions/test_withheld/csubset_native.jsonl:data/cache/giga_c_test_withheld" \
  --eval "d_test_withheld=${WH}:data/cache/giga_d_test_withheld" \
  --eval "e_test_withheld=${WH}:data/cache/giga_e_test_withheld"
