#!/usr/bin/env bash
#SBATCH --job-name=rff_p3adump
#SBATCH --partition=g-3090-2
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/home/junxuanliu/projects/RNA-Folding-Frontier/logs/p3adump_%j.out
#
# Dump the retrained P3 condition-a arm on id_test and withheld.
#   sbatch scripts/jobs/p3_dump_a.sbatch <17|23|41>
set -euo pipefail

ROOT="/home/junxuanliu/projects/RNA-Folding-Frontier"
PY="/home/junxuanliu/.conda/envs/rff/bin/python"
SEED="${1:?seed}"
CKPT="runs/p3_a_seed${SEED}/best.pt"
OUT="experiments/002_e1prime_runs/dump_p3_a_seed${SEED}"

cd "${ROOT}"
test -f "${CKPT}"
mkdir -p "${OUT}"

"${PY}" scripts/e1prime_dump.py \
  --checkpoint "${CKPT}" \
  --out-dir "${OUT}" \
  --eval "id_test=splits/s1prime_filtered/id_test.jsonl:data/cache/giga_a_id_test" \
  --eval "test_withheld=splits/s1prime_filtered/test_withheld.jsonl:data/cache/giga_a_test_withheld"
