#!/usr/bin/env bash
#SBATCH --job-name=rff_p2eval
#SBATCH --partition=g-3090-2
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/home/junxuanliu/projects/RNA-Folding-Frontier/logs/p2eval_%j.out
#
# Evaluate corrected P2-central dumps for one seed.
#   sbatch scripts/jobs/p2central_eval_one.sbatch <17|23|41>
set -euo pipefail

ROOT="/home/junxuanliu/projects/RNA-Folding-Frontier"
SEED="${1:?seed}"
PY="/home/junxuanliu/.conda/envs/rff/bin/python"
cd "${ROOT}"

"${PY}" scripts/p2central_eval.py \
  --dump-prefix experiments/002_e1prime_runs/dump_p2fix_seed \
  --seeds "${SEED}" \
  --threshold 0.5 --min-true-pairs 4 \
  --out "experiments/002_e1prime_runs/p2central_fixed_seed${SEED}.json"
