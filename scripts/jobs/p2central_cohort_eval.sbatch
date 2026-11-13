#!/usr/bin/env bash
#SBATCH --job-name=rff_p2coh
#SBATCH --partition=g-3090-2
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/home/junxuanliu/projects/RNA-Folding-Frontier/logs/p2coh_%j.out
#
# Evaluate corrected P2-central head dumps on the frozen common cohort,
# one seed per job (the head row of the four-row external-baseline table).
#   sbatch scripts/jobs/p2central_cohort_eval.sbatch <17|23|41>
set -euo pipefail
ROOT="/home/junxuanliu/projects/RNA-Folding-Frontier"
SEED="${1:?seed}"
PY="/home/junxuanliu/.conda/envs/rff/bin/python"
cd "${ROOT}"
"${PY}" scripts/p2central_eval.py \
  --dump-prefix experiments/002_e1prime_runs/dump_p2fix_seed \
  --seeds "${SEED}" \
  --threshold 0.5 --min-true-pairs 4 \
  --restrict-ids experiments/002_e1prime_runs/cohort/all_splits_ids.jsonl \
  --out "experiments/002_e1prime_runs/cohort/p2central_cohort_seed${SEED}.json"
