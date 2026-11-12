#!/usr/bin/env bash
#SBATCH --job-name=rff_fit
#SBATCH --partition=g-3090-2
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/home/junxuanliu/projects/RNA-Folding-Frontier/logs/fit_energy_%j.out
#
# Amendment 7 coarse (s,b) fit on the frozen training split.
set -euo pipefail

ROOT="/home/junxuanliu/projects/RNA-Folding-Frontier"
PY="/home/junxuanliu/.conda/envs/rff/bin/python"
cd "${ROOT}"

"${PY}" scripts/fit_energy_calibration.py \
  --train-dump experiments/002_e1prime_runs/dump_p2fix_seed17/train.jsonl \
  --train-split splits/s1prime_filtered/train.jsonl \
  --grid-s=0.25,0.5,1.0,2.0,4.0 \
  --grid-b=-6,-3,0,3,6 \
  --folds 5 --seed 17 --workers 8 \
  --label-version v3 --min-true-pairs 4 --truth-mode nested \
  --out experiments/002_e1prime_runs/energy_fit_coarse_p2fix.json
