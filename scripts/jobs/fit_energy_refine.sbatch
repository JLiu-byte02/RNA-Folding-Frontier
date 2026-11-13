#!/usr/bin/env bash
#SBATCH --job-name=rff_fitref
#SBATCH --partition=g-3090-2
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=/home/junxuanliu/projects/RNA-Folding-Frontier/logs/fit_energy_refine_%j.out
#
# Amendment 7 A7.4 refinement: the coarse optimum was s=0.25 (lower edge),
# b=0. Span +/- one coarse step; s is kept positive.
set -euo pipefail

ROOT="/home/junxuanliu/projects/RNA-Folding-Frontier"
PY="/home/junxuanliu/.conda/envs/rff/bin/python"
cd "${ROOT}"

"${PY}" scripts/fit_energy_calibration.py \
  --train-dump experiments/002_e1prime_runs/dump_p2fix_seed17/train.jsonl \
  --train-split splits/s1prime_filtered/train.jsonl \
  --grid-s=0.125,0.1875,0.25,0.3125,0.5 \
  --grid-b=-3,-1.5,0,1.5,3 \
  --folds 5 --seed 17 --workers 8 \
  --label-version v3 --min-true-pairs 4 --truth-mode nested \
  --out experiments/002_e1prime_runs/energy_fit_refined_p2fix.json
