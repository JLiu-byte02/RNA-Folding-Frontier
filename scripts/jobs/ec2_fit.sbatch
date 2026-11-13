#!/usr/bin/env bash
#SBATCH --job-name=rff_ec2fit
#SBATCH --partition=g-3090-2
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=/home/junxuanliu/projects/RNA-Folding-Frontier/logs/ec2_fit_%j.out
#
# Amendment 13 A13.2: fit the E-C2 coverage-gated decoder on the TRAIN split
# only (clan-grouped 5-fold CV), freeze + hash before any withheld scoring.
# CPU-only: no --gres is requested, so the job lands on a GPU node but uses no GPU.
set -euo pipefail

ROOT="/home/junxuanliu/projects/RNA-Folding-Frontier"
PY="/home/junxuanliu/.conda/envs/rff/bin/python"
cd "${ROOT}"

"${PY}" scripts/fit_ec2_gate.py \
  --train-dump experiments/002_e1prime_runs/dump_p2fix_seed17/train.jsonl \
  --train-split splits/s1prime_filtered/train.jsonl \
  --pppl experiments/002_e1prime_runs/ec2_train_pppl.json \
  --clans data/raw/rfam/clan_membership.txt \
  --folds 5 --seed 17 --min-true-pairs 4 \
  --workers "${SLURM_CPUS_PER_TASK}" \
  --out experiments/002_e1prime_runs/ec2_fit.json
