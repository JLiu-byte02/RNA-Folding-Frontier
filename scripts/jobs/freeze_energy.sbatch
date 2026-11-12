#!/usr/bin/env bash
#SBATCH --job-name=rff_freeze
#SBATCH --partition=g-3090-2
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=/home/junxuanliu/projects/RNA-Folding-Frontier/logs/freeze_energy_%j.out
#
# Validate and lock the (s,b) fit before H2.4.
set -euo pipefail

ROOT="/home/junxuanliu/projects/RNA-Folding-Frontier"
PY="/home/junxuanliu/.conda/envs/rff/bin/python"
cd "${ROOT}"

"${PY}" scripts/freeze_energy_fit.py \
  --fit experiments/002_e1prime_runs/energy_fit_coarse_p2fix.json \
  --lockfile state/amendments.lock.json \
  --write-lock

