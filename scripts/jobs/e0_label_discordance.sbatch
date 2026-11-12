#!/usr/bin/env bash
#SBATCH --job-name=rff_e0disc
#SBATCH --partition=g-3090-2
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/home/junxuanliu/projects/RNA-Folding-Frontier/logs/e0disc_%j.out
#
# Recompute E0 label-discordance F1 after the pair_prf correction.
set -euo pipefail

ROOT="/home/junxuanliu/projects/RNA-Folding-Frontier"
PY="/home/junxuanliu/.conda/envs/rff/bin/python"
cd "${ROOT}"

"${PY}" scripts/label_discordance.py \
  --archiveii data/raw/archiveII/data.json \
  --bpseq-dir data/raw/bpfold/BPfold_data/bpRNA/TR0 \
  --out experiments/000_e0_leakage/label_discordance_bpfold_archiveII

"${PY}" scripts/label_discordance.py \
  --archiveii data/raw/archiveII/data.json \
  --bpseq-dir data/raw/bpfold/BPfold_data/RNAStrAlign_bpseq \
  --out experiments/000_e0_leakage/label_discordance_rnastralign
