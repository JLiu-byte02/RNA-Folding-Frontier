#!/usr/bin/env bash
#SBATCH --job-name=rff_cache
#SBATCH --partition=g-3090-2
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
# Absolute: Slurm starts batch jobs in $HOME, not the submit directory.
#SBATCH --output=/home/junxuanliu/projects/RNA-Folding-Frontier/logs/cache_%j.out
#
# Build a frozen-backbone embedding cache for one dataset file.
#   sbatch scripts/jobs/cache_emb.sbatch <backbone> <k|full> <data.jsonl> <outdir>
set -euo pipefail

# Slurm spools this script, so BASH_SOURCE does not point into the repo; and
# batch jobs start in $HOME, where relative paths fail. Use absolute paths.
ROOT="/home/junxuanliu/projects/RNA-Folding-Frontier"
BACKBONE="${1:?backbone}"
K="${2:?k or full}"
DATA="${3:?dataset}"
OUT="${4:?outdir}"
MODEL="${RFF_MODEL_DIR:-data/models/rinalmo-giga}"

cd "${ROOT}"
if [ "${BACKBONE}" != "rinalmo-giga" ]; then
  MODEL="data/models/${BACKBONE}"
fi

/home/junxuanliu/.conda/envs/rff/bin/python scripts/build_cache.py \
  --backbone "${BACKBONE}" --model-path "${MODEL}" \
  --k "${K}" --data "${DATA}" --out "${OUT}" --max-len 1024
