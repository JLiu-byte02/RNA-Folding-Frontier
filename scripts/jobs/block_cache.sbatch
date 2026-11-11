#!/usr/bin/env bash
#SBATCH --job-name=rff_block
#SBATCH --partition=g-3090-2
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=08:00:00
# Absolute: Slurm starts batch jobs in $HOME, not the submit directory.
#SBATCH --output=/home/junxuanliu/projects/RNA-Folding-Frontier/logs/block_%j.out
#
# Build a P2-block embedding cache (Amendment 2 A2.2).
#   sbatch scripts/jobs/block_cache.sbatch <condition> <data.jsonl> <donors.json> <outdir> [extra args]
#
# Condition a needs no donors: with donor == native every pass returns the
# native sequence, so the existing full cache (giga_a_*) is the same object.
set -euo pipefail

ROOT="/home/junxuanliu/projects/RNA-Folding-Frontier"
PY="/home/junxuanliu/.conda/envs/rff/bin/python"
COND="${1:?condition}"
DATA="${2:?data jsonl}"
DONORS="${3:?donors json}"
OUT="${4:?out dir}"
shift 4

cd "${ROOT}"
"${PY}" scripts/build_block_cache.py \
  --condition "${COND}" --data "${DATA}" --donors "${DONORS}" --out "${OUT}" \
  --backbone rinalmo-giga --model-path data/models/rinalmo-giga \
  --block 16 --flank 32 --max-len 1024 "$@"
