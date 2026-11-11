#!/usr/bin/env bash
#SBATCH --job-name=rff_bench
#SBATCH --partition=g-3090-2
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --time=00:30:00
#SBATCH --output=logs/bench_%j.out

set -euo pipefail
cd /home/junxuanliu/projects/RNA-Folding-Frontier

PY=/home/junxuanliu/.conda/envs/rff/bin/python
export HF_HOME=$PWD/data/models/.hf
export TOKENIZERS_PARALLELISM=false

$PY scripts/bench_backbone.py \
  --backbone rinalmo-giga \
  --path data/models/rinalmo-giga \
  --n-seqs 4 --length 300 \
  --out experiments/001_e1prime/bench_rinalmo_giga.json
