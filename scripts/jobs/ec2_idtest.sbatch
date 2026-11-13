#!/usr/bin/env bash
#SBATCH --job-name=rff_ec2id
#SBATCH --partition=g-3090-2
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=/home/junxuanliu/projects/RNA-Folding-Frontier/logs/ec2_idtest_%j.out
#
# E-C2 id_test replication: compute the A12 masked pppl on id_test (never done
# for this split), then re-run the frozen selector on id_test and the
# exploratory withheld scale-control arm. No label-driven choices here.
set -euo pipefail

ROOT="/home/junxuanliu/projects/RNA-Folding-Frontier"
PY="/home/junxuanliu/.conda/envs/rff/bin/python"
cd "${ROOT}"
R="experiments/002_e1prime_runs"

echo "[ec2id] step 1: A12 masked pppl on id_test"
"${PY}" scripts/ec1_p2_pppl.py \
  --split splits/s1prime_filtered/id_test.jsonl \
  --all-split \
  --checkpoint data/models/rinalmo-giga \
  --masks 4 --mask-frac 0.15 --seeds 17 23 41 53 \
  --batch-size 4 \
  --out "${R}/ec2_pppl_id_test.json"

echo "[ec2id] step 2: id_test confirmatory (frozen selector)"
"${PY}" scripts/score_ec2_gate.py \
  --dump "${R}/dump_p2fix_seed17/a_id_test.jsonl" \
  --split splits/s1prime_filtered/id_test.jsonl \
  --pppl "${R}/ec2_pppl_id_test.json" \
  --ref-preds "${R}/rnafold_id_test.jsonl" \
  --workers "${SLURM_CPUS_PER_TASK:-8}" --n-boot 2000 \
  --out "${R}/ec2_id_test_score.json"

echo "[ec2id] step 3: withheld EXPLORATORY scale-confound control"
"${PY}" scripts/score_ec2_gate.py \
  --dump "${R}/dump_p2fix_seed17/a_test_withheld.jsonl" \
  --split splits/s1prime_filtered/test_withheld.jsonl \
  --pppl "${R}/ec1_p2_pppl.json" \
  --fixed-scale 0.09375 \
  --workers "${SLURM_CPUS_PER_TASK:-8}" --n-boot 2000 \
  --out "${R}/ec2_test_withheld_scalecontrol.json"

echo "[ec2id] done"