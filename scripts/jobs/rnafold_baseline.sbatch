#!/bin/bash
#SBATCH -J rff_rnafold
#SBATCH -p g-3090-2
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH -t 8:00:00
#SBATCH -o logs/rnafold_%j.out
set -euo pipefail
cd /home/junxuanliu/projects/RNA-Folding-Frontier
export PATH=/home/junxuanliu/.conda/envs/rff/bin:$PATH
export OMP_NUM_THREADS=1
echo "[start] $(date -Iseconds)"
python scripts/score_external_baseline.py predict \
  --split splits/s1prime_filtered/test_withheld.jsonl \
  --out experiments/002_e1prime_runs/rnafold_test_withheld.jsonl \
  --workers 16
echo "[predict withheld done] $(date -Iseconds)"
python scripts/score_external_baseline.py predict \
  --split splits/s1prime_filtered/id_test.jsonl \
  --out experiments/002_e1prime_runs/rnafold_id_test.jsonl \
  --workers 16
echo "[predict id_test done] $(date -Iseconds)"
python scripts/score_external_baseline.py score \
  --split splits/s1prime_filtered/test_withheld.jsonl \
  --preds experiments/002_e1prime_runs/rnafold_test_withheld.jsonl \
  --restrict-ids experiments/002_e1prime_runs/dump_p2fix_seed17/a_test_withheld.jsonl \
  --out experiments/002_e1prime_runs/rnafold_test_withheld_score.json \
  --per-family-out experiments/002_e1prime_runs/rnafold_test_withheld_by_family.json
python scripts/score_external_baseline.py score \
  --split splits/s1prime_filtered/id_test.jsonl \
  --preds experiments/002_e1prime_runs/rnafold_id_test.jsonl \
  --restrict-ids experiments/002_e1prime_runs/dump_p2fix_seed17/a_id_test.jsonl \
  --out experiments/002_e1prime_runs/rnafold_id_test_score.json \
  --per-family-out experiments/002_e1prime_runs/rnafold_id_test_by_family.json
echo "[all done] $(date -Iseconds)"