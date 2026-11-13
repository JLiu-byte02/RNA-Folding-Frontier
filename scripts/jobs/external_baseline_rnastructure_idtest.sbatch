#!/usr/bin/env bash
#SBATCH -J rff_rs_idtest
#SBATCH -p g-3090-2
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH -t 4:00:00
#SBATCH -o logs/rs_idtest_%j.out
# CPU-only companion to external_baseline_rnastructure.sbatch: the withheld
# fold is rate-limited by long sequences, so the much smaller id_test split is
# folded in parallel into its own cache file (no shared append target).
set -euo pipefail
cd /home/junxuanliu/projects/RNA-Folding-Frontier
export PATH=/home/junxuanliu/.conda/envs/rff/bin:$PATH
echo "[start] $(date -Iseconds)"
python scripts/predict_rnastructure.py \
  --split splits/s1prime_filtered/id_test.jsonl \
  --out experiments/002_e1prime_runs/rnastructure_id_test.jsonl \
  --workers 16 --skip-ambiguous --allow-partial \
  --cache experiments/002_e1prime_runs/rnastructure_id_test.partial.jsonl
echo "[predict id_test done] $(date -Iseconds)"
python scripts/score_external_baseline.py score \
  --split splits/s1prime_filtered/id_test.jsonl \
  --preds experiments/002_e1prime_runs/rnastructure_id_test.jsonl \
  --restrict-ids experiments/002_e1prime_runs/dump_p2fix_seed17/a_id_test.jsonl \
  --out experiments/002_e1prime_runs/rnastructure_id_test_score.json \
  --per-family-out experiments/002_e1prime_runs/rnastructure_id_test_by_family.json
echo "[all done] $(date -Iseconds)"
