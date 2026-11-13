#!/bin/bash
#SBATCH -J rff_ef_base
#SBATCH -p g-3090-2
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH -t 8:00:00
#SBATCH -o logs/ef_baseline_%j.out
# CPU-only (no --gres): EternaFold is a CPU folding engine.
#
# EternaFold (CONTRAfold-SE, EternaFoldParams.v1) on both evaluation splits,
# then the shared scorer with the same --restrict-ids dump and
# --min-true-pairs 4 / v3-label policy as every other baseline.
#
# Length cap 1200 nt (job 15619, no cap, stalled at 5,600/17,371 for 55 min on
# the 574 sequences above 2,000 nt; CONTRAfold is O(L^3)).  The cap removes
# 578/17,371 withheld rows and is declared in <out>.excluded.jsonl; every
# baseline is rescored on the resulting common cohort, never on a per-method
# subset.
set -euo pipefail
cd /home/junxuanliu/projects/RNA-Folding-Frontier
export PATH=/home/junxuanliu/.conda/envs/rff/bin:$PATH
echo "[start] $(date -Iseconds)"

python scripts/predict_eternafold.py \
  --split splits/s1prime_filtered/test_withheld.jsonl \
  --out experiments/002_e1prime_runs/eternafold_test_withheld.jsonl \
  --workers 16 --max-length 1200 --timeout 1800
echo "[predict withheld done] $(date -Iseconds)"

python scripts/predict_eternafold.py \
  --split splits/s1prime_filtered/id_test.jsonl \
  --out experiments/002_e1prime_runs/eternafold_id_test.jsonl \
  --workers 16 --max-length 1200 --timeout 1800
echo "[predict id_test done] $(date -Iseconds)"

for pair in "test_withheld:test_withheld" "id_test:id_test"; do
  split="${pair%%:*}"; tag="${pair##*:}"
  python scripts/score_external_baseline.py score \
    --split "splits/s1prime_filtered/${split}.jsonl" \
    --preds "experiments/002_e1prime_runs/eternafold_${tag}.jsonl" \
    --restrict-ids "experiments/002_e1prime_runs/dump_p2fix_seed17/a_${tag}.jsonl" \
    --out "experiments/002_e1prime_runs/eternafold_${tag}_score.json" \
    --per-family-out "experiments/002_e1prime_runs/eternafold_${tag}_by_family.json"
  echo "[score ${tag} done] $(date -Iseconds)"
done
echo "[all done] $(date -Iseconds)"