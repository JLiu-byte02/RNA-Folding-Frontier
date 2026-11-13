#!/bin/bash
#SBATCH -J rff_rs_base
#SBATCH -p g-3090-2
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH -t 8:00:00
#SBATCH -o logs/rs_baseline_%j.out
# CPU-only (no --gres): RNAstructure Fold is a CPU folding engine.
#
# Fold rejects IUPAC ambiguity codes outright (Invalid nucleobase), so 111/17,371
# withheld rows are dropped by declared rule (--skip-ambiguous) and logged to
# <out>.skipped.jsonl. --allow-partial keeps the ACGU rows when a long sequence
# times out; every failure lands in <out>.errors.jsonl so the exclusion is
# auditable rather than silent.
#
# RNAstructure 6.6 (Turner 2004 MFE) on both evaluation splits, then the
# shared scorer with the same --restrict-ids dump and --min-true-pairs 4 /
# v3-label policy as every other baseline. Fold needs the conda libstdc++
# (GLIBCXX_3.4.21) and DATAPATH; both are set inside the wrapper.
set -euo pipefail
cd /home/junxuanliu/projects/RNA-Folding-Frontier
export PATH=/home/junxuanliu/.conda/envs/rff/bin:$PATH
echo "[start] $(date -Iseconds)"

python scripts/predict_rnastructure.py \
  --split splits/s1prime_filtered/test_withheld.jsonl \
  --out experiments/002_e1prime_runs/rnastructure_test_withheld.jsonl \
  --workers 16 --skip-ambiguous --allow-partial \
  --cache experiments/002_e1prime_runs/rnastructure_test_withheld.partial.jsonl
echo "[predict withheld done] $(date -Iseconds)"

python scripts/predict_rnastructure.py \
  --split splits/s1prime_filtered/id_test.jsonl \
  --out experiments/002_e1prime_runs/rnastructure_id_test.jsonl \
  --workers 16 --skip-ambiguous --allow-partial \
  --cache experiments/002_e1prime_runs/rnastructure_id_test.partial.jsonl
echo "[predict id_test done] $(date -Iseconds)"

for pair in "test_withheld:test_withheld" "id_test:id_test"; do
  split="${pair%%:*}"; tag="${pair##*:}"
  python scripts/score_external_baseline.py score \
    --split "splits/s1prime_filtered/${split}.jsonl" \
    --preds "experiments/002_e1prime_runs/rnastructure_${tag}.jsonl" \
    --restrict-ids "experiments/002_e1prime_runs/dump_p2fix_seed17/a_${tag}.jsonl" \
    --out "experiments/002_e1prime_runs/rnastructure_${tag}_score.json" \
    --per-family-out "experiments/002_e1prime_runs/rnastructure_${tag}_by_family.json"
  echo "[score ${tag} done] $(date -Iseconds)"
done
echo "[all done] $(date -Iseconds)"