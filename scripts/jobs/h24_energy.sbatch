#!/usr/bin/env bash
#SBATCH --job-name=rff_h24
#SBATCH --partition=g-3090-2
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=/home/junxuanliu/projects/RNA-Folding-Frontier/logs/h24_%j.out
#
# H2.4 with A-energy-fixed using the frozen (s,b) fit. Exploratory on id_test
# and withheld; test_novel remains embargoed.
set -euo pipefail

ROOT="/home/junxuanliu/projects/RNA-Folding-Frontier"
PY="/home/junxuanliu/.conda/envs/rff/bin/python"
FIT="experiments/002_e1prime_runs/energy_fit_coarse_p2fix.json"
DUMP="experiments/002_e1prime_runs/dump_p2fix_seed17"

cd "${ROOT}"
S=$("${PY}" -c "import json; print(json.load(open('${FIT}'))['selected']['s'])")
B=$("${PY}" -c "import json; print(json.load(open('${FIT}'))['selected']['b'])")
echo "frozen s=${S} b=${B}"

for SPLIT in id_test test_withheld; do
  "${PY}" scripts/h24_decoder_compare.py \
    --dump "${DUMP}/a_${SPLIT}.jsonl" \
    --split "splits/s1prime_filtered/${SPLIT}.jsonl" \
    --name "a_${SPLIT}" \
    --energy vienna --energy-scale "${S}" --energy-offset "${B}" \
    --label-version v3 --min-true-pairs 4 \
    --out "experiments/002_e1prime_runs/h24_seed17_${SPLIT}_v3primary.json"
done
