#!/usr/bin/env bash
#SBATCH --job-name=rff_4tool
#SBATCH --partition=g-3090-2
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=/home/junxuanliu/projects/RNA-Folding-Frontier/logs/rebuild_4tools_%j.out
#
# One-shot rebuild of the four-tool common cohort and every downstream number.
#
# Run only after experiments/013_mxfold2/mxfold2_<split>.jsonl exists. The
# script refuses to continue if a method is missing a prediction for any cohort
# id, so a silently smaller population cannot reach the table.
#
# The Windows workstation has no bash, so validate once with `bash -n` on the
# server before the first submission.
set -euo pipefail
ROOT="/home/junxuanliu/projects/RNA-Folding-Frontier"
PY="/home/junxuanliu/.conda/envs/rff/bin/python"
R="experiments/002_e1prime_runs"
MX="experiments/013_mxfold2"
cd "${ROOT}"

for split in test_withheld id_test; do
  echo "=== ${split}: building cohort"
  "${PY}" scripts/build_common_cohort.py \
    --split "splits/s1prime_filtered/${split}.jsonl" \
    --dump "${R}/dump_p2fix_seed17/a_${split}.jsonl" \
    --preds "${R}/rnafold_${split}.jsonl" \
            "${R}/eternafold_${split}.jsonl" \
            "${R}/rnastructure_${split}.jsonl" \
            "${MX}/mxfold2_${split}.jsonl" \
    --tool-names rnafold eternafold rnastructure mxfold2 \
    --out "${R}/cohort/${split}.jsonl" \
    --report "${R}/cohort/${split}.report.json" \
    --excluded-out "${R}/cohort/${split}.excluded.jsonl"

  echo "=== ${split}: scoring every method on the frozen cohort"
  declare -A PREDS=(
    [rnafold]="${R}/rnafold_${split}.jsonl"
    [eternafold]="${R}/eternafold_${split}.jsonl"
    [rnastructure]="${R}/rnastructure_${split}.jsonl"
    [mxfold2]="${MX}/mxfold2_${split}.jsonl"
    [head17]="${R}/headpreds/seed17_${split}.jsonl"
    [head23]="${R}/headpreds/seed23_${split}.jsonl"
    [head41]="${R}/headpreds/seed41_${split}.jsonl"
  )
  for name in rnafold eternafold rnastructure mxfold2 head17 head23 head41; do
    "${PY}" scripts/score_external_baseline.py score \
      --split "splits/s1prime_filtered/${split}.jsonl" \
      --preds "${PREDS[$name]}" \
      --restrict-ids "${R}/cohort/${split}.jsonl" \
      --min-true-pairs 4 \
      --tool "${name}" \
      --out "${R}/${name}_${split}_score.json" \
      --per-family-out "${R}/${name}_${split}_by_family.json"
  done

  echo "=== ${split}: paired family bootstrap"
  "${PY}" scripts/cohort_paired_eval.py \
    --split "splits/s1prime_filtered/${split}.jsonl" \
    --cohort "${R}/cohort/${split}.jsonl" \
    --method RNAfold="${R}/rnafold_${split}.jsonl" \
    --method EternaFold="${R}/eternafold_${split}.jsonl" \
    --method RNAstructure="${R}/rnastructure_${split}.jsonl" \
    --method MXfold2="${MX}/mxfold2_${split}.jsonl" \
    --method head17="${R}/headpreds/seed17_${split}.jsonl" \
    --method head23="${R}/headpreds/seed23_${split}.jsonl" \
    --method head41="${R}/headpreds/seed41_${split}.jsonl" \
    --mean-of head_mean=head17,head23,head41 \
    --compare head_mean:RNAfold \
    --compare head_mean:EternaFold \
    --compare head_mean:RNAstructure \
    --compare head_mean:MXfold2 \
    --compare head17:RNAfold \
    --compare head17:EternaFold \
    --out "${R}/cohort/paired_${split}.json"

  "${PY}" - "${split}" <<'PY'
import json, sys
split = sys.argv[1]
base = "experiments/002_e1prime_runs/cohort"
report = json.load(open(f"{base}/{split}.report.json"))
paired = json.load(open(f"{base}/paired_{split}.json"))
missing = paired["excluded"]["missing_prediction"]
assert missing == 0, f"{split}: {missing} cohort rows lack a prediction"
missing_by_tool = {
    name: row.get("n_missing_from_label_eligible", 0)
    for name, row in report.get("per_tool", {}).items()
}
methods = paired["methods"]
gap = max(abs(m["family_macro_f1"] - m["micro_f1"]) for m in methods.values())
print(json.dumps({
    "split": split,
    "n_sequences": paired["n_sequences"],
    "n_families": paired["n_families"],
    "max_macro_minus_micro": round(gap, 4),
    "missing_from_label_eligible_by_tool": missing_by_tool,
    "methods": {k: round(v["family_macro_f1"], 4) for k, v in methods.items()},
}, indent=1))
assert gap <= 0.15, f"{split}: macro/micro gap {gap:.3f} exceeds the .15 invariant"
PY
done

echo "[rebuild done] $(date -Iseconds)"
