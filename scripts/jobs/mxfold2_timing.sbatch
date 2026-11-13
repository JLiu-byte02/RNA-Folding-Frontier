#!/usr/bin/env bash
#SBATCH --job-name=rff_mx_timing
#SBATCH --partition=g-3090-2
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --gres=gpu:1
#SBATCH --exclude=gpu0004
#SBATCH --output=/home/junxuanliu/projects/RNA-Folding-Frontier/logs/mx_timing_%j.out
#
# Purpose: measure MXfold2 per-sequence training cost at the real 600-nt cap on a
# stride sample that spans the whole train split.  The smoke test used max_len 120,
# which is not representative, and MXfold2 trains with batch_size=1 so the wall
# time scales with sum(L^3) over sequences.  The resulting JSON decides how many
# epochs fit in one 24 h job; any reduction is recorded as an explicit deviation.
set -euo pipefail

ROOT="/home/junxuanliu/projects/RNA-Folding-Frontier"
PY="/home/junxuanliu/.conda/envs/rff/bin/python"
MXP="/home/junxuanliu/.conda/envs/mxfold2/bin/python"
export PYTHONPATH="/home/junxuanliu/projects/vendor/MXfold2"
cd "${ROOT}"

OUT="experiments/013_mxfold2/timing"
mkdir -p "${OUT}" logs
nvidia-smi --query-gpu=name,memory.total,memory.used --format=csv,noheader

if [ ! -f "${OUT}/data/train.lst" ]; then
  "${PY}" scripts/export_mxfold2_bpseq.py \
    --split splits/s1prime_filtered/train.jsonl \
    --out-dir "${OUT}/data" --max-len 600 --stride 40 --limit 400 \
    --manifest "${OUT}/export_manifest.json"
fi

start=$(date +%s)
"${MXP}" -m mxfold2 train \
  --model MixC --epochs 1 --gpu 0 --disable-progress-bar \
  --param "${OUT}/model.pth" --save-config "${OUT}/model.conf" \
  "${OUT}/data/train.lst"
end=$(date +%s)

"${PY}" - <<PY
import json, pathlib
out = pathlib.Path("${OUT}")
man = json.loads((out / "export_manifest.json").read_text())
wall = ${end} - ${start}
n = man["n_exported"]
per_seq = wall / max(n, 1)
rec = {
    "job_id": "${SLURM_JOB_ID}",
    "gpu": "RTX 3090",
    "max_len": man["max_len"],
    "stride": man["stride"],
    "n_train_sample": n,
    "exported_len_mean": man["exported_len_mean"],
    "exported_len_p50": man["exported_len_p50"],
    "exported_len_max": man["exported_len_max"],
    "wall_seconds_one_epoch": wall,
    "seconds_per_sequence_one_epoch": round(per_seq, 4),
    "train_rows_full": man["n_rows_in_split"],
    "full_export_estimate": int(round(n * man["stride"])),
    "full_epoch_seconds_estimate": int(round(per_seq * n * man["stride"])),
    "note": "sample is a stride subset; per-sequence cost is length-dependent, "
            "so the full-split estimate assumes the sample length distribution "
            "matches the split. Check exported_len_mean of the full export.",
}
(out / "timing.json").write_text(json.dumps(rec, indent=2, sort_keys=True) + "\n")
print(json.dumps(rec, indent=2, sort_keys=True))
PY
