#!/usr/bin/env bash
#SBATCH --job-name=rff_mx_train
#SBATCH --partition=g-3090-2
#SBATCH --cpus-per-task=12
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --gres=gpu:1
#SBATCH --exclude=gpu0004
#SBATCH --output=/home/junxuanliu/projects/RNA-Folding-Frontier/logs/mx_train_%j.out
#
# Amendment 13 A13.1 secondary control: retrain MXfold2 (MixC, README defaults)
# on OUR train split so the comparison does not mix training data.  Sequences
# longer than --max-len are dropped because MXfold2's O(L^2) nested DP and
# memory make them impractical; the drop rate is recorded in the manifest.
set -euo pipefail

ROOT="/home/junxuanliu/projects/RNA-Folding-Frontier"
PY="/home/junxuanliu/.conda/envs/rff/bin/python"
MXP="/home/junxuanliu/.conda/envs/mxfold2/bin/python"
export PYTHONPATH="/home/junxuanliu/projects/vendor/MXfold2"
cd "${ROOT}"

OUT="experiments/013_mxfold2/run"
mkdir -p "${OUT}" logs
MAXLEN="${MAXLEN:-600}"

nvidia-smi --query-gpu=name,memory.total,memory.used --format=csv,noheader

if [ ! -f "${OUT}/data/train.lst" ]; then
  "${PY}" scripts/export_mxfold2_bpseq.py \
    --split splits/s1prime_filtered/train.jsonl \
    --out-dir "${OUT}/data" --max-len "${MAXLEN}" \
    --manifest "${OUT}/export_manifest.json"
fi

git -C "${ROOT}" rev-parse HEAD > "${OUT}/repo_head.txt"
git -C /home/junxuanliu/projects/vendor/MXfold2 rev-parse HEAD > "${OUT}/mxfold2_commit.txt"

start=$(date +%s)
"${MXP}" -m mxfold2 train \
  --model MixC --epochs 10 --gpu 0 --disable-progress-bar \
  --param "${OUT}/model.pth" --save-config "${OUT}/model.conf" \
  --log-dir "${OUT}/ckpt" \
  "${OUT}/data/train.lst"
end=$(date +%s)

"${PY}" - <<PY
import hashlib, json, pathlib, subprocess, os
out = pathlib.Path("${OUT}")
def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()
man = {
    "mxfold2_commit": "51b2136",
    "repo_head": (out / "repo_head.txt").read_text().strip(),
    "model": "MixC",
    "epochs": 10,
    "gpu": 0,
    "batch_size": 1,
    "max_len": int("${MAXLEN}"),
    "train_wall_seconds": ${end} - ${start},
    "param_sha256": sha(out / "model.pth"),
    "config_sha256": sha(out / "model.conf"),
    "lst_sha256": sha(out / "data/train.lst"),
    "export_manifest": json.loads((out / "export_manifest.json").read_text()),
}
(out / "train_manifest.json").write_text(json.dumps(man, indent=2, sort_keys=True) + "\n")
print(json.dumps(man, indent=2, sort_keys=True))
PY