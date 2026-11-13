#!/usr/bin/env bash
#SBATCH --job-name=rff_ec2gate
#SBATCH --partition=g-3090-2
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=/home/junxuanliu/projects/RNA-Folding-Frontier/logs/ec2_gate_%j.out
#
# Amendment 13 A13.2: confirmatory read of the FROZEN E-C2 coverage gate on the
# held-out splits.  Arms fixed(0.1875) / g1 / g2 are exactly the registered
# three.  A fourth exploratory call with --fixed-scale 0.09375 is run afterwards
# to separate the gate's information content from a pure constant-scale change
# (best_g1 is a constant, so g1 == fixed-at-0.09375 by construction); it is
# labelled exploratory and must not enter the confirmatory table.
# CPU-only: no --gres is requested.
set -euo pipefail

ROOT="/home/junxuanliu/projects/RNA-Folding-Frontier"
PY="/home/junxuanliu/.conda/envs/rff/bin/python"
cd "${ROOT}"
R=experiments/002_e1prime_runs

echo "[ec2gate] smoke on 300 dumped rows"
head -300 "${R}/dump_p2fix_seed17/a_test_withheld.jsonl" > "${R}/tmp_ec2_smoke_dump.jsonl"
"${PY}" scripts/score_ec2_gate.py \
  --dump "${R}/tmp_ec2_smoke_dump.jsonl" \
  --split splits/s1prime_filtered/test_withheld.jsonl \
  --pppl "${R}/ec1_p2_pppl.json" \
  --workers 4 --n-boot 200 --out "${R}/tmp_ec2_smoke.json"
"${PY}" - <<'EOF'
import json
d = json.load(open("experiments/002_e1prime_runs/tmp_ec2_smoke.json"))
n = d["n_items"]
print("[ec2gate] smoke n_items =", n)
assert n >= 250, n
EOF
rm -f "${R}/tmp_ec2_smoke_dump.jsonl" "${R}/tmp_ec2_smoke.json"

echo "[ec2gate] medium-scale timing self-check on 2000 dumped rows (O(n^2) guard)"
head -2000 "${R}/dump_p2fix_seed17/a_test_withheld.jsonl" > "${R}/tmp_ec2_med_dump.jsonl"
t0=$(date +%s)
"${PY}" scripts/score_ec2_gate.py \
  --dump "${R}/tmp_ec2_med_dump.jsonl" \
  --split splits/s1prime_filtered/test_withheld.jsonl \
  --pppl "${R}/ec1_p2_pppl.json" \
  --workers 8 --n-boot 200 --out "${R}/tmp_ec2_med.json" 2>&1 | tee "${R}/tmp_ec2_med.log"
t1=$(date +%s)
"${PY}" - <<'EOF'
import json
log = open("experiments/002_e1prime_runs/tmp_ec2_med.log").read()
d = json.load(open("experiments/002_e1prime_runs/tmp_ec2_med.json"))
n = d["n_items"]
assert n >= 1500, n
for arm in ("fixed", "g1", "g2"):
    line = f"[ec2] arm {arm}: {n} decodes ({n} items x 1)"
    assert line in log, f"missing/non-linear plan line for {arm}: {line!r}"
print("[ec2gate] medium n_items =", n, "decodes lines verified linear")
EOF
elapsed=$((t1 - t0))
echo "[ec2gate] medium-scale wall = ${elapsed}s"
if [ "${elapsed}" -gt 900 ]; then
  echo "[ec2gate] FATAL: 3 arms on the medium split took ${elapsed}s (>900s); refusing to run the full split"
  exit 1
fi
rm -f "${R}/tmp_ec2_med_dump.jsonl" "${R}/tmp_ec2_med.json" "${R}/tmp_ec2_med.log"

echo "[ec2gate] withheld confirmatory"
"${PY}" scripts/score_ec2_gate.py \
  --dump "${R}/dump_p2fix_seed17/a_test_withheld.jsonl" \
  --split splits/s1prime_filtered/test_withheld.jsonl \
  --pppl "${R}/ec1_p2_pppl.json" \
  --ref-preds "${R}/rnafold_test_withheld.jsonl" \
  --workers "${SLURM_CPUS_PER_TASK}" --n-boot 2000 \
  --out "${R}/ec2_test_withheld_score.json"

echo "[ec2gate] id_test confirmatory"
"${PY}" scripts/score_ec2_gate.py \
  --dump "${R}/dump_p2fix_seed17/a_id_test.jsonl" \
  --split splits/s1prime_filtered/id_test.jsonl \
  --pppl "${R}/ec1_p2_pppl.json" \
  --ref-preds "${R}/rnafold_id_test.jsonl" \
  --workers "${SLURM_CPUS_PER_TASK}" --n-boot 2000 \
  --out "${R}/ec2_id_test_score.json"

echo "[ec2gate] withheld EXPLORATORY scale-confound control (fixed-scale = best_g1 constant)"
"${PY}" scripts/score_ec2_gate.py \
  --dump "${R}/dump_p2fix_seed17/a_test_withheld.jsonl" \
  --split splits/s1prime_filtered/test_withheld.jsonl \
  --pppl "${R}/ec1_p2_pppl.json" \
  --fixed-scale 0.09375 \
  --workers "${SLURM_CPUS_PER_TASK}" --n-boot 2000 \
  --out "${R}/ec2_test_withheld_scalecontrol.json"

echo "[ec2gate] done"
