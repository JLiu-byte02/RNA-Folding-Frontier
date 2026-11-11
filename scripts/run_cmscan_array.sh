#!/usr/bin/env bash
# Shard a FASTA and run cmscan against Rfam 15.0 CMs as a Slurm array.
#
#   bash scripts/run_cmscan_array.sh <input.fasta> <outdir> [shard_size] [jobs]
#
# Each shard writes <outdir>/shard_XXXX.tbl (Infernal tabular output) and a
# .done marker. The CM database must already be pressed (scripts/server_download_rfam.sh).
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
IN="${1:?input fasta}"
OUT="${2:?output dir}"
SHARD="${3:-500}"
JOBS="${4:-16}"
CM="${ROOT}/data/raw/rfam/Rfam.cm"
CMS="${ROOT}/tools/infernal/bin/cmscan"

mkdir -p "${OUT}/shards" "${OUT}/tbl" "${OUT}/logs"

/home/junxuanliu/.conda/envs/rff/bin/python - "$IN" "${OUT}/shards" "$SHARD" <<'PY'
import sys, pathlib
src, outdir, shard = sys.argv[1], pathlib.Path(sys.argv[2]), int(sys.argv[3])
outdir.mkdir(parents=True, exist_ok=True)
name, buf, idx = None, [], 0
def flush(idx, name, buf):
    if name is None:
        return
    with open(outdir / f"shard_{idx:05d}.fa", "w") as fh:
        fh.write(f">{name}\n{''.join(buf)}\n")
n = 0
with open(src) as fh:
    for line in fh:
        if line.startswith(">"):
            flush(idx, name, buf)
            if name is not None:
                n += 1
            name, buf = line[1:].strip().split()[0], []
            idx = n // shard
        else:
            buf.append(line.strip())
flush(idx, name, buf)
print("shards:", idx + 1)
PY

N=$(ls "${OUT}/shards"/shard_*.fa | wc -l)
echo "submitting ${N} shards, ${JOBS} at a time"

sbatch --partition=c-80-255 --cpus-per-task=2 --mem=8G --time=24:00:00 \
  --job-name=rff_cmscan --array=0-$((N-1))%"${JOBS}" \
  --output="${OUT}/logs/cmscan_%A_%a.out" \
  --wrap="i=\$(printf '%05d' \$SLURM_ARRAY_TASK_ID); \
    ${CMS} --rfam --cut_ga --nohmmonly --noali --cpu 2 \
      --tblout ${OUT}/tbl/shard_\${i}.tbl \
      -o /dev/null ${CM} ${OUT}/shards/shard_\${i}.fa \
    && touch ${OUT}/tbl/shard_\${i}.done"
