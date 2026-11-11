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
PY="${RFF_PY:-/home/junxuanliu/.conda/envs/rff/bin/python}"

mkdir -p "${OUT}/shards" "${OUT}/tbl" "${OUT}/logs"
rm -f "${OUT}"/tbl/*.tbl "${OUT}"/tbl/*.done

"${PY}" "${ROOT}/scripts/shard_fasta.py" \
  "${IN}" "${OUT}/shards" --size "${SHARD}"

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
