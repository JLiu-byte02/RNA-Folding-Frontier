#!/usr/bin/env bash
# Fallback when the CPU partition is fully allocated: run cmscan shards on the
# login node at low priority with bounded parallelism.
#
#   bash scripts/jobs/cmscan_local.sh <outdir> [nproc]
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
OUT="${1:?outdir}"
NPROC="${2:-12}"
FIRST="${3:-0}"
LAST="${4:-9999}"
CMS="${ROOT}/tools/infernal/bin/cmscan"
CM="${ROOT}/data/raw/rfam/Rfam.cm"

cd "${ROOT}"
ls "${OUT}"/shards/shard_*.fa |
  awk -v a="${FIRST}" -v b="${LAST}" '{n=substr($0, length($0)-7, 5)+0; if (n>=a && n<=b) print}' |
  xargs -P "${NPROC}" -I{} bash -c '
    shard="{}"
    base=$(basename "${shard}" .fa)
    tbl="'"${OUT}"'/tbl/${base}.tbl"
    if [ -s "${tbl}" ] && [ -f "${OUT}/tbl/${base}.done" ]; then exit 0; fi
    nice -n 19 '"${CMS}"' --rfam --cut_ga --nohmmonly --noali --cpu 1 \
      --tblout "${tbl}" -o /dev/null '"${CM}"' "${shard}" \
      && touch "'"${OUT}"'/tbl/${base}.done"
  '
echo "done: $(ls "${OUT}"/tbl/*.done | wc -l) shards"
