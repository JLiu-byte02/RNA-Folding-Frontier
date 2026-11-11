#!/usr/bin/env bash
# Relaxed-cmscan audit required by PREREGISTRATION Amendment 1.5.
#
# The frozen scan used --cut_ga, which is a curated per-model threshold: a hit
# below it is not reported at all, so the frozen output cannot answer "would a
# looser threshold have hit a test clan?". This re-scans the same shards with
# -E 10 and no GA cutoff and writes to a separate directory.
#
#   bash scripts/jobs/cmscan_relaxed.sh <outdir> [nproc]
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
OUT="${1:?outdir}"
NPROC="${2:-10}"
SHARDS="${RFF_CMS_SHARDS:-${ROOT}/data/processed/cmscan/shards}"
CMS="${ROOT}/tools/infernal/bin/cmscan"
CM="${ROOT}/data/raw/rfam/Rfam.cm"

cd "${ROOT}"
mkdir -p "${OUT}/tbl" "${OUT}/logs"
ls "${SHARDS}"/shard_*.fa |
  xargs -P "${NPROC}" -I{} bash -c '
    shard="{}"
    base=$(basename "${shard}" .fa)
    tbl="'"${OUT}"'/tbl/${base}.tbl"
    if [ -s "${tbl}" ] && [ -f "'"${OUT}"'/tbl/${base}.done" ]; then exit 0; fi
    nice -n 19 '"${CMS}"' --rfam -E 10 --nohmmonly --noali --cpu 1 \
      --tblout "${tbl}" -o /dev/null '"${CM}"' "${shard}" \
      && touch "'"${OUT}"'/tbl/${base}.done"
  '
echo "done: $(ls "${OUT}"/tbl/*.done | wc -l) shards"
