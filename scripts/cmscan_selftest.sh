#!/usr/bin/env bash
# End-to-end smoke test of the cmscan -> parse -> validate chain.
#
# Runs cmscan on a tiny, family-labelled sample so the parser and the
# prefix-agreement validator can be checked in minutes instead of after the
# full 94-shard scan.  Writes under data/processed/cmscan/selftest/.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="${ROOT}/data/processed/cmscan/selftest"
PY="${RFF_PY:-/home/junxuanliu/.conda/envs/rff/bin/python}"
CMS="${ROOT}/tools/infernal/bin/cmscan"
CM="${ROOT}/data/raw/rfam/Rfam.cm"

mkdir -p "${OUT}/shards" "${OUT}/tbl"

"${PY}" - "${ROOT}/data/processed/cmscan/input.fasta" "${OUT}/shards/shard_00000.fa" <<'PY'
import sys
from collections import Counter

src, dst = sys.argv[1], sys.argv[2]
want_per_prefix = 2
want_other = {"ts0": 6, "bprna": 6}
seen = Counter()
out = []
keep = False
for line in open(src):
    if line.startswith(">"):
        rid = line[1:].strip()
        source, _, body = rid.partition(":")
        if source == "archiveii":
            prefix = body.split("_", 1)[0]
            keep = seen[prefix] < want_per_prefix
            if keep:
                seen[prefix] += 1
        else:
            keep = seen[source] < want_other.get(source, 0)
            if keep:
                seen[source] += 1
    if keep:
        out.append(line)
open(dst, "w").writelines(out)
print(f"selected {sum(1 for l in out if l.startswith('>'))} sequences")
PY

"${CMS}" --rfam --cut_ga --nohmmonly --noali --cpu 2 \
  --tblout "${OUT}/tbl/shard_00000.tbl" -o /dev/null "${CM}" \
  "${OUT}/shards/shard_00000.fa"
touch "${OUT}/tbl/shard_00000.done"

"${PY}" "${ROOT}/scripts/parse_cmscan.py" \
  --tbl-dir "${OUT}/tbl" --out "${OUT}/assignments.tsv" --only-done

"${PY}" "${ROOT}/scripts/validate_cmscan.py" \
  --assignments "${OUT}/assignments.tsv" --out "${OUT}/validation"
