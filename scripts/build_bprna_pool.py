#!/usr/bin/env python3
"""Convert bpRNA-1m (HuggingFace mirror) into the project's record format.

The authors' ``bpRNA_1m_90.zip`` is 4.02 GB and serves at ~8 KB/s, which makes it
unusable. The ``rouskinlab/bpRNA-1m`` mirror carries the same 66,715 sequences
that survived their conversion of the original 102,318 entries, with base pairs.
We therefore re-derive a 90%-identity subset ourselves with MMseqs2 and record
that this is a re-derivation, not the authors' file.

Outputs
  data/processed/bprna/pool.jsonl      one record per line: id, seq, pairs
  data/processed/bprna/pool.fasta      sequences for clustering
  data/processed/bprna/convert.json    provenance + counts
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", default="data/raw/bprna/bpRNA_1m.json")
    ap.add_argument("--out", default="data/processed/bprna")
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    data = json.loads(Path(args.src).read_text())
    n_in = len(data)
    n_written = 0
    n_bad = 0
    lengths = []

    with open(out / "pool.jsonl", "w") as jf, open(out / "pool.fasta", "w") as ff:
        for rid, payload in data.items():
            seq = str(payload.get("sequence", "")).upper().replace("T", "U")
            seq = re.sub(r"[^ACGU]", "", seq)
            if not seq:
                n_bad += 1
                continue
            pairs = [
                [int(i), int(j)]
                for i, j in payload.get("structure", [])
                if 0 <= int(i) < len(seq) and 0 <= int(j) < len(seq)
            ]
            jf.write(json.dumps({"id": rid, "seq": seq, "pairs": pairs}) + "\n")
            ff.write(f">{rid}\n{seq}\n")
            lengths.append(len(seq))
            n_written += 1

    lengths.sort()
    meta = {
        "source": str(args.src),
        "n_input": n_in,
        "n_written": n_written,
        "n_dropped_empty": n_bad,
        "length_min": lengths[0],
        "length_median": lengths[len(lengths) // 2],
        "length_max": lengths[-1],
        "note": (
            "Re-derived from the HuggingFace mirror of bpRNA-1m; not the authors' "
            "bpRNA_1m_90.zip. 90% identity subset is produced by MMseqs2 "
            "easy-cluster (see scripts/cluster_pool.sh)."
        ),
    }
    (out / "convert.json").write_text(json.dumps(meta, indent=2))
    print(json.dumps(meta, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
