#!/usr/bin/env python3
"""Build a stand-in dump + split on local ArchiveII so the confirmatory chain can run.

PI Round 4g §3 asks for a rehearsal of `fit_energy_calibration` ->
greedy / A-energy-fixed -> `h24_decoder_compare` on real sequences, using
stand-in head scores. EternaFold or CONTRAfold probabilities are the preferred
stand-in; neither is installed here, so this uses ViennaRNA's partition function.
That is circular — ViennaRNA probabilities fed back into ViennaRNA — and the
report says so. The point is to exercise the *chain*, not to produce a result.

Outputs (under ``--out-dir``):

* ``dump.jsonl`` — the per-sequence dump schema D1/H2.4 read:
  ``{id, family, length, kept, truth, scored}`` with ``scored`` keyed ``"i,j"``;
* ``split.jsonl`` — ``{id, seq, family, pairs}``;
* ``standin_meta.json`` — provenance: archive hash, sample ids and seeds, and the
  circularity caveat.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import statistics
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

CLIP = 1e-6


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def logit(p: float) -> float:
    p = min(max(p, CLIP), 1.0 - CLIP)
    return math.log(p / (1.0 - p))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--archive", default="data/raw/archiveII/data.json")
    ap.add_argument("--limit", type=int, default=300)
    ap.add_argument("--max-length", type=int, default=400)
    ap.add_argument("--seed", type=int, default=17)
    ap.add_argument("--out-dir", default="experiments/010_rehearsal")
    args = ap.parse_args()

    import numpy as np
    import RNA

    archive = ROOT / args.archive
    data = json.loads(archive.read_text(encoding="utf-8"))
    items = [
        (rid, rec)
        for rid, rec in data.items()
        if rec.get("sequence") and len(rec["sequence"]) <= args.max_length
    ]
    rng = np.random.default_rng(args.seed)
    idx = rng.choice(len(items), size=min(args.limit, len(items)), replace=False)
    picked = [items[i] for i in sorted(idx)]

    out_dir = ROOT / args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    dump_rows = []
    split_rows = []
    for rid, rec in picked:
        seq = rec["sequence"]
        truth = sorted({(int(i), int(j)) for i, j in rec.get("structure", [])})
        fc = RNA.fold_compound(seq)
        fc.pf()
        bpp = fc.bpp()
        n = len(seq)
        scored = {}
        for i in range(1, n + 1):
            for j in range(i + 1, n + 1):
                p = float(bpp[i][j])
                if p > CLIP:
                    scored[f"{i - 1},{j - 1}"] = round(logit(p), 6)
        dump_rows.append(
            {
                "id": rid,
                "family": rid.split("_")[0],
                "length": n,
                "kept": [0, n],
                "truth": [list(p) for p in truth],
                "scored": scored,
            }
        )
        split_rows.append(
            {
                "id": rid,
                "seq": seq,
                "family": rid.split("_")[0],
                "pairs": [list(p) for p in truth],
            }
        )

    (out_dir / "dump.jsonl").write_text(
        "".join(json.dumps(r) + "\n" for r in dump_rows), encoding="utf-8"
    )
    (out_dir / "split.jsonl").write_text(
        "".join(json.dumps(r) + "\n" for r in split_rows), encoding="utf-8"
    )
    meta = {
        "purpose": "pipeline rehearsal only; no scientific claim",
        "archive": args.archive,
        "archive_sha256": sha256_file(archive),
        "n_sequences": len(dump_rows),
        "seed": args.seed,
        "max_length": args.max_length,
        "score_source": "ViennaRNA partition function base-pair probabilities, logit",
        "caveat": (
            "ViennaRNA probabilities fed back into ViennaRNA: circular by construction. "
            "EternaFold/CONTRAfold were preferred but are not installed locally."
        ),
        "median_length": statistics.median(r["length"] for r in dump_rows),
    }
    (out_dir / "standin_meta.json").write_text(
        json.dumps(meta, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({"out_dir": str(out_dir), "n": len(dump_rows)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
