#!/usr/bin/env python3
"""Provenance audit for the HuggingFace bpRNA-1m mirror.

The mirror is only usable as a training pool if the splits we already trust
(TR0/VL0/TS0) are contained in it with identical labels. This script measures
that containment, tabulates the composition of the mirror, and records the
non-canonical-pair rate so the pair policy can be stated explicitly.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rff.data.io import load_path  # noqa: E402

CANONICAL = {("A", "U"), ("U", "A"), ("G", "C"), ("C", "G"), ("G", "U"), ("U", "G")}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pool", default="data/processed/bprna/pool.jsonl")
    ap.add_argument("--raw", default="data/raw/bprna/bpRNA_1m.json")
    ap.add_argument("--out", default="experiments/000_e0_leakage/bprna_mirror")
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    by_seq: dict[str, list[tuple[str, list[list[int]]]]] = {}
    prefixes: Counter[str] = Counter()
    n_pairs = 0
    n_noncanon = 0
    n_gu = 0
    lengths = []

    with open(ROOT / args.pool) as fh:
        for line in fh:
            rec = json.loads(line)
            seq = rec["seq"]
            pairs = rec["pairs"]
            by_seq.setdefault(seq, []).append((rec["id"], pairs))
            parts = rec["id"].split("_")
            prefixes["_".join(parts[:2])] += 1
            lengths.append(len(seq))
            for i, j in pairs:
                a, b = seq[i], seq[j]
                n_pairs += 1
                if (a, b) == ("G", "U") or (a, b) == ("U", "G"):
                    n_gu += 1
                if (a, b) not in CANONICAL:
                    n_noncanon += 1

    lengths.sort()
    report: dict = {
        "pool_records": sum(prefixes.values()),
        "unique_sequences": len(by_seq),
        "duplicate_sequences": sum(prefixes.values()) - len(by_seq),
        "prefixes": prefixes.most_common(12),
        "length": {
            "min": lengths[0],
            "median": lengths[len(lengths) // 2],
            "p99": lengths[int(0.99 * len(lengths))],
            "max": lengths[-1],
        },
        "pairs": {
            "total": n_pairs,
            "gu_wobble": n_gu,
            "non_canonical": n_noncanon,
            "non_canonical_rate": round(n_noncanon / n_pairs, 5) if n_pairs else 0.0,
        },
        "containment": {},
    }

    for name, sub in (
        ("TR0", "data/raw/bpfold/BPfold_data/bpRNA/TR0"),
        ("VL0", "data/raw/bpfold/BPfold_data/bpRNA/VL0"),
        ("TS0", "data/raw/bpfold/BPfold_data/bpRNA/TS0"),
    ):
        recs = load_path(ROOT / sub)
        found = 0
        identical = 0
        superset = 0
        identical_canonical = 0
        extra_pairs = 0
        extra_noncanon = 0
        missing_pairs = 0
        for rec in recs:
            cands = by_seq.get(rec.seq)
            if not cands:
                continue
            found += 1
            want = {tuple(sorted(p)) for p in rec.pairs}
            best = max(
                ({tuple(sorted(p)) for p in pr} for _, pr in cands),
                key=lambda s: len(want & s),
            )
            if best == want:
                identical += 1
            best_canon = {
                (i, j)
                for i, j in best
                if (rec.seq[i], rec.seq[j]) in CANONICAL
            }
            if best_canon == want:
                identical_canonical += 1
            if want <= best:
                superset += 1
            missing_pairs += len(want - best)
            for i, j in best - want:
                extra_pairs += 1
                if (rec.seq[i], rec.seq[j]) not in CANONICAL:
                    extra_noncanon += 1
        report["containment"][name] = {
            "n": len(recs),
            "found_in_mirror": found,
            "found_rate": round(found / len(recs), 4),
            "identical_pairs": identical,
            "identical_rate_of_found": round(identical / found, 4) if found else 0.0,
            "identical_after_canonical_filter": identical_canonical,
            "identical_after_canonical_filter_rate": (
                round(identical_canonical / found, 4) if found else 0.0
            ),
            "bpseq_is_subset_of_mirror": superset,
            "subset_rate": round(superset / found, 4) if found else 0.0,
            "extra_mirror_pairs": extra_pairs,
            "extra_mirror_pairs_noncanonical": extra_noncanon,
            "extra_noncanonical_rate": (
                round(extra_noncanon / extra_pairs, 4) if extra_pairs else 0.0
            ),
            "missing_bpseq_pairs": missing_pairs,
        }

    raw = ROOT / args.raw
    if raw.exists():
        report["file"] = {
            "path": str(args.raw),
            "sha256": sha256_file(raw),
            "bytes": raw.stat().st_size,
        }
    report["hf_revision"] = "8b973adeba0218113a696dd1c22ada682f26324a"
    report["hf_source"] = "https://huggingface.co/datasets/rouskinlab/bpRNA-1m"

    (out / "summary.json").write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
