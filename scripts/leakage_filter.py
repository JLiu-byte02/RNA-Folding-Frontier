#!/usr/bin/env python3
"""Screen candidate S1' splits for residual train/test homology.

Two independent screens, both applied before the split is frozen:

``family_homology``
    A training sequence may be assigned to family F by its best cmscan hit yet
    still hit a *test* family at E <= ``--evalue-max`` with a weaker hit.  Such
    a sequence is weak but real homology to the held-out family, so it is
    removed from training.  This is the "reverse cmscan" screen; it needs no
    extra alignment run because the full scan already reports every hit.

``kmers``
    Exact 15-mer containment of each test sequence against all training
    sequences.  Also reports exact duplicate sequences.  K-mers whose document
    frequency exceeds ``--max-train-df`` are dropped as uninformative
    (poly-A/poly-G and other low-complexity kmers otherwise dominate the
    statistic).

The output directory mirrors the input splits with filtered files plus
``filter_report.json``, which records every dropped id and the containment
distribution of the withheld clans so the threshold can be calibrated on data
rather than guessed.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

SPLITS = ("train", "id_test", "test_novel", "test_withheld")


def read_jsonl(path: Path) -> list[dict]:
    if not path.exists():
        return []
    with open(path) as fh:
        return [json.loads(line) for line in fh if line.strip()]


def write_jsonl(path: Path, records: list[dict]) -> str:
    with open(path, "w") as fh:
        for rec in records:
            fh.write(json.dumps(rec) + "\n")
    h = hashlib.sha256()
    for rec in sorted(records, key=lambda r: r["id"]):
        h.update(rec["id"].encode())
    return h.hexdigest()


def seq_hash(seq: str) -> str:
    return hashlib.blake2b(seq.encode(), digest_size=16).hexdigest()


def kmer_codes(seq: str, k: int) -> np.ndarray:
    """Unique k-mer set packed as base-4 integers (injective up to k=32)."""
    if len(seq) < k:
        return np.empty(0, dtype=np.uint64)
    lut = np.full(256, 255, dtype=np.uint8)
    for code, base in enumerate("ACGU"):
        lut[ord(base)] = code
        lut[ord(base.lower())] = code
    arr = np.frombuffer(seq.encode(), dtype=np.uint8)
    codes = lut[arr]
    valid = codes != 255
    if arr.size >= k:
        win = np.lib.stride_tricks.sliding_window_view(codes, k)
        ok = np.lib.stride_tricks.sliding_window_view(valid, k).all(axis=1)
        win = win[ok]
    else:
        return np.empty(0, dtype=np.uint64)
    if win.shape[0] == 0:
        return np.empty(0, dtype=np.uint64)
    packed = np.zeros(win.shape[0], dtype=np.uint64)
    for j in range(k):
        packed = (packed << np.uint64(2)) | win[:, j].astype(np.uint64)
    return np.unique(packed)


def build_index(train_seqs: list[dict], k: int, max_df: float):
    n_train = max(1, len(train_seqs))
    df: Counter = Counter()
    for rec in train_seqs:
        df.update(kmer_codes(rec["seq"], k).tolist())
    cut = max(5, int(round(max_df * n_train)))
    codes = np.array(sorted(c for c, n in df.items() if n <= cut), dtype=np.uint64)
    return np.unique(codes), df, cut


def _present(codes: np.ndarray, sorted_index: np.ndarray) -> np.ndarray:
    """Membership in a sorted uint64 index without rebuilding it per query."""
    if sorted_index.size == 0 or codes.size == 0:
        return np.zeros(codes.shape, dtype=bool)
    pos = np.searchsorted(sorted_index, codes)
    pos = np.clip(pos, 0, sorted_index.size - 1)
    return sorted_index[pos] == codes


def containment(seq: str, k: int, index, df: Counter, cut: int) -> tuple[float, float, int]:
    kmers = kmer_codes(seq, k)
    if kmers.size == 0:
        return 0.0, 0.0, 0
    raw = float(_present(kmers, index).sum()) / kmers.size
    info = kmers[[df.get(int(c), 0) <= cut for c in kmers]]
    if info.size == 0:
        return raw, 0.0, 0
    return raw, float(_present(info, index).sum()) / info.size, int(info.size)


def load_hits_by_query(path: Path, evalue_max: float) -> dict[str, set[str]]:
    out: dict[str, set[str]] = defaultdict(set)
    if not path.exists():
        return out
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != len(header):
                continue
            try:
                ev = float(parts[idx["evalue"]])
            except ValueError:
                continue
            if ev <= evalue_max:
                out[parts[idx["query"]]].add(parts[idx["family_accession"]])
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--split-dir", default="splits/s1prime")
    ap.add_argument("--all-hits", default="data/processed/cmscan/all_hits.tsv")
    ap.add_argument("--out", default="splits/s1prime_filtered")
    ap.add_argument("--evalue-max", type=float, default=1.0)
    ap.add_argument("--k", type=int, default=15)
    ap.add_argument("--containment-max", type=float, default=0.90)
    ap.add_argument("--max-train-df", type=float, default=0.05)
    ap.add_argument(
        "--k-scan",
        default="8,12,15,20",
        help="k values for the containment curve fitted on withheld clans; empty to skip",
    )
    args = ap.parse_args()

    split_dir = Path(args.split_dir)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    splits = {name: read_jsonl(split_dir / f"{name}.jsonl") for name in SPLITS}

    # ---- screen 1: homology to a held-out family -------------------------
    test_fams: set[str] = set()
    for name in ("test_novel", "test_withheld"):
        test_fams |= {rec["family"] for rec in splits[name]}
    hits = load_hits_by_query(Path(args.all_hits), args.evalue_max)
    dropped_homology: dict[str, list[str]] = {}
    for name in ("train", "id_test"):
        keep, drop = [], []
        for rec in splits[name]:
            fams = hits.get(rec["id"], set())
            bad = fams & test_fams
            if bad:
                drop.append({"id": rec["id"], "hits": sorted(bad)})
            else:
                keep.append(rec)
        splits[name] = keep
        dropped_homology[name] = drop

    # ---- screen 2: exact duplicates and k-mer containment ----------------
    train_seqs = splits["train"]
    train_hashes = Counter(seq_hash(rec["seq"]) for rec in train_seqs)
    n_train = max(1, len(train_seqs))

    train_kmers, df, df_keep_cut = build_index(train_seqs, args.k, args.max_train_df)

    containment_stats: dict[str, dict] = {}
    dropped_dup: dict[str, list[str]] = {}
    for name in ("id_test", "test_novel", "test_withheld"):
        keep, drop = [], []
        vals = []
        n_unscreened = 0
        for rec in splits[name]:
            if train_hashes.get(seq_hash(rec["seq"])):
                drop.append({"id": rec["id"], "reason": "exact_duplicate"})
                continue
            raw, info, n_info = containment(
                rec["seq"], args.k, train_kmers, df, df_keep_cut
            )
            vals.append({"id": rec["id"], "clan": rec["clan"], "raw": raw, "info": info, "n_info": n_info})
            if n_info < 20:
                n_unscreened += 1
            if n_info >= 20 and info > args.containment_max:
                drop.append(
                    {"id": rec["id"], "reason": "kmer_containment", "info": round(info, 4)}
                )
                continue
            keep.append(rec)
        splits[name] = keep
        dropped_dup[name] = drop
        by_clan: dict[str, list[float]] = defaultdict(list)
        for v in vals:
            by_clan[v["clan"]].append(v["info"])
        containment_stats[name] = {
            "n": len(vals),
            "quantiles": {
                str(q): round(float(np.quantile([v["info"] for v in vals], q)), 4)
                for q in (0.5, 0.9, 0.99, 1.0)
            }
            if vals
            else {},
            "n_above_max": sum(1 for v in vals if v["info"] > args.containment_max),
            "n_unscreened_short": n_unscreened,
            "per_clan_withheld": {
                clan: {
                    "n": len(v),
                    "max": round(max(v), 4),
                    "median": round(float(np.median(v)), 4),
                }
                for clan, v in sorted(by_clan.items())
            }
            if name == "test_withheld"
            else {},
        }

    k_curve: dict[str, dict] = {}
    scan = [int(v) for v in args.k_scan.split(",") if v.strip()]
    withheld = splits["test_withheld"]
    for k in scan:
        if k == args.k:
            idx, dfk, cutk = train_kmers, df, df_keep_cut
        else:
            idx, dfk, cutk = build_index(train_seqs, k, args.max_train_df)
        per_clan: dict[str, list[float]] = defaultdict(list)
        for rec in withheld:
            _raw, info, n_info = containment(rec["seq"], k, idx, dfk, cutk)
            if n_info >= 20:
                per_clan[rec["clan"]].append(info)
        allv = [v for vals in per_clan.values() for v in vals]
        k_curve[str(k)] = {
            "n": len(allv),
            "median": round(float(np.median(allv)), 4) if allv else None,
            "p99": round(float(np.quantile(allv, 0.99)), 4) if allv else None,
            "max": round(max(allv), 4) if allv else None,
            "frac_above_0.9": round(float(np.mean([v > 0.9 for v in allv])), 4)
            if allv
            else None,
        }

    hashes = {name: write_jsonl(out / f"{name}.jsonl", splits[name]) for name in SPLITS}
    report = {
        "params": {
            "evalue_max": args.evalue_max,
            "k": args.k,
            "containment_max": args.containment_max,
            "max_train_df": args.max_train_df,
            "df_cut_count": df_keep_cut,
            "train_kmers": int(train_kmers.size),
        },
        "counts": {name: len(splits[name]) for name in SPLITS},
        "dropped_family_homology": {
            k: {"n": len(v), "examples": v[:20]} for k, v in dropped_homology.items()
        },
        "dropped_duplicate_or_containment": {
            k: {"n": len(v), "by_reason": dict(Counter(d["reason"] for d in v)), "examples": v[:20]}
            for k, v in dropped_dup.items()
        },
        "containment": containment_stats,
        "k_curve_on_withheld": k_curve,
        "hashes": hashes,
    }
    (out / "filter_report.json").write_text(json.dumps(report, indent=2))
    print(json.dumps({k: report[k] for k in ("counts",)}, indent=2))
    print(
        json.dumps(
            {
                "dropped_homology": {
                    k: v["n"] for k, v in report["dropped_family_homology"].items()
                },
                "dropped_dup_containment": {
                    k: v["by_reason"]
                    for k, v in report["dropped_duplicate_or_containment"].items()
                },
                "containment_quantiles": {
                    k: v["quantiles"] for k, v in containment_stats.items()
                },
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
