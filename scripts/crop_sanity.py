#!/usr/bin/env python3
"""Positional sanity check for the crop/full-context comparison (E1').

The E1' mechanism claim is that a +/-k crop removes *distal context* and nothing
else. That is only true if the frozen model's positional encoding is relative:
otherwise a crop also moves the residue's absolute position, and any measured
effect is confounded.

Two measurements, both on real sequences:

``agreement``
    cosine(crop_k(i), full(i)) versus cosine(crop_k(i), full(j != i)). The
    second number is the floor; the gap is how much position-specific
    information survives the crop.
``translation``
    take a (2k+1) window, embed it alone, then embed the same window inside 100
    nt of unrelated flank and read the same residue. Cosine ~1 means the model
    is translation invariant (relative positional encoding); a visible drop means
    absolute positions leak into the embedding.
"""

from __future__ import annotations

import argparse
import json
import random
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import torch  # noqa: E402
import torch.nn.functional as F  # noqa: E402

from rff.models.backbone import load_backbone  # noqa: E402


def read_fasta(path: Path, limit: int) -> list[tuple[str, str]]:
    out: list[tuple[str, str]] = []
    name, buf = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    out.append((name, "".join(buf)))
                    if len(out) >= limit:
                        return out
                name, buf = line[1:].strip().split()[0], []
            else:
                buf.append(line.strip())
    if name and len(out) < limit:
        out.append((name, "".join(buf)))
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--backbone", default="rinalmo-giga")
    ap.add_argument("--path", default="data/models/rinalmo-giga")
    ap.add_argument("--fasta", default="data/processed/cmscan/selftest/shards/shard_00000.fa")
    ap.add_argument("--n-seqs", type=int, default=6)
    ap.add_argument("--min-len", type=int, default=180)
    ap.add_argument("--max-len", type=int, default=400)
    ap.add_argument("--ks", default="8,32,128,256")
    ap.add_argument("--n-positions", type=int, default=24)
    ap.add_argument("--flank", type=int, default=100)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    rng = random.Random(17)
    seqs = [
        (name, seq)
        for name, seq in read_fasta(Path(args.fasta), 4000)
        if args.min_len <= len(seq) <= args.max_len
    ][: args.n_seqs]
    if not seqs:
        print("no sequences in the length window")
        return 1

    backbone = load_backbone(args.backbone, path=args.path)
    ks = [int(v) for v in args.ks.split(",") if v.strip()]
    report: dict = {
        "backbone": args.backbone,
        "device": torch.cuda.get_device_name(0) if torch.cuda.is_available() else "cpu",
        "n_seqs": len(seqs),
        "lengths": [len(s) for _n, s in seqs],
        "ks": ks,
        "agreement": {},
        "translation": {},
    }

    fulls = backbone.embed_full([s for _n, s in seqs])
    t0 = time.time()
    for k in ks:
        rows = []
        for (name, seq), full in zip(seqs, fulls):
            n = len(seq)
            lo, hi = k + 1, n - k - 1
            if hi - lo < 4:
                continue
            idxs = sorted(rng.sample(range(lo, hi), min(args.n_positions, hi - lo)))
            crops = backbone.embed_crops([seq], k, positions=idxs)[0]
            for i, vec in zip(idxs, crops):
                own = F.cosine_similarity(vec, full[i], dim=0).item()
                j = rng.choice([p for p in range(lo, hi) if p != i])
                other = F.cosine_similarity(vec, full[j], dim=0).item()
                rows.append({"own": own, "other": other, "len": n})
        if rows:
            report["agreement"][str(k)] = {
                "n": len(rows),
                "own_mean": round(sum(r["own"] for r in rows) / len(rows), 4),
                "other_mean": round(sum(r["other"] for r in rows) / len(rows), 4),
                "own_min": round(min(r["own"] for r in rows), 4),
            }

    # translation invariance: same window content at absolute position 0 vs 100
    rows = []
    for name, seq in seqs[:3]:
        n = len(seq)
        centre = n // 2
        for k in ks:
            if centre - k < 0 or centre + k + 2 > n:
                continue
            window = seq[centre - k : centre + k + 1]
            alone = backbone.embed_strings([window])[0][k]
            flank_l = "".join(rng.choice("ACGU") for _ in range(args.flank))
            flank_r = "".join(rng.choice("ACGU") for _ in range(args.flank))
            shifted = backbone.embed_strings([flank_l + window + flank_r])[0][
                args.flank + k
            ]
            rows.append(
                {
                    "k": k,
                    "cos": F.cosine_similarity(alone, shifted, dim=0).item(),
                }
            )
    for k in ks:
        vals = [r["cos"] for r in rows if r["k"] == k]
        if vals:
            report["translation"][str(k)] = {
                "n": len(vals),
                "cos_mean": round(sum(vals) / len(vals), 4),
                "cos_min": round(min(vals), 4),
            }
    report["elapsed_seconds"] = round(time.time() - t0, 1)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
