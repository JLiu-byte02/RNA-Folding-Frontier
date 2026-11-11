#!/usr/bin/env python3
"""Measure frozen-backbone throughput so the E1' grid can be costed honestly.

Reports seconds per residue for (a) one full-context pass and (b) crop passes at
several k, then extrapolates to the planned grid.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import torch  # noqa: E402

from rff.models.backbone import load_backbone  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--backbone", default="rinalmo-giga")
    ap.add_argument("--path", default="data/models/rinalmo-giga")
    ap.add_argument("--n-seqs", type=int, default=4)
    ap.add_argument("--length", type=int, default=300)
    ap.add_argument("--ks", default="8,32,128")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    seq = ("ACGU" * (args.length // 4 + 1))[: args.length]
    seqs = [seq] * args.n_seqs

    t0 = time.time()
    backbone = load_backbone(args.backbone, path=args.path)
    load_s = time.time() - t0

    report: dict = {
        "backbone": args.backbone,
        "device": torch.cuda.get_device_name(0) if torch.cuda.is_available() else "cpu",
        "load_seconds": round(load_s, 1),
        "n_seqs": args.n_seqs,
        "length": args.length,
        "full": {},
        "crops": {},
    }

    torch.cuda.reset_peak_memory_stats()
    t0 = time.time()
    emb = backbone.embed_full(seqs)
    dt = time.time() - t0
    report["full"] = {
        "seconds": round(dt, 3),
        "residues": args.n_seqs * args.length,
        "seconds_per_residue": round(dt / (args.n_seqs * args.length), 7),
        "dim": int(emb[0].shape[-1]),
        "peak_gb": round(torch.cuda.max_memory_allocated() / 1e9, 2),
    }

    for k in [int(x) for x in args.ks.split(",")]:
        t0 = time.time()
        backbone.embed_crops(seqs, k)
        dt = time.time() - t0
        residues = args.n_seqs * args.length
        report["crops"][str(k)] = {
            "seconds": round(dt, 3),
            "seconds_per_residue": round(dt / residues, 6),
            "window": 2 * k + 1,
        }

    Path(args.out).write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
