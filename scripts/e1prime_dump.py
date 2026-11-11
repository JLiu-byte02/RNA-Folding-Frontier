#!/usr/bin/env python3
"""Dump per-sequence pair probabilities so the D1/D2 diagnostics need no GPU.

The P2 table only kept F1 per sequence, which is enough to compare conditions
but not to ask *why* a condition changed: whether recall or precision moved,
whether the logits were rescaled, whether the errors are long-range, or whether
a false pair links the kept window to a complementary distal stretch.

This writes, per sequence, every candidate cell whose probability exceeds
``--floor`` together with its probability, the true pairs, and the kept window.
Everything downstream is then pure numpy.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402
import torch  # noqa: E402

from rff.models.pairhead import PairMapHead  # noqa: E402


def read_jsonl(path: Path) -> list[dict]:
    with open(path) as fh:
        return [json.loads(line) for line in fh if line.strip()]


def load_cache(path: Path) -> dict[str, torch.Tensor]:
    embs: dict[str, torch.Tensor] = {}
    for shard in sorted(path.glob("shard_*.pt")):
        blob = torch.load(shard, map_location="cpu", weights_only=False)
        for rid, emb in zip(blob["ids"], blob["embs"]):
            embs[rid] = emb
    return embs


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--checkpoint", required=True)
    ap.add_argument("--eval", action="append", required=True, help="name=data.jsonl:cache_dir")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--floor", type=float, default=0.02)
    ap.add_argument("--window", type=int, default=32)
    ap.add_argument("--pair-budget", type=int, default=400 * 400 * 4)
    ap.add_argument("--device", default="cuda")
    args = ap.parse_args()

    device = args.device if torch.cuda.is_available() else "cpu"
    ckpt = torch.load(ROOT / args.checkpoint, map_location="cpu", weights_only=False)
    state = ckpt["model"]
    d_in = int(state["proj.1.weight"].shape[1])
    model = PairMapHead(d_in=d_in).to(device)
    model.load_state_dict(state)
    model.eval()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    for spec in args.eval:
        name, rest = spec.split("=", 1)
        data, cache = rest.split(":", 1)
        embs = load_cache(ROOT / cache)
        recs = [r for r in read_jsonl(ROOT / data) if r["id"] in embs]
        recs.sort(key=lambda r: len(r["seq"]))
        path = out_dir / f"{name}.jsonl"
        n_kept = 0
        with open(path, "w") as fh, torch.no_grad():
            start = 0
            while start < len(recs):
                cur, cur_max = [], 0
                while start < len(recs):
                    L = len(recs[start]["seq"])
                    m = max(cur_max, L)
                    if cur and (len(cur) + 1) * m * m > args.pair_budget:
                        break
                    cur.append(recs[start])
                    cur_max = m
                    start += 1
                n = max(len(r["seq"]) for r in cur)
                h = torch.zeros(len(cur), n, d_in)
                for k, r in enumerate(cur):
                    e = embs[r["id"]]
                    h[k, : e.shape[0]] = e.float()
                prob = torch.sigmoid(model(h.to(device))).cpu().numpy()
                for k, r in enumerate(cur):
                    L = len(r["seq"])
                    p = prob[k][:L, :L]
                    ii, jj = np.where(np.triu(p, 1) >= args.floor)
                    lo = max(0, L // 2 - args.window)
                    hi = min(L, L // 2 + args.window + 1)
                    row = {
                        "id": r["id"],
                        "family": r["family"],
                        "source": r.get("source"),
                        "length": L,
                        "kept": [lo, hi],
                        "truth": [[int(i), int(j)] for i, j in r["pairs"]],
                        "scored": {
                            f"{int(i)},{int(j)}": round(float(p[i, j]), 5)
                            for i, j in zip(ii, jj)
                        },
                    }
                    fh.write(json.dumps(row) + "\n")
                    n_kept += 1
        print(json.dumps({"eval": name, "sequences": n_kept, "out": str(path)}), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
