#!/usr/bin/env python3
"""Cache frozen-backbone embeddings so many small heads can be trained cheaply.

The backbone is frozen, so its output depends only on (backbone, k, sequence).
Caching turns a k-sweep from "re-embed 8k sequences per run" into "embed once,
train many heads". Storage is fp16: 1280 dims x 400 nt x 2 B ~ 1 MB per
sequence, ~8 GB for the pre-registered 8k cap.

Reads a JSONL dataset (``id``, ``seq``) or ArchiveII ``data.json`` and writes
shards of ``{ids, lengths, embs}`` plus a manifest with the measured cost.
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


def load_records(path: Path, limit: int | None) -> list[dict]:
    recs: list[dict] = []
    if path.suffix == ".json":
        data = json.loads(path.read_text())
        for rid, rec in data.items():
            recs.append({"id": rid, "seq": rec["sequence"]})
    else:
        with open(path) as fh:
            for line in fh:
                if line.strip():
                    recs.append(json.loads(line))
    return recs[:limit] if limit else recs


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True, help="jsonl dataset or ArchiveII data.json")
    ap.add_argument("--out", required=True, help="output directory")
    ap.add_argument("--backbone", default="rinalmo-giga")
    ap.add_argument("--model-path", default="data/models/rinalmo-giga")
    ap.add_argument("--k", default="full", help="'full' or an integer half-width")
    ap.add_argument("--max-len", type=int, default=1024)
    ap.add_argument("--min-len", type=int, default=0)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--shard-size", type=int, default=256)
    ap.add_argument("--dtype", default="float16", choices=["float16", "float32"])
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    recs = load_records(Path(args.data), args.limit)
    total_in = len(recs)
    recs = [r for r in recs if args.min_len <= len(r["seq"]) <= args.max_len]
    if not recs:
        print("no records in the length window")
        return 1

    backbone = load_backbone(args.backbone, path=args.model_path)
    dtype = torch.float16 if args.dtype == "float16" else torch.float32
    k = None if args.k == "full" else int(args.k)

    t0 = time.time()
    residues = 0
    dim = 0
    shard_ids: list[str] = []
    shard_embs: list[torch.Tensor] = []
    n_shards = 0

    def flush() -> None:
        nonlocal shard_ids, shard_embs, n_shards
        if not shard_ids:
            return
        torch.save(
            {
                "ids": shard_ids,
                "lengths": torch.tensor([e.shape[0] for e in shard_embs]),
                "embs": [e.to(dtype) for e in shard_embs],
            },
            out / f"shard_{n_shards:04d}.pt",
        )
        n_shards += 1
        shard_ids, shard_embs = [], []

    for n, rec in enumerate(recs, 1):
        seq = rec["seq"][: args.max_len]
        residues += len(seq)
        if k is None:
            emb = backbone.embed_full([seq])[0]
        else:
            emb = backbone.embed_crops([seq], k)[0]
        shard_ids.append(rec["id"])
        shard_embs.append(emb)
        dim = int(emb.shape[-1])
        if len(shard_ids) >= args.shard_size:
            flush()
        if n % 200 == 0:
            dt = time.time() - t0
            print(f"{n}/{len(recs)}  {dt:.0f}s  {dt / n:.3f}s/seq", flush=True)
    flush()

    elapsed = time.time() - t0
    manifest = {
        "backbone": args.backbone,
        "k": args.k,
        "data": str(args.data),
        "dtype": args.dtype,
        "max_len": args.max_len,
        "min_len": args.min_len,
        "records_in_file": total_in,
        "records_cached": len(recs),
        "residues": residues,
        "shards": n_shards,
        "shard_size": args.shard_size,
        "elapsed_seconds": round(elapsed, 1),
        "seconds_per_sequence": round(elapsed / len(recs), 4),
        "seconds_per_residue": round(elapsed / max(1, residues), 7),
        "dim": dim,
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2))
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
