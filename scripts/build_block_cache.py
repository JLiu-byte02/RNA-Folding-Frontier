#!/usr/bin/env python3
"""Cache P2-block embeddings (Amendment 2 A2.2).

Round 1 defined the context intervention per position: embed position ``i`` with
its own native ``+/-32`` nt and everything else replaced by the condition's
donor content. A single forward pass per sequence cannot express that, which is
why P2-central is only an exploratory readout. This script implements the
per-position version in blocks of ``--block`` nt, which costs
``ceil(L/block)`` passes instead of ``L``:

    for each block [s, e):
        composite = donor everywhere except [s - flank, e + flank), native there
        embedding = frozen_backbone(composite)[s:e]

Positions at the sequence ends have their flank clipped by the sequence
boundary, the same convention as ``FrozenBackbone.embed_crops``. Every position
is therefore embedded with at least ``min(flank, distance to end)`` nt of its
own native flank, and each position is produced by exactly one pass.

The donor strings come from ``build_block_donors.py`` and are full-length, so
the same donor file can be reused if the block/flank geometry changes.

Output format matches ``build_cache.py`` (shards of ``ids/lengths/embs``) so
``train_head_s1.py`` reads it without changes.
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


def block_spans(length: int, block: int) -> list[tuple[int, int]]:
    return [(s, min(s + block, length)) for s in range(0, length, block)]


def block_context(
    native: str, donor: str, start: int, end: int, flank: int
) -> str:
    """Composite sequence: native within ``[start-flank, end+flank)``."""
    if len(native) != len(donor):
        raise ValueError("native and donor must have the same length")
    lo = max(0, start - flank)
    hi = min(len(native), end + flank)
    return donor[:lo] + native[lo:hi] + donor[hi:]


def assemble_sequence(
    backbone,
    native: str,
    donor: str,
    block: int,
    flank: int,
    batch_size: int = 64,
) -> "torch.Tensor":
    """One [L, d] embedding built from ``ceil(L/block)`` passes."""
    spans = block_spans(len(native), block)
    composites = [block_context(native, donor, s, e, flank) for s, e in spans]
    embedded = backbone.embed_strings(composites, batch_size=batch_size)
    rows = []
    for (s, e), emb in zip(spans, embedded):
        if emb.shape[0] != len(native):
            raise AssertionError(
                f"pass for block [{s},{e}) returned {emb.shape[0]} positions, "
                f"expected {len(native)}"
            )
        rows.append(emb[s:e])
    assembled = torch.cat(rows, dim=0)
    if assembled.shape[0] != len(native):
        raise AssertionError("block assembly did not cover the sequence once")
    return assembled


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True, help="jsonl of target records")
    ap.add_argument("--donors", required=True, help="json from build_block_donors.py")
    ap.add_argument("--out", required=True)
    ap.add_argument("--condition", required=True, help="label written to the manifest")
    ap.add_argument("--backbone", default="rinalmo-giga")
    ap.add_argument("--model-path", default="data/models/rinalmo-giga")
    ap.add_argument("--block", type=int, default=16)
    ap.add_argument("--flank", type=int, default=32)
    ap.add_argument("--max-len", type=int, default=1024)
    ap.add_argument("--min-len", type=int, default=0)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--shard-size", type=int, default=128)
    ap.add_argument("--batch-size", type=int, default=64)
    ap.add_argument("--dtype", default="float16", choices=["float16", "float32"])
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    donors: dict[str, str] = json.loads(Path(args.donors).read_text())
    recs = [json.loads(line) for line in open(Path(args.data)) if line.strip()]
    total_in = len(recs)
    recs = [
        r
        for r in recs
        if args.min_len <= len(r["seq"]) <= args.max_len and r["id"] in donors
    ]
    if args.limit:
        recs = recs[: args.limit]
    if not recs:
        print("no records with donors inside the length window")
        return 1

    backbone = load_backbone(args.backbone, path=args.model_path)
    dtype = torch.float16 if args.dtype == "float16" else torch.float32

    t0 = time.time()
    residues = 0
    passes = 0
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
        native = rec["seq"][: args.max_len]
        donor = donors[rec["id"]][: args.max_len]
        assembled = assemble_sequence(
            backbone, native, donor, args.block, args.flank, args.batch_size
        )
        dim = int(assembled.shape[-1])
        residues += len(native)
        passes += len(block_spans(len(native), args.block))
        shard_ids.append(rec["id"])
        shard_embs.append(assembled)
        if len(shard_ids) >= args.shard_size:
            flush()
        if n % 50 == 0:
            dt = time.time() - t0
            print(f"{n}/{len(recs)}  {dt:.0f}s  {dt / n:.2f}s/seq", flush=True)
    flush()

    elapsed = time.time() - t0
    manifest = {
        "backbone": args.backbone,
        "condition": args.condition,
        "geometry": "p2-block",
        "block": args.block,
        "flank": args.flank,
        "data": str(args.data),
        "donors": str(args.donors),
        "dtype": args.dtype,
        "max_len": args.max_len,
        "min_len": args.min_len,
        "records_in_file": total_in,
        "records_cached": len(recs),
        "residues": residues,
        "forward_passes": passes,
        "passes_per_sequence": round(passes / len(recs), 3),
        "shards": n_shards,
        "shard_size": args.shard_size,
        "elapsed_seconds": round(elapsed, 1),
        "seconds_per_sequence": round(elapsed / len(recs), 4),
        "dim": dim,
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2))
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
