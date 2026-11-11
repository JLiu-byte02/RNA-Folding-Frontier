#!/usr/bin/env python3
"""Train the pair-map head on cached embeddings and evaluate base-pair F1.

Only the head is trained; the backbone is a cache. Length-bucketed batches keep
padding waste small, and evaluation reports both sequence-micro F1 and
family-macro F1, because the two answer different questions ("how many pairs
are right" vs "is every family handled").
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402
import torch  # noqa: E402

from rff.analysis.labels import pair_prf  # noqa: E402
from rff.models.dp import pair_map_to_pairs  # noqa: E402
from rff.models.pairhead import PairMapHead, pair_bce, pairs_to_matrix  # noqa: E402


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


def batches(records: list[dict], size: int, shuffle: bool, seed: int):
    order = sorted(records, key=lambda r: len(r["seq"]))
    chunks = [order[i : i + size] for i in range(0, len(order), size)]
    if shuffle:
        rng = np.random.default_rng(seed)
        rng.shuffle(chunks)
    yield from chunks


def collate(batch: list[dict], embs: dict[str, torch.Tensor], device: str):
    n = max(len(rec["seq"]) for rec in batch)
    b = len(batch)
    h = torch.zeros(b, n, next(iter(embs.values())).shape[-1])
    lens = torch.tensor([len(rec["seq"]) for rec in batch])
    target = torch.zeros(b, n, n)
    for k, rec in enumerate(batch):
        e = embs[rec["id"]]
        h[k, : e.shape[0]] = e.float()
        target[k] = pairs_to_matrix(rec["pairs"], n)
    return h.to(device), target.to(device), lens.to(device)


@torch.no_grad()
def evaluate(
    model: PairMapHead,
    records: list[dict],
    embs: dict[str, torch.Tensor],
    device: str,
    threshold: float,
    batch_size: int,
) -> dict:
    model.eval()
    per_seq = []
    per_family: dict[str, list[dict]] = defaultdict(list)
    for batch in batches(records, batch_size, shuffle=False, seed=0):
        h, target, lens = collate(batch, embs, device)
        logits = model(h, mask=None)
        prob = torch.sigmoid(logits).cpu().numpy()
        for k, rec in enumerate(batch):
            n = len(rec["seq"])
            pred = pair_map_to_pairs(prob[k][:n, :n], threshold=threshold)
            truth = [tuple(p) for p in rec["pairs"]]
            m = pair_prf(truth, pred)
            per_seq.append(m)
            per_family[rec["family"]].append(m)
    macro = {
        fam: float(np.mean([m["f1"] for m in ms]))
        for fam, ms in per_family.items()
    }
    return {
        "n": len(per_seq),
        "micro_f1": round(float(np.mean([m["f1"] for m in per_seq])), 4),
        "precision": round(float(np.mean([m["precision"] for m in per_seq])), 4),
        "recall": round(float(np.mean([m["recall"] for m in per_seq])), 4),
        "family_macro_f1": round(float(np.mean(list(macro.values()))), 4)
        if macro
        else 0.0,
        "per_family": {k: round(v, 4) for k, v in sorted(macro.items())},
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="data/processed/archiveii")
    ap.add_argument("--cache", default="data/cache/archiveii_full")
    ap.add_argument("--out", required=True)
    ap.add_argument("--epochs", type=int, default=30)
    ap.add_argument("--batch-size", type=int, default=4)
    ap.add_argument("--lr", type=float, default=1e-3)
    ap.add_argument("--weight-decay", type=float, default=0.01)
    ap.add_argument("--pos-weight", type=float, default=20.0)
    ap.add_argument("--eval-every", type=int, default=5)
    ap.add_argument("--threshold", type=float, default=0.5)
    ap.add_argument("--seed", type=int, default=17)
    ap.add_argument("--device", default="cuda")
    args = ap.parse_args()

    torch.manual_seed(args.seed)
    np.random.seed(args.seed)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    train = read_jsonl(ROOT / args.data_dir / "train.jsonl")
    val = read_jsonl(ROOT / args.data_dir / "val.jsonl")
    embs = load_cache(ROOT / args.cache)
    train = [r for r in train if r["id"] in embs]
    val = [r for r in val if r["id"] in embs]
    if not train or not val:
        print("cache does not cover the dataset")
        return 1

    d_in = next(iter(embs.values())).shape[-1]
    model = PairMapHead(d_in=d_in).to(args.device)
    opt = torch.optim.AdamW(
        model.parameters(), lr=args.lr, weight_decay=args.weight_decay
    )
    steps = args.epochs * math.ceil(len(train) / args.batch_size)
    sched = torch.optim.lr_scheduler.OneCycleLR(
        opt, max_lr=args.lr, total_steps=steps, pct_start=0.1
    )

    log_path = out / "log.jsonl"
    best = -1.0
    t0 = time.time()
    with open(log_path, "w") as log:
        for epoch in range(1, args.epochs + 1):
            model.train()
            losses = []
            for batch in batches(train, args.batch_size, True, args.seed + epoch):
                h, target, lens = collate(batch, embs, args.device)
                logits = model(h)
                loss = pair_bce(logits, target, lens, args.pos_weight)
                opt.zero_grad(set_to_none=True)
                loss.backward()
                torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
                opt.step()
                sched.step()
                losses.append(float(loss))
            entry = {
                "epoch": epoch,
                "loss": round(float(np.mean(losses)), 4),
                "lr": round(float(sched.get_last_lr()[0]), 6),
                "elapsed": round(time.time() - t0, 1),
            }
            if epoch % args.eval_every == 0 or epoch == args.epochs:
                entry["val"] = evaluate(
                    model, val, embs, args.device, args.threshold, args.batch_size
                )
                if entry["val"]["family_macro_f1"] > best:
                    best = entry["val"]["family_macro_f1"]
                    torch.save(
                        {"model": model.state_dict(), "args": vars(args)},
                        out / "best.pt",
                    )
                    entry["saved"] = True
            log.write(json.dumps(entry) + "\n")
            log.flush()
            print(json.dumps(entry), flush=True)

    torch.save({"model": model.state_dict(), "args": vars(args)}, out / "last.pt")
    print(json.dumps({"best_family_macro_f1": best}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
