#!/usr/bin/env python3
"""Train the decoder-B pair head on S1' and evaluate it under swap conditions.

This is the P2 harness of PREREGISTRATION.md Amendment 1.2: **one** set of
weights is trained on native inputs (condition a) and then evaluated on the same
sequences with their distal context replaced (conditions b/d/e). The only thing
that changes between rows of the P2 table is the input, so a difference in F1 is
a difference caused by distal context.

Model selection never sees ``id_test`` or ``test_withheld``: a family-disjoint
slice of the training split is held out for checkpoint selection
(``--val-frac``, stratified by clan, seed fixed). ``id_test``/``test_withheld``
are only ever used by ``--eval``.

Batch size is chosen per bucket so that ``batch * max_len^2`` stays under a fixed
budget: the head materialises a dense [B, L, L, 3d] tensor, so a fixed batch size
either wastes memory on short sequences or OOMs on the long ones (p99 is ~1.5 kb).
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


def length_band(length: int) -> str:
    for upper, name in ((80, "L<80"), (150, "80-149"), (300, "150-299"), (10**9, "L>=300")):
        if length < upper:
            return name
    return "L>=300"


def bucket_batches(records: list[dict], budget: int, shuffle: bool, seed: int):
    """Length-sorted buckets with batch size set by the L^2 memory budget."""
    order = sorted(records, key=lambda r: len(r["seq"]))
    chunks, cur, cur_max = [], [], 0
    for rec in order:
        m = max(cur_max, len(rec["seq"]))
        if cur and (len(cur) + 1) * m * m > budget:
            chunks.append(cur)
            cur, cur_max = [rec], len(rec["seq"])
        else:
            cur.append(rec)
            cur_max = m
    if cur:
        chunks.append(cur)
    if shuffle:
        np.random.default_rng(seed).shuffle(chunks)
    yield from chunks


def collate(batch: list[dict], embs: dict[str, torch.Tensor], device: str, d_in: int):
    n = max(len(rec["seq"]) for rec in batch)
    b = len(batch)
    h = torch.zeros(b, n, d_in)
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
    budget: int,
    per_seq_out: list[dict] | None = None,
) -> dict:
    model.eval()
    d_in = int(model.proj[1].in_features)
    per_seq, per_family = [], defaultdict(list)
    bands: dict[str, list[float]] = defaultdict(list)
    sources: dict[str, list[float]] = defaultdict(list)
    for batch in bucket_batches(records, budget, False, 0):
        h, _target, _lens = collate(batch, embs, device, d_in)
        prob = torch.sigmoid(model(h)).cpu().numpy()
        for k, rec in enumerate(batch):
            n = len(rec["seq"])
            pred = pair_map_to_pairs(prob[k][:n, :n], threshold=threshold)
            m = pair_prf([tuple(p) for p in rec["pairs"]], pred)
            per_seq.append(m)
            per_family[rec["family"]].append(m["f1"])
            bands[length_band(n)].append(m["f1"])
            sources[rec.get("source", "?")].append(m["f1"])
            if per_seq_out is not None:
                per_seq_out.append({"id": rec["id"], "f1": round(m["f1"], 6)})
    fam_macro = {f: float(np.mean(v)) for f, v in per_family.items()}
    return {
        "n": len(per_seq),
        "micro_f1": round(float(np.mean([m["f1"] for m in per_seq])), 4),
        "precision": round(float(np.mean([m["precision"] for m in per_seq])), 4),
        "recall": round(float(np.mean([m["recall"] for m in per_seq])), 4),
        "family_macro_f1": round(float(np.mean(list(fam_macro.values()))), 4)
        if fam_macro
        else 0.0,
        "by_band": {b: round(float(np.mean(v)), 4) for b, v in sorted(bands.items())},
        "by_source": {s: round(float(np.mean(v)), 4) for s, v in sorted(sources.items())},
        "n_source": {s: len(v) for s, v in sorted(sources.items())},
        "n_band": {b: len(v) for b, v in sorted(bands.items())},
        "per_family": {k: round(v, 4) for k, v in sorted(fam_macro.items())},
    }


def split_families(records: list[dict], frac: float, seed: int) -> tuple[list[dict], list[dict]]:
    """Hold out whole families (not sequences) so val cannot leak sequence variants."""
    fams = sorted({r["family"] for r in records})
    rng = np.random.default_rng(seed)
    rng.shuffle(fams)
    n_val = max(1, int(round(frac * len(fams))))
    val_fams = set(fams[:n_val])
    fit = [r for r in records if r["family"] not in val_fams]
    val = [r for r in records if r["family"] in val_fams]
    return fit, val


def parse_eval_specs(specs: list[str]) -> list[tuple[str, Path, Path]]:
    out = []
    for spec in specs:
        name, rest = spec.split("=", 1)
        data, cache = rest.split(":", 1)
        out.append((name, Path(data), Path(cache)))
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--train", default="splits/s1prime_filtered/train.jsonl")
    ap.add_argument(
        "--val",
        default=None,
        help="explicit validation jsonl (already cluster-disjoint); "
        "default holds out whole families from --train",
    )
    ap.add_argument("--train-cache", default="data/cache/giga_a_train")
    ap.add_argument("--eval", nargs="*", default=[], help="name=data.jsonl:cache_dir")
    ap.add_argument("--out", required=True)
    ap.add_argument("--epochs", type=int, default=12)
    ap.add_argument("--pair-budget", type=int, default=400 * 400 * 4)
    ap.add_argument("--lr", type=float, default=2e-3)
    ap.add_argument("--weight-decay", type=float, default=0.01)
    ap.add_argument("--pos-weight", type=float, default=20.0)
    ap.add_argument("--threshold", type=float, default=0.5)
    ap.add_argument("--val-frac", type=float, default=0.05)
    ap.add_argument("--seed", type=int, default=17)
    ap.add_argument("--device", default="cuda")
    ap.add_argument("--init", default=None, help="checkpoint to start from")
    ap.add_argument("--eval-only", action="store_true")
    args = ap.parse_args()

    torch.manual_seed(args.seed)
    np.random.seed(args.seed)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    device = args.device if torch.cuda.is_available() else "cpu"

    train_all = read_jsonl(ROOT / args.train)
    embs = load_cache(ROOT / args.train_cache)
    train_all = [r for r in train_all if r["id"] in embs]
    if not train_all:
        print("training cache does not cover the split")
        return 1
    d_in = int(next(iter(embs.values())).shape[-1])

    if args.val:
        fit = train_all
        val = [r for r in read_jsonl(ROOT / args.val) if r["id"] in embs]
        if not val:
            print("validation file does not overlap the training cache")
            return 1
        shared = {r["family"] for r in fit} & {r["family"] for r in val}
        print(
            json.dumps({"val_source": args.val, "families_shared_with_train": len(shared)}),
            flush=True,
        )
    else:
        fit, val = split_families(train_all, args.val_frac, args.seed)
    print(json.dumps({"fit": len(fit), "val": len(val), "d_in": d_in}), flush=True)

    eval_sets = []
    for name, data, cache in parse_eval_specs(args.eval):
        ce = load_cache(ROOT / cache)
        recs = [r for r in read_jsonl(ROOT / data) if r["id"] in ce]
        missing = len(read_jsonl(ROOT / data)) - len(recs)
        print(json.dumps({"eval": name, "n": len(recs), "missing": missing}), flush=True)
        eval_sets.append((name, recs, ce))

    model = PairMapHead(d_in=d_in).to(device)
    results: dict = {"args": vars(args), "fit": len(fit), "val": len(val)}

    if args.eval_only:
        ckpt = torch.load(ROOT / (args.init or Path(args.out) / "best.pt"), map_location="cpu")
        model.load_state_dict(ckpt["model"])
        model.to(device)
    else:
        opt = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
        steps = args.epochs * len(list(bucket_batches(fit, args.pair_budget, False, 0)))
        sched = torch.optim.lr_scheduler.OneCycleLR(
            opt, max_lr=args.lr, total_steps=max(1, steps), pct_start=0.15
        )
        best, best_val = -1.0, None
        t0 = time.time()
        log_path = out / "log.jsonl"
        with open(log_path, "w") as log:
            for epoch in range(1, args.epochs + 1):
                model.train()
                losses = []
                for batch in bucket_batches(fit, args.pair_budget, True, args.seed + epoch):
                    h, target, lens = collate(batch, embs, device, d_in)
                    logits = model(h)
                    loss = pair_bce(logits, target, lens, args.pos_weight)
                    opt.zero_grad(set_to_none=True)
                    loss.backward()
                    torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
                    opt.step()
                    sched.step()
                    losses.append(float(loss))
                vm = evaluate(model, val, embs, device, args.threshold, args.pair_budget)
                entry = {
                    "epoch": epoch,
                    "loss": round(float(np.mean(losses)), 4),
                    "val_micro_f1": vm["micro_f1"],
                    "val_family_macro_f1": vm["family_macro_f1"],
                    "elapsed": round(time.time() - t0, 1),
                }
                if vm["family_macro_f1"] > best:
                    best, best_val = vm["family_macro_f1"], entry
                    torch.save({"model": model.state_dict(), "args": vars(args)}, out / "best.pt")
                    entry["saved"] = True
                log.write(json.dumps(entry) + "\n")
                log.flush()
                print(json.dumps(entry), flush=True)
        results["best_val"] = best_val
        ckpt = torch.load(out / "best.pt", map_location="cpu")
        model.load_state_dict(ckpt["model"])
        model.to(device)

    evals = {}
    seq_dir = out / "per_seq"
    seq_dir.mkdir(parents=True, exist_ok=True)
    for name, recs, ce in eval_sets:
        rows: list[dict] = []
        evals[name] = evaluate(
            model, recs, ce, device, args.threshold, args.pair_budget, per_seq_out=rows
        )
        with open(seq_dir / f"{name}.jsonl", "w") as fh:
            for row in rows:
                fh.write(json.dumps(row) + "\n")
        print(
            json.dumps(
                {
                    "eval": name,
                    "micro_f1": evals[name]["micro_f1"],
                    "family_macro_f1": evals[name]["family_macro_f1"],
                    "n": evals[name]["n"],
                }
            ),
            flush=True,
        )
    results["evals"] = evals
    (out / "results.json").write_text(json.dumps(results, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
