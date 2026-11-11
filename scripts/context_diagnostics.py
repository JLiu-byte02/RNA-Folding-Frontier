#!/usr/bin/env python3
"""Measure what a distal swap did to the representation (Amendment 1.3).

Two numbers per condition, as pre-registered:

*distortion*  ``1 - cos(h_swap[i], h_native[i])``.  Reported separately for the
  positions that were rewritten (*distal*) and for the positions that were kept
  (*local*).  Local distortion is not zero even though the local tokens are
  untouched: attention is global, so the encoder's activations at kept positions
  move when their context changes.  That is the honest measure of how much of
  the intervention is indirect.

*clan decodability*  balanced accuracy of a linear probe trained on native
  training embeddings (mean-pooled) and applied to the swapped embeddings of the
  evaluated split.  The probe never sees the swapped split's labels during
  fitting, so a drop in accuracy is a loss of clan identity carried by the
  representation.

The point of computing both is that a swap can lower F1 either because it removed
family identity (the hypothesis) or because it moved the input off the
pretraining distribution (a nuisance).  Distortion measures the nuisance; clan
decodability measures the hypothesis.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402
import torch  # noqa: E402
from sklearn.linear_model import LogisticRegression  # noqa: E402
from sklearn.metrics import balanced_accuracy_score  # noqa: E402


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


def distal_mask(length: int, window: int) -> np.ndarray:
    mask = np.ones(length, dtype=bool)
    lo, hi = max(0, length // 2 - window), min(length, length // 2 + window + 1)
    mask[lo:hi] = False
    return mask


def cos_rows(a: torch.Tensor, b: torch.Tensor) -> torch.Tensor:
    a = a.float()
    b = b.float()
    na = a.norm(dim=-1).clamp(min=1e-6)
    nb = b.norm(dim=-1).clamp(min=1e-6)
    return (a * b).sum(-1) / (na * nb)


def pool(embs: dict[str, torch.Tensor], ids: list[str]) -> np.ndarray:
    return np.stack([embs[i].float().mean(0).numpy() for i in ids])


def fit_probe(embs, records, min_per_clan: int, C: float):
    ids = [r["id"] for r in records if r["id"] in embs]
    lookup = {r["id"]: r["clan"] for r in records}
    labels = [lookup[i] for i in ids]
    counts = Counter(labels)
    keep = [k for k, v in counts.items() if v >= min_per_clan]
    idx = [n for n, y in enumerate(labels) if y in keep]
    X = pool(embs, [ids[n] for n in idx])
    y = [labels[n] for n in idx]
    X = (X - X.mean(0)) / (X.std(0) + 1e-6)
    clf = LogisticRegression(max_iter=2000, C=C, n_jobs=-1)
    clf.fit(X, y)
    return clf, (X.mean(0), X.std(0) + 1e-6), len(keep)


def probe_acc(clf, scaler, embs, records) -> dict:
    mean, std = scaler
    lookup = {r["id"]: r["clan"] for r in records}
    known = set(clf.classes_)
    ids = [r["id"] for r in records if r["id"] in embs and lookup[r["id"]] in known]
    if not ids:
        return {"n": 0, "balanced_acc": None, "acc": None}
    X = (pool(embs, ids) - mean) / std
    y = [lookup[i] for i in ids]
    pred = clf.predict(X)
    return {
        "n": len(ids),
        "clans": len(set(y)),
        "balanced_acc": round(float(balanced_accuracy_score(y, pred)), 4),
        "acc": round(float((pred == np.array(y)).mean()), 4),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--train", default="splits/s1prime_filtered/train.jsonl")
    ap.add_argument("--train-cache", default="data/cache/giga_a_train")
    ap.add_argument("--split", required=True, help="jsonl of the evaluated split")
    ap.add_argument("--native-cache", required=True)
    ap.add_argument("--condition", action="append", required=True, help="name:cache_dir")
    ap.add_argument("--window", type=int, default=32)
    ap.add_argument("--min-per-clan", type=int, default=5)
    ap.add_argument("--C", type=float, default=1.0)
    ap.add_argument("--out", required=True)
    ap.add_argument(
        "--per-seq-out",
        default=None,
        help="jsonl with one row per (sequence, condition) for the F1 regression",
    )
    args = ap.parse_args()

    train = read_jsonl(ROOT / args.train)
    train_embs = load_cache(ROOT / args.train_cache)
    clf, scaler, n_clans = fit_probe(train_embs, train, args.min_per_clan, args.C)
    print(json.dumps({"probe_clans": n_clans, "train_n": len(train_embs)}), flush=True)
    del train_embs

    recs = read_jsonl(ROOT / args.split)
    native = load_cache(ROOT / args.native_cache)
    recs = [r for r in recs if r["id"] in native]

    rows = {
        "split": str(args.split),
        "window": args.window,
        "probe_clans": n_clans,
        "n": len(recs),
        "conditions": {},
    }
    native_probe = probe_acc(clf, scaler, native, recs)
    rows["conditions"]["a"] = {"native": True, "probe": native_probe}
    print(json.dumps({"cond": "a", **native_probe}), flush=True)

    keep = set(clf.classes_)
    lookup = {r["id"]: r["clan"] for r in recs}
    seq_rows: list[dict] = []
    classes = list(clf.classes_)
    class_index = {c: n for n, c in enumerate(classes)}

    def per_seq_probe(embs, rec: dict) -> dict:
        """p(true clan) and p(best clan) for one sequence's pooled embedding."""
        clan = lookup[rec["id"]]
        if clan not in keep:
            return {"p_true": None, "p_max": None}
        mean, std = scaler
        x = (embs[rec["id"]].float().mean(0).numpy() - mean) / std
        p = clf.predict_proba(x[None, :])[0]
        return {"p_true": float(p[class_index[clan]]), "p_max": float(p.max())}

    for rec in recs:
        n = native[rec["id"]].shape[0]
        seq_rows.append(
            {
                "id": rec["id"],
                "family": rec["family"],
                "clan": rec["clan"],
                "condition": "a",
                "length": n,
                "distortion_all": 0.0,
                "distortion_distal": 0.0 if n > 2 * args.window + 1 else None,
                "distortion_local": 0.0,
                **per_seq_probe(native, rec),
            }
        )

    for spec in args.condition:
        name, cache_dir = spec.split(":", 1)
        cond = load_cache(ROOT / cache_dir)
        d_dist, d_loc, d_all = [], [], []
        for rec in recs:
            h0, h1 = native[rec["id"]], cond.get(rec["id"])
            if h1 is None:
                continue
            n = min(h0.shape[0], h1.shape[0])
            c = 1.0 - cos_rows(h0[:n], h1[:n]).numpy()
            mask = distal_mask(n, args.window)
            row = {
                "id": rec["id"],
                "family": rec["family"],
                "clan": rec["clan"],
                "condition": name,
                "length": n,
                "distortion_all": round(float(c.mean()), 5),
            }
            d_all.append(float(c.mean()))
            if mask.any():
                d_dist.append(float(c[mask].mean()))
                row["distortion_distal"] = round(float(c[mask].mean()), 5)
            if (~mask).any():
                d_loc.append(float(c[~mask].mean()))
                row["distortion_local"] = round(float(c[~mask].mean()), 5)
            probe = per_seq_probe(cond, rec)
            row.update(probe)
            seq_rows.append(row)
        pr = probe_acc(clf, scaler, cond, recs)
        rows["conditions"][name] = {
            "cache": cache_dir,
            "n": len(d_all),
            "distortion_all": round(float(np.mean(d_all)), 5),
            "distortion_distal": round(float(np.mean(d_dist)), 5) if d_dist else None,
            "distortion_local": round(float(np.mean(d_loc)), 5) if d_loc else None,
            "distortion_local_p90": round(float(np.percentile(d_loc, 90)), 5) if d_loc else None,
            "probe": pr,
        }
        print(json.dumps({"cond": name, **rows["conditions"][name]}), flush=True)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(rows, indent=2) + "\n")
    print(f"wrote {out}")
    if args.per_seq_out:
        p = Path(args.per_seq_out)
        p.parent.mkdir(parents=True, exist_ok=True)
        with open(p, "w") as fh:
            for row in seq_rows:
                fh.write(json.dumps(row) + "\n")
        print(f"wrote {len(seq_rows)} per-sequence rows -> {p}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
