#!/usr/bin/env python3
"""Is the clan probe weak because of the method, or because of the data?

The diagnostics report a clan probe at roughly chance level on id_test, which
would be a strong negative result for the whole mechanism story. Before that is
believed, the probe itself has to be shown to work somewhere:

*train-holdout*  a random 80/20 split of the training sequences, same probe,
  evaluated on held-out sequences of the same clans. If this is also chance, the
  probe is broken. If it is high, the representation carries clan identity and
  the id_test number is telling us something real about the split.

*per-family*  clan labels include families that are singleton clans (mapped to
  their accession), so the same measurement is repeated with family labels.
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
from sklearn.metrics import balanced_accuracy_score, top_k_accuracy_score  # noqa: E402


def read_jsonl(path: Path) -> list[dict]:
    with open(path) as fh:
        return [json.loads(line) for line in fh if line.strip()]


def load_cache(path: Path) -> dict[str, torch.Tensor]:
    embs = {}
    for shard in sorted(path.glob("shard_*.pt")):
        blob = torch.load(shard, map_location="cpu", weights_only=False)
        for rid, emb in zip(blob["ids"], blob["embs"]):
            embs[rid] = emb
    return embs


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--split", default="splits/s1prime_filtered/train.jsonl")
    ap.add_argument("--cache", default="data/cache/giga_a_train")
    ap.add_argument("--label", default="clan", choices=["clan", "family"])
    ap.add_argument("--min-per-class", type=int, default=5)
    ap.add_argument("--holdout", type=float, default=0.2)
    ap.add_argument("--C", type=float, default=1.0)
    ap.add_argument("--seed", type=int, default=17)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    recs = read_jsonl(ROOT / args.split)
    embs = load_cache(ROOT / args.cache)
    recs = [r for r in recs if r["id"] in embs]
    counts = Counter(r[args.label] for r in recs)
    keep = {k for k, v in counts.items() if v >= args.min_per_class}
    recs = [r for r in recs if r[args.label] in keep]

    X = np.stack([embs[r["id"]].float().mean(0).numpy() for r in recs])
    y = np.array([r[args.label] for r in recs])
    mean, std = X.mean(0), X.std(0) + 1e-6
    X = (X - mean) / std

    rng = np.random.default_rng(args.seed)
    idx = rng.permutation(len(X))
    cut = int((1 - args.holdout) * len(X))
    tr, te = idx[:cut], idx[cut:]
    clf = LogisticRegression(max_iter=2000, C=args.C)
    clf.fit(X[tr], y[tr])
    pred = clf.predict(X[te])

    # Nearest-centroid baseline: a probe-free measure of how separable the classes
    # are in mean-pooled embedding space.
    classes = np.unique(y[tr])
    cents = np.stack([X[tr][y[tr] == c].mean(0) for c in classes])
    d = ((X[te][:, None, :] - cents[None, :, :]) ** 2).sum(-1)
    nn = classes[d.argmin(1)]

    proba = clf.predict_proba(X[te])
    report = {
        "split": args.split,
        "cache": args.cache,
        "label": args.label,
        "n_classes": int(len(classes)),
        "n_train": int(len(tr)),
        "n_test": int(len(te)),
        "chance": round(1.0 / len(classes), 4),
        "probe_acc": round(float((pred == y[te]).mean()), 4),
        "probe_balanced_acc": round(float(balanced_accuracy_score(y[te], pred)), 4),
        "probe_top5": round(
            float(top_k_accuracy_score(y[te], proba, k=min(5, len(classes)), labels=classes)), 4
        ),
        "centroid_acc": round(float((nn == y[te]).mean()), 4),
    }
    out = ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
