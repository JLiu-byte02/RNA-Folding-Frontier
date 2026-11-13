#!/usr/bin/env python3
"""Fit the E-C2 coverage-gated decoder (Amendment 13 section A13.2).

The gate turns masked pseudo-perplexity into a per-sequence soft-constraint
scale,

    z(x) = (pppl(x) - mean) / sd
    G1:  s(x) = s_max * sigmoid(alpha - beta * z)
    G2:  s(x) = s_max * sigmoid(alpha - beta * z - gamma * log L)

and the decode is the A-energy-fixed arm with ``b = 0`` and ``scale = s(x)``,
i.e. the gate rescales the whole score matrix. This is the same object the
fixed-s comparator uses, so the comparison isolates the gate.

Folds are clan-grouped (Amendment 7 section A7.4): families in the same Rfam
clan never straddle folds; unclanned families act as singleton clans. Clans are
assigned largest-first to the currently smallest fold, which balances sizes
without looking at any outcome.

The head that produced the score matrices was trained on all of the train
split, so the CV is a *selector*-level estimate: it ranks gate settings without
holding out the head. Absolute F1 here is therefore optimistic; the withheld
number is the honest one. That limitation is recorded in the output JSON.

Counts, not predictions, cross the process boundary: 81 configurations times
17k sequences of pair sets would not fit in memory. Each worker returns
``tp/fp/fn`` per configuration and the parent aggregates family-macro F1.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import multiprocessing
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rff.analysis.labels import max_nested_subset  # noqa: E402
from rff.analysis.scoring import normalise_scored  # noqa: E402
from rff.decoders.nn_dp import NEG, logit_score, mask_allowed  # noqa: E402


def read_jsonl(path: Path) -> list[dict]:
    with open(path) as fh:
        return [json.loads(line) for line in fh if line.strip()]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_clans(path: Path) -> dict[str, str]:
    clans: dict[str, str] = {}
    with open(path) as fh:
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) < 2 or not parts[0].startswith("CL"):
                continue
            for family in parts[1:]:
                if family:
                    clans[family] = parts[0]
    return clans


def clan_folds(clans: dict[str, str], families: list[str], k: int) -> dict[str, int]:
    """Assign clans to folds, largest clan first into the smallest fold."""
    counts: dict[str, int] = {}
    for family in families:
        key = clans.get(family, family)
        counts[key] = counts.get(key, 0) + 1
    order = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))
    load = [0] * k
    assign: dict[str, int] = {}
    for key, n in order:
        fold = min(range(k), key=lambda f: (load[f], f))
        assign[key] = fold
        load[fold] += n
    return assign


def dense_scores(scored: dict, seq: str) -> np.ndarray:
    n = len(seq)
    prob = np.zeros((n, n))
    present = np.zeros((n, n), dtype=bool)
    for (i, j), p in scored.items():
        if i < j:
            prob[i, j] = prob[j, i] = p
            present[i, j] = present[j, i] = True
    score = logit_score(prob)
    allowed = mask_allowed(seq)
    tri = np.triu(np.ones((n, n), dtype=bool), 1)
    return np.where(present & allowed & tri, score, NEG)


def gate_scale(cfg: dict, z: float, log_len: float) -> float:
    arg = cfg["alpha"] - cfg["beta"] * z
    if cfg.get("gamma"):
        arg -= cfg["gamma"] * log_len
    return cfg["s_max"] / (1.0 + float(np.exp(-arg)))


_ITEMS: list[dict] = []
_CFGS: list[dict] = []


def _init(items: list[dict], cfgs: list[dict]) -> None:
    global _ITEMS, _CFGS
    _ITEMS = items
    _CFGS = cfgs


def _decode_item(idx: int) -> tuple[int, np.ndarray]:
    from rff.decoders.vienna_soft import ViennaSoftDecoder

    item = _ITEMS[idx]
    seq = item["seq"]
    scores = dense_scores(item["scored"], seq)
    truth = item["truth"]
    n_cfg = len(_CFGS)
    counts = np.zeros((n_cfg, 3), dtype=np.int32)
    for c, cfg in enumerate(_CFGS):
        scale = gate_scale(cfg, item["z"], item["log_len"])
        pred = set(ViennaSoftDecoder(scale=scale, offset=0.0).decode(seq, scores)[0])
        counts[c, 0] = len(pred & truth)
        counts[c, 1] = len(pred - truth)
        counts[c, 2] = len(truth - pred)
    return idx, counts


def f1_from_counts(tp: int, fp: int, fn: int) -> float:
    if tp + fp + fn == 0:
        return 1.0
    if tp == 0:
        return 0.0
    return 2 * tp / (2 * tp + fp + fn)


def family_macro_f1(counts: np.ndarray, families: list[str], idxs: list[int], c: int) -> float:
    per: dict[str, list[float]] = {}
    for i in idxs:
        sub = counts[i, c]
        per.setdefault(families[i], []).append(
            f1_from_counts(int(sub[0]), int(sub[1]), int(sub[2]))
        )
    means = [float(np.mean(v)) for v in per.values()]
    return float(np.mean(means)) if means else float("nan")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--train-dump", required=True)
    ap.add_argument("--train-split", required=True)
    ap.add_argument("--pppl", required=True)
    ap.add_argument("--clans", default="data/raw/rfam/clan_membership.txt")
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--max-items", type=int, default=0, help="0 = all; used for smoke tests")
    ap.add_argument("--seed", type=int, default=17)
    ap.add_argument("--min-true-pairs", type=int, default=4)
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    t_start = time.time()
    pppl_raw = json.loads(Path(args.pppl).read_text())
    pppl = pppl_raw["values"]
    all_pppl = np.array([v["pppl"] for v in pppl.values()], dtype=float)
    pppl_mean = float(all_pppl.mean())
    pppl_sd = float(all_pppl.std(ddof=0))

    seqs = {r["id"]: r["seq"] for r in read_jsonl(ROOT / args.train_split)}
    rows = read_jsonl(Path(args.train_dump))
    items: list[dict] = []
    n_no_seq = n_no_pppl = n_parse = n_below = 0
    for row in rows:
        seq = seqs.get(row["id"])
        if seq is None:
            n_no_seq += 1
            continue
        rec = pppl.get(row["id"])
        if rec is None:
            n_no_pppl += 1
            continue
        raw = row.get("truth_v3")
        if raw is None:
            n_parse += 1
            continue
        if len(raw) < args.min_true_pairs:
            n_below += 1
            continue
        truth = max_nested_subset({tuple(p) for p in raw})
        z = (rec["pppl"] - pppl_mean) / pppl_sd if pppl_sd > 0 else 0.0
        items.append(
            {
                "id": row["id"],
                "family": row["family"],
                "seq": seq,
                "scored": normalise_scored(row["scored"]),
                "truth": truth,
                "z": z,
                "log_len": float(np.log(len(seq))),
            }
        )
    if not items:
        raise SystemExit("no usable records")
    if args.max_items and len(items) > args.max_items:
        rng = np.random.default_rng(args.seed)
        keep = sorted(rng.choice(len(items), size=args.max_items, replace=False))
        items = [items[i] for i in keep]

    families = [it["family"] for it in items]
    clans = load_clans(ROOT / args.clans)
    assign = clan_folds(clans, families, args.folds)
    fold_of_item = [assign[clans.get(f, f)] for f in families]
    n_unclanned = sum(1 for f in set(families) if f not in clans)

    grid = []
    for s_max in (0.1875, 0.5, 1.0):
        for alpha in (-1.0, 0.0, 1.0):
            for beta in (0.0, 1.0, 2.0):
                grid.append({"gate": "G1", "s_max": s_max, "alpha": alpha, "beta": beta, "gamma": 0.0})
    for s_max in (0.1875, 0.5, 1.0):
        for alpha in (-1.0, 0.0, 1.0):
            for beta in (0.0, 1.0, 2.0):
                for gamma in (0.25, 0.5):
                    grid.append({"gate": "G2", "s_max": s_max, "alpha": alpha, "beta": beta, "gamma": gamma})

    counts = np.zeros((len(items), len(grid), 3), dtype=np.int32)
    global _ITEMS, _CFGS
    _ITEMS, _CFGS = items, grid
    try:
        ctx = multiprocessing.get_context("fork")
        executor = ProcessPoolExecutor(max_workers=args.workers, mp_context=ctx)
    except ValueError:
        executor = ProcessPoolExecutor(
            max_workers=args.workers, initializer=_init, initargs=(items, grid)
        )
    done = 0
    with executor as ex:
        for idx, c in ex.map(_decode_item, range(len(items)), chunksize=4):
            counts[idx] = c
            done += 1
            if done % 2000 == 0:
                print(f"[fit] {done}/{len(items)} items {time.time()-t_start:.0f}s", flush=True)

    folds = [[i for i in range(len(items)) if fold_of_item[i] == k] for k in range(args.folds)]
    results = []
    for c, cfg in enumerate(grid):
        fold_f1 = [family_macro_f1(counts, families, idxs, c) for idxs in folds if idxs]
        results.append(
            {
                **cfg,
                "folds": [round(f, 4) for f in fold_f1],
                "mean_f1": round(float(np.mean(fold_f1)), 4),
                "sd_f1": round(float(np.std(fold_f1)), 4),
                "pooled_f1": round(
                    family_macro_f1(counts, families, list(range(len(items))), c), 4
                ),
            }
        )

    g1 = [r for r in results if r["gate"] == "G1"]
    g2 = [r for r in results if r["gate"] == "G2"]
    best_g1 = max(g1, key=lambda r: r["mean_f1"])
    best_g2 = max(g2, key=lambda r: r["mean_f1"])

    out = {
        "status": "fitted; freeze before scoring withheld (A13.2)",
        "amendment": "AMENDMENT_13.md",
        "label_policy": {
            "version": "v3",
            "truth_mode": "nested",
            "min_true_pairs": args.min_true_pairs,
        },
        "head_source": str(args.train_dump),
        "head_source_sha256": sha256_file(Path(args.train_dump)),
        "pppl_source": str(args.pppl),
        "pppl_source_sha256": sha256_file(Path(args.pppl)),
        "n_sequences": len(items),
        "n_dropped": {
            "no_sequence": n_no_seq,
            "no_pppl": n_no_pppl,
            "parse_failed": n_parse,
            "below_min_pairs": n_below,
        },
        "train_pppl_mean": round(pppl_mean, 6),
        "train_pppl_sd": round(pppl_sd, 6),
        "folds": {
            "grouping": "clan (Rfam.clan_membership.txt; unclanned family = singleton clan)",
            "n_folds": args.folds,
            "n_clans": len(set(assign.values())) and len(assign),
            "n_unclanned_families": n_unclanned,
            "sizes": [int(len(f)) for f in folds],
        },
        "selection_objective": "mean over folds of family-macro F1",
        "caveat": (
            "score matrices come from a head trained on all folds, so this CV "
            "selects the gate without holding out the head; absolute F1 here is "
            "optimistic and the withheld number is the confirmatory one"
        ),
        "scale_semantics": "gate rescales the whole score matrix; offset b = 0",
        "grid": results,
        "best_g1": best_g1,
        "best_g2": best_g2,
        "elapsed_seconds": round(time.time() - t_start, 1),
        "workers": args.workers,
    }
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"best_g1": best_g1, "best_g2": best_g2}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
