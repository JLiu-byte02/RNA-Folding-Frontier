#!/usr/bin/env python3
"""Score the frozen E-C2 coverage gate on a held-out split (Amendment 13 A13.2).

The fit in ``experiments/002_e1prime_runs/ec2_fit.json`` selects gate parameters
on the training split only; this script is the confirmatory read. The same text
dump is decoded with the same soft decoder under three pre-registered arms:

    fixed   s = --fixed-scale                      (A-energy-fixed comparator)
    g1      s = s_max * sigmoid(alpha - beta*z)    (frozen best_g1)
    g2      s = s_max * sigmoid(alpha - beta*z - gamma*log L)

and we report family-macro F1 with a family-block bootstrap interval plus
sequence-micro F1. Deltas are paired at the family level, so a sequence never
contributes to the two sides of a comparison with different weights.

Label policy is the one used by every other table (``score_external_baseline``):
seed rows use ``labels_v3``, a ``null`` entry is a parse failure and is excluded,
bpRNA rows use their frozen ``pairs``, ``--min-true-pairs`` applies to the raw
pair set, and truth is then reduced with ``max_nested_subset`` (Amendment 8/9).

``test_novel`` is refused. No predictions are written; only counts.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import multiprocessing
import os
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from rff.analysis.labels import max_nested_subset  # noqa: E402
from rff.analysis.scoring import normalise_scored  # noqa: E402
from rff.data.io import pairs_from_dbn  # noqa: E402
from fit_ec2_gate import dense_scores, gate_scale, sha256_file  # noqa: E402

LABELS_V3 = "experiments/011_labels_v3/labels_v3_pairs.json.gz"


def resolve(path: str) -> Path:
    p = Path(path)
    return p if p.is_absolute() else ROOT / p


def refuse_embargoed(path: Path) -> None:
    if "test_novel" in str(path).lower():
        raise SystemExit(f"refusing {path}: test_novel is embargoed")


def read_jsonl(path: Path) -> list[dict]:
    with open(path, encoding="utf-8") as fh:
        return [json.loads(line) for line in fh if line.strip()]


def load_labels(path: Path) -> dict[str, set[tuple[int, int]] | None]:
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        blob = json.load(fh)
    if isinstance(blob, dict) and "pairs" in blob and isinstance(blob["pairs"], dict):
        blob = blob["pairs"]
    return {
        k: None if v is None else {(int(a), int(b)) for a, b in v}
        for k, v in blob.items()
    }


def f1_from_counts(tp: int, fp: int, fn: int) -> float:
    if tp + fp + fn == 0:
        return 1.0
    if tp == 0:
        return 0.0
    return 2 * tp / (2 * tp + fp + fn)


def family_means(per_family: dict[str, list[float]]) -> tuple[np.ndarray, list[str]]:
    keys = sorted(per_family)
    return np.asarray([np.mean(per_family[k]) for k in keys], dtype=float), keys


def bootstrap_ci(values: np.ndarray, n_boot: int, seed: int) -> tuple[float, float]:
    if values.size < 2:
        return (float("nan"), float("nan"))
    rng = np.random.default_rng(seed)
    draws = values[rng.integers(0, values.size, size=(n_boot, values.size))].mean(axis=1)
    lo, hi = np.percentile(draws, [2.5, 97.5])
    return (float(lo), float(hi))


def paired_delta_ci(
    a: np.ndarray, b: np.ndarray, n_boot: int, seed: int
) -> tuple[float, float, float]:
    """Mean(a - b) over families with a family-block bootstrap interval."""
    diff = a - b
    lo, hi = bootstrap_ci(diff, n_boot, seed)
    return (float(diff.mean()), lo, hi)


_ITEMS: list[dict] = []
_ARMS: list[float] = []


def _init(items: list[dict], scales: list[float]) -> dict[str, int]:
    """Bind one decode scale per item and report the planned loop cost.

    ``scales`` is per *sequence*, not per arm: item ``i`` is decoded exactly
    once, with ``scales[i]``.  The returned plan lets the caller assert
    linearity instead of silently decoding ``n_items`` times per item -- the
    O(n^2) stall that burned 56 CPU-h on job 15627 (the smoke split, n=298,
    passed because the squared cost was still small).
    """
    global _ITEMS, _ARMS
    if len(scales) != len(items):
        raise SystemExit(
            f"scale list has {len(scales)} entries for {len(items)} items; "
            "the scoring loop requires exactly one scale per sequence"
        )
    _ITEMS = items
    _ARMS = scales
    return {"n_items": len(items), "decodes_per_item": 1, "total_decodes": len(items)}


def _decode_item(idx: int) -> tuple[int, np.ndarray]:
    from rff.decoders.vienna_soft import ViennaSoftDecoder

    item = _ITEMS[idx]
    seq, truth = item["seq"], item["truth"]
    scores = dense_scores(item["scored"], seq)
    pred = set(
        ViennaSoftDecoder(scale=float(_ARMS[idx]), offset=0.0).decode(seq, scores)[0]
    )
    tp = len(pred & truth)
    out = np.zeros((1, 4), dtype=np.int32)
    out[0, 0] = tp
    out[0, 1] = len(pred) - tp
    out[0, 2] = len(truth) - tp
    out[0, 3] = len(pred)
    return idx, out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump", required=True)
    ap.add_argument("--split", required=True)
    ap.add_argument("--pppl", required=True)
    ap.add_argument("--fit", default="experiments/002_e1prime_runs/ec2_fit.json")
    ap.add_argument("--labels", default=LABELS_V3)
    ap.add_argument("--fixed-scale", type=float, default=0.1875)
    ap.add_argument("--min-true-pairs", type=int, default=4)
    ap.add_argument("--n-boot", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=17)
    ap.add_argument("--ref-preds", default="", help="optional rival prediction JSONL (id, dbn)")
    ap.add_argument("--baseline-macro", type=float, default=float("nan"))
    ap.add_argument("--min-gain", type=float, default=0.02)
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    dump_path, split_path = resolve(args.dump), resolve(args.split)
    pppl_path, fit_path = resolve(args.pppl), resolve(args.fit)
    for p in (dump_path, split_path, pppl_path, fit_path):
        refuse_embargoed(p)

    fit = json.loads(fit_path.read_text(encoding="utf-8"))
    if fit.get("status", "").startswith("fitted") is False:
        raise SystemExit(f"fit at {fit_path} is not a frozen selector fit")
    mean = float(fit["train_pppl_mean"])
    sd = float(fit["train_pppl_sd"])
    g1, g2 = fit["best_g1"], fit["best_g2"]

    pppl = json.loads(pppl_path.read_text(encoding="utf-8"))["values"]
    labels = load_labels(resolve(args.labels))
    split_rows = read_jsonl(split_path)
    seq_of = {r["id"]: r for r in split_rows}

    items: list[dict] = []
    dropped = {"no_split_row": 0, "no_pppl": 0, "parse_failed": 0, "below_min_pairs": 0}
    for row in read_jsonl(dump_path):
        rid = row["id"]
        src = seq_of.get(rid)
        if src is None:
            dropped["no_split_row"] += 1
            continue
        rec = pppl.get(rid)
        if rec is None:
            dropped["no_pppl"] += 1
            continue
        key = rid[5:] if rid.startswith("rfam:") else rid
        if key in labels:
            raw = labels[key]
            if raw is None:
                dropped["parse_failed"] += 1
                continue
        else:
            raw = {(int(a), int(b)) for a, b in src["pairs"]}
        if len(raw) < args.min_true_pairs:
            dropped["below_min_pairs"] += 1
            continue
        seq = src["seq"].upper().replace("T", "U")
        z = (rec["pppl"] - mean) / sd if sd > 0 else 0.0
        items.append(
            {
                "id": rid,
                "family": src["family"],
                "seq": seq,
                "scored": normalise_scored(row["scored"]),
                "truth": max_nested_subset(raw),
                "z": z,
                "log_len": float(np.log(len(seq))),
            }
        )
    if not items:
        raise SystemExit("no usable records")

    arm_names = ["fixed", "g1", "g2"]
    scales = {
        "fixed": [args.fixed_scale] * len(items),
        "g1": [gate_scale(g1, it["z"], it["log_len"]) for it in items],
        "g2": [gate_scale(g2, it["z"], it["log_len"]) for it in items],
    }
    # Decode one arm at a time so peak memory stays at len(items) x 4 ints.
    counts: dict[str, np.ndarray] = {}
    ref_ok = bool(args.ref_preds)
    ref_counts = np.zeros((len(items), 4), dtype=np.int32)
    for name in arm_names:
        # Set the globals in the parent and let fork() copy them into each
        # worker.  Passing them as Pool initargs would pickle the multi-GB item
        # list once per worker.
        plan = _init(items, scales[name])
        print(
            f"[ec2] arm {name}: {plan['total_decodes']} decodes "
            f"({plan['n_items']} items x {plan['decodes_per_item']})",
            flush=True,
        )
        out = np.zeros((len(items), 4), dtype=np.int32)
        with multiprocessing.get_context("fork").Pool(args.workers) as ex:
            for idx, c in ex.imap_unordered(_decode_item, range(len(items)), chunksize=16):
                # _ARMS has exactly one entry per call, so each item returns a
                # single 4-count row (tp, fp, fn, npred).
                out[idx] = c[0]
        counts[name] = out
        print(f"[ec2] decoded arm {name} (n={len(items)})", flush=True)

    if ref_ok:
        preds = {r["id"]: r["dbn"] for r in read_jsonl(resolve(args.ref_preds))}
        n_missing = 0
        for i, it in enumerate(items):
            dbn = preds.get(it["id"])
            if dbn is None:
                n_missing += 1
                continue
            pred = set(pairs_from_dbn(dbn))
            tp = len(pred & it["truth"])
            ref_counts[i] = (tp, len(pred) - tp, len(it["truth"]) - tp, len(pred))
        print(f"[ec2] reference preds: missing={n_missing}", flush=True)

    families = [it["family"] for it in items]
    fam_set = sorted(set(families))
    fam_index = {f: k for k, f in enumerate(fam_set)}

    def aggregate(c: np.ndarray) -> dict:
        per: dict[str, list[float]] = defaultdict(list)
        totals = np.zeros(3, dtype=np.int64)
        n_pred = n_empty = 0
        for i, fam in enumerate(families):
            tp, fp, fn, npred = (int(v) for v in c[i])
            per[fam].append(f1_from_counts(tp, fp, fn))
            totals += (tp, fp, fn)
            n_pred += npred
            n_empty += int(npred == 0)
        means, keys = family_means(per)
        lo, hi = bootstrap_ci(means, args.n_boot, args.seed)
        tp, fp, fn = (int(v) for v in totals)
        return {
            "family_macro_f1": round(float(means.mean()), 4),
            "macro_ci95": [round(lo, 4), round(hi, 4)],
            "micro_f1": round(f1_from_counts(tp, fp, fn), 4),
            "n_families": len(keys),
            "n_sequences": len(items),
            "tp": tp,
            "fp": fp,
            "fn": fn,
            "mean_pred_pairs": round(n_pred / len(items), 2),
            "n_empty_pred": n_empty,
            "_means": {k: means[i] for i, k in enumerate(keys)},
        }

    results = {name: aggregate(counts[name]) for name in arm_names}
    if ref_ok:
        results["reference"] = aggregate(ref_counts)
    means_by_arm = {k: v.pop("_means") for k, v in results.items()}

    def means_array(d: dict) -> np.ndarray:
        return np.asarray([d[f] for f in fam_set], dtype=float)

    deltas = {}
    for name in arm_names[1:]:
        d, lo, hi = paired_delta_ci(
            means_array(means_by_arm[name]),
            means_array(means_by_arm["fixed"]),
            args.n_boot,
            args.seed,
        )
        deltas[f"{name}_minus_fixed"] = {
            "delta": round(d, 4),
            "ci95": [round(lo, 4), round(hi, 4)],
            "paired_over_families": True,
        }
    if ref_ok:
        for name in arm_names:
            d, lo, hi = paired_delta_ci(
                means_array(means_by_arm[name]),
                means_array(means_by_arm["reference"]),
                args.n_boot,
                args.seed,
            )
            deltas[f"{name}_minus_reference"] = {
                "delta": round(d, 4),
                "ci95": [round(lo, 4), round(hi, 4)],
                "paired_over_families": True,
            }

    baseline = args.baseline_macro
    if ref_ok:
        baseline = results["reference"]["family_macro_f1"]
    verdict = {}
    for name in arm_names[1:]:
        macro = results[name]["family_macro_f1"]
        d_fixed = deltas[f"{name}_minus_fixed"]
        row = {
            "macro_f1": macro,
            "gain_over_fixed": d_fixed["delta"],
            "gain_over_fixed_ci": d_fixed["ci95"],
            "needs_gain": args.min_gain,
            "beats_fixed_by_min_gain": bool(d_fixed["delta"] >= args.min_gain),
            "fixed_gain_ci_excludes_zero": bool(d_fixed["ci95"][0] > 0),
        }
        if np.isfinite(baseline):
            key = f"{name}_minus_reference" if ref_ok else None
            if key:
                row["baseline_macro"] = baseline
                row["gain_over_baseline"] = deltas[key]["delta"]
                row["gain_over_baseline_ci"] = deltas[key]["ci95"]
                row["beats_baseline_by_min_gain"] = bool(deltas[key]["delta"] >= args.min_gain)
                row["baseline_gain_ci_excludes_zero"] = bool(deltas[key]["ci95"][0] > 0)
            else:
                row["baseline_macro"] = baseline
                row["gain_over_baseline"] = round(macro - baseline, 4)
        row["passes"] = bool(
            row.get("beats_fixed_by_min_gain")
            and row.get("fixed_gain_ci_excludes_zero")
            and row.get("beats_baseline_by_min_gain", True)
            and row.get("baseline_gain_ci_excludes_zero", True)
        )
        verdict[name] = row

    def s_summary(seq_scales: list[float]) -> dict:
        arr = np.asarray(seq_scales, dtype=float)
        return {
            "min": round(float(arr.min()), 6),
            "mean": round(float(arr.mean()), 6),
            "max": round(float(arr.max()), 6),
            "sd": round(float(arr.std()), 6),
            "n_unique_rounded": int(np.unique(np.round(arr, 5)).size),
        }

    out = {
        "protocol": "E-C2 coverage gate, frozen selector fit, withheld read (A13.2)",
        "split": args.split,
        "split_sha256": sha256_file(split_path),
        "dump": args.dump,
        "dump_sha256": sha256_file(dump_path),
        "pppl": args.pppl,
        "pppl_sha256": sha256_file(pppl_path),
        "fit": args.fit,
        "fit_sha256": sha256_file(fit_path),
        "labels": args.labels,
        "label_policy": {
            "version": "v3",
            "truth_mode": "nested",
            "min_true_pairs": args.min_true_pairs,
        },
        "n_items": len(items),
        "n_dropped": dropped,
        "arms": {
            "fixed": {"scale": args.fixed_scale, "scale_sd": 0.0},
            "g1": dict(g1),
            "g2": dict(g2),
        },
        "scale_summary": {
            "g1": s_summary(scales["g1"]),
            "g2": s_summary(scales["g2"]),
        },
        "pppl_uses": {
            "g1_beta": g1["beta"],
            "g2_beta": g2["beta"],
            "note": (
                "beta = 0 means the frozen gate does not read the confidence "
                "signal at all; g1 then collapses to a constant scale and g2 "
                "varies with length only. The grid only allowed beta >= 0, so "
                "a negative dependence is untested, not excluded."
            ),
        },
        "results": results,
        "deltas": deltas,
        "verdict": verdict,
        "baseline_macro_used": None if not np.isfinite(baseline) else round(float(baseline), 4),
        "n_boot": args.n_boot,
        "seed": args.seed,
    }
    dst = resolve(args.out)
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text(json.dumps(out, indent=1) + "\n", encoding="utf-8")
    print(json.dumps({k: out[k] for k in ("n_items", "results", "deltas", "verdict")}, indent=1)[:4000])
    print(f"[ec2] wrote {dst}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())