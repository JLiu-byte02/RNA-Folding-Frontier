#!/usr/bin/env python3
"""Predict and score an external secondary-structure baseline on one split.

``predict``
    Run ViennaRNA RNAfold (Turner 2004 MFE) over a split and write one
    ``{"id", "dbn", "mfe"}`` line per sequence. ``--workers`` folds in
    parallel processes; ``--shard N M`` is kept for cluster-array use.

``score``
    Join a prediction file with the label policy used by every other table in
    this project and report family-macro base-pair F1 with a family-block
    bootstrap interval plus sequence-micro F1.

Label policy (must match ``rff.analysis.block_eval.score_dump``):
  * rows whose id is ``rfam:<ACC>|<seed>`` use ``labels_v3``; an explicit
    ``null`` entry (parse failure) is excluded and counted;
  * every other row (bpRNA) uses its own frozen ``pairs`` field, which is
    already canonicalised by ``build_s1prime``;
  * ``--min-true-pairs`` is applied to the **raw** pair set before the nested
    projection, exactly as the dump scorer does;
  * truth is then reduced with ``max_nested_subset`` (Amendment 8).

``--restrict-ids`` limits scoring to the id set of a dump file so the bar is
computed on the same sequences our model is scored on.

Both subcommands refuse any path containing ``test_novel``. The confirmatory
set stays embargoed until Amendment 14 fixes the one-shot analysis list. The
runner is tool-agnostic: other predictors write ``{"id", "dbn"}`` and reuse
``score``.
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402

from rff.analysis.labels import max_nested_subset  # noqa: E402
from rff.data.io import pairs_from_dbn  # noqa: E402

LABELS_V3 = "experiments/011_labels_v3/labels_v3_pairs.json.gz"


def _resolve(path: str) -> Path:
    p = Path(path)
    return p if p.is_absolute() else ROOT / p


def refuse_embargoed(path: Path) -> None:
    if "test_novel" in str(path).lower():
        raise SystemExit(
            f"refusing {path}: test_novel is embargoed until Amendment 14"
        )


def read_jsonl(path: Path) -> list[dict]:
    with open(path, encoding="utf-8") as fh:
        return [json.loads(line) for line in fh if line.strip()]


def load_labels_v3(path: Path) -> dict[str, set[tuple[int, int]] | None]:
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        blob = json.load(fh)
    if isinstance(blob, dict) and "pairs" in blob and isinstance(blob["pairs"], dict):
        blob = blob["pairs"]
    return {
        key: None if pairs is None else {(int(a), int(b)) for a, b in pairs}
        for key, pairs in blob.items()
    }


def f1_from_counts(tp: int, fp: int, fn: int) -> float:
    if tp + fp + fn == 0:
        return 1.0
    if tp == 0:
        return 0.0
    return 2 * tp / (2 * tp + fp + fn)


def family_macro(per_family: dict[str, list[float]]) -> tuple[float, int]:
    if not per_family:
        return (float("nan"), 0)
    means = [float(np.mean(v)) for v in per_family.values() if v]
    return (float(np.mean(means)), len(means))


def family_bootstrap(
    per_family: dict[str, list[float]], n_boot: int, seed: int
) -> tuple[float, float]:
    fam_means = np.asarray(
        [float(np.mean(v)) for v in per_family.values() if v], dtype=float
    )
    if len(fam_means) < 2:
        return (float("nan"), float("nan"))
    rng = np.random.default_rng(seed)
    draws = np.empty(n_boot, dtype=float)
    for k in range(n_boot):
        draws[k] = fam_means[rng.integers(0, len(fam_means), len(fam_means))].mean()
    lo, hi = np.percentile(draws, [2.5, 97.5])
    return (float(lo), float(hi))


def _fold_one(item: tuple[str, str, float | None]) -> tuple[str, str, float]:
    import RNA

    rid, seq, temperature = item
    if temperature is not None:
        RNA.cvar.temperature = temperature
    dbn, mfe = RNA.fold(seq)
    return (rid, dbn, float(mfe))


def cmd_predict(args: argparse.Namespace) -> int:
    split = _resolve(args.split)
    refuse_embargoed(split)
    out = _resolve(args.out)
    rows = read_jsonl(split)

    shard, nshard = args.shard
    if not (0 <= shard < nshard):
        raise SystemExit(f"bad shard {shard}/{nshard}")
    rows = [r for k, r in enumerate(rows) if k % nshard == shard]

    import RNA  # noqa: E402  (validate import in the parent too)

    items = [
        (r["id"], r["seq"].upper().replace("T", "U"), args.temperature) for r in rows
    ]
    out.parent.mkdir(parents=True, exist_ok=True)
    n_written = 0
    with open(out, "w", encoding="utf-8") as fh:
        if args.workers <= 1:
            for rid, dbn, mfe in map(_fold_one, items):
                fh.write(
                    json.dumps(
                        {"id": rid, "dbn": dbn, "mfe": round(mfe, 4)},
                        separators=(",", ":"),
                    )
                    + "\n"
                )
                n_written += 1
        else:
            import multiprocessing as mp

            with mp.Pool(args.workers) as pool:
                for rid, dbn, mfe in pool.imap(_fold_one, items, chunksize=4):
                    fh.write(
                        json.dumps(
                            {"id": rid, "dbn": dbn, "mfe": round(mfe, 4)},
                            separators=(",", ":"),
                        )
                        + "\n"
                    )
                    n_written += 1
                    if n_written % 500 == 0:
                        print(f"[predict] {n_written}/{len(items)}", flush=True)
    meta = {
        "tool": args.tool,
        "split": args.split,
        "n_rows": len(rows),
        "n_written": n_written,
        "shard": [shard, nshard],
        "workers": args.workers,
        "temperature": args.temperature,
        "vienna_version": getattr(RNA, "__version__", None),
    }
    print(json.dumps(meta, indent=2))
    return 0


def cmd_score(args: argparse.Namespace) -> int:
    split_path = _resolve(args.split)
    preds_path = _resolve(args.preds)
    refuse_embargoed(split_path)
    refuse_embargoed(preds_path)

    labels = load_labels_v3(_resolve(args.labels))
    rows = read_jsonl(split_path)
    preds: dict[str, set[tuple[int, int]]] = {}
    for row in read_jsonl(preds_path):
        if "dbn" in row:
            preds[row["id"]] = set(pairs_from_dbn(row["dbn"]))
        elif "pairs" in row:
            preds[row["id"]] = {(int(a), int(b)) for a, b in row["pairs"]}
        else:
            raise SystemExit(f"prediction row for {row.get('id')!r} has neither dbn nor pairs")

    restrict: set[str] | None = None
    if args.restrict_ids:
        restrict = {row["id"] for row in read_jsonl(_resolve(args.restrict_ids))}

    per_family: dict[str, list[float]] = defaultdict(list)
    totals = [0, 0, 0]
    n_used = n_missing_pred = n_lt_min = n_parse_fail = n_empty_truth = 0
    n_not_restricted = 0
    excluded_pred_pairs = 0
    label_source: Counter[str] = Counter()
    n_pred_pairs = 0

    for row in rows:
        rid = row["id"]
        if restrict is not None and rid not in restrict:
            n_not_restricted += 1
            continue
        key = rid[5:] if rid.startswith("rfam:") else rid
        if key in labels:
            truth_raw = labels[key]
            if truth_raw is None:
                n_parse_fail += 1
                label_source["seed_v3_parse_failed"] += 1
                continue
            label_source["seed_v3"] += 1
        else:
            truth_raw = {(int(a), int(b)) for a, b in row["pairs"]}
            label_source["bprna_v1"] += 1

        pred = preds.get(rid)
        if not truth_raw:
            n_empty_truth += 1
        if len(truth_raw) < args.min_true_pairs:
            n_lt_min += 1
            if pred is not None:
                excluded_pred_pairs += len(pred)
            continue
        if pred is None:
            n_missing_pred += 1
            continue
        truth = max_nested_subset(truth_raw)
        tp = len(truth & pred)
        fp = len(pred - truth)
        fn = len(truth - pred)
        totals[0] += tp
        totals[1] += fp
        totals[2] += fn
        n_pred_pairs += len(pred)
        per_family[row["family"]].append(f1_from_counts(tp, fp, fn))
        n_used += 1

    macro, n_fam = family_macro(per_family)
    lo, hi = family_bootstrap(per_family, args.n_boot, args.seed)
    result = {
        "tool": args.tool,
        "split": args.split,
        "preds": args.preds,
        "label_set": args.labels,
        "label_source": dict(label_source),
        "min_true_pairs": args.min_true_pairs,
        "restrict_ids": args.restrict_ids,
        "family_macro_f1": round(macro, 6),
        "family_macro_f1_ci95": [round(lo, 6), round(hi, 6)],
        "n_families": n_fam,
        "micro_f1": round(f1_from_counts(*totals), 6),
        "tp": totals[0],
        "fp": totals[1],
        "fn": totals[2],
        "n_sequences": n_used,
        "mean_pred_pairs": round(n_pred_pairs / n_used, 3) if n_used else None,
        "n_excluded_lt_min_pairs": n_lt_min,
        "n_excluded_empty_truth": n_empty_truth,
        "n_excluded_parse_failure": n_parse_fail,
        "excluded_predicted_pairs": excluded_pred_pairs,
        "n_missing_prediction": n_missing_pred,
        "n_outside_restrict_set": n_not_restricted,
        "n_boot": args.n_boot,
        "bootstrap_seed": args.seed,
    }
    out = _resolve(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(result, indent=2, sort_keys=True) + "\n"
    out.write_text(payload, encoding="utf-8")
    if args.per_family_out:
        pf = _resolve(args.per_family_out)
        pf.parent.mkdir(parents=True, exist_ok=True)
        pf.write_text(
            json.dumps(
                {k: round(float(np.mean(v)), 6) for k, v in per_family.items()},
                indent=2,
                sort_keys=True,
            )
            + "\n",
            encoding="utf-8",
        )
    print(payload)
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("predict")
    p.add_argument("--split", required=True)
    p.add_argument("--tool", default="rnafold")
    p.add_argument("--out", required=True)
    p.add_argument("--shard", type=int, nargs=2, default=[0, 1])
    p.add_argument("--workers", type=int, default=1)
    p.add_argument("--temperature", type=float, default=None)
    p.set_defaults(func=cmd_predict)

    s = sub.add_parser("score")
    s.add_argument("--split", required=True)
    s.add_argument("--preds", required=True)
    s.add_argument("--tool", default="rnafold")
    s.add_argument("--labels", default=LABELS_V3)
    s.add_argument("--restrict-ids", default=None)
    s.add_argument("--out", required=True)
    s.add_argument("--per-family-out", default=None)
    s.add_argument("--min-true-pairs", type=int, default=4)
    s.add_argument("--n-boot", type=int, default=10000)
    s.add_argument("--seed", type=int, default=20261113)
    s.set_defaults(func=cmd_score)

    args = ap.parse_args()
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
