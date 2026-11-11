#!/usr/bin/env python3
"""H2.4 — does exact DP remove long-range false positives faster than true pairs?

Pre-registered wording, **restated in Amendment 3 §A3.1** because the production
decoder already enforced one partner and non-crossing:

    H2.4  Decoder A — exact nearest-neighbour DP over loop energies plus
          learned pair corrections — reduces long-range false positives
          relative to greedy nested matching AT MATCHED RECALL, and the
          reduction is larger than any accompanying loss of true positives.

So there are two readouts, and only the second one confirms or refutes H2.4:

* ``at_0_5`` — both decoders at the fixed 0.5 threshold (the old wording;
  exploratory, and mostly a calibration contrast);
* ``matched_recall`` — a single global greedy threshold ``t_g`` chosen by grid
  search so that pooled micro recall equals decoder A's, within
  ``--recall-tol``. This is the confirmatory readout: it asks whether energetic
  consistency helps once the operating point is taken away as an explanation.

Three decoders are compared on the same probability map, at the same 0.5
threshold:

* ``bare``        every cell >= 0.5 (no constraints at all);
* ``greedy``      the P2 production rule: one partner per position, non-crossing;
* ``dp``          decoder A: exact nearest-neighbour DP over the log-odds map
                  (``rff.decoders.nn_dp``), with structural priors at their
                  documented default of zero.

The zero-penalty DP is the honest default here: any loop penalty would be a
knob that could be tuned until the H2.4 contrast appears, and this run is a
test, not a tuning session. Penalties are exposed as flags so a later run can
state its settings explicitly; the settings used are recorded in the output.

Statistic. Per sequence, ``dTP = TP_dp - TP_greedy`` and
``dLRFP = LRFP_dp - LRFP_greedy`` (LRFP = false positive with ``|i-j| > 128``).
The hypothesis is that LRFP falls more than TP, i.e.
``drop_lrfp = -dLRFP`` exceeds ``drop_tp = -dTP``. The contrast
``drop_lrfp - drop_tp`` is bootstrapped over families, because sequences inside
a family are not independent — the same unit as the P2 table and D1/K1.

Matchability is not guaranteed. Greedy nested decoding accepts pairs in
decreasing probability, so a pair that blocks true positives always outranks
them; removing it by raising the threshold removes the true pairs as well. When
no threshold reaches decoder A's recall, the script reports
``H2_4.status = "inconclusive_no_recall_match"`` rather than inventing a
comparison: the DP recall advantage is then itself the finding, and the LRFP
contrast must be restated by the PI before it is computed.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rff.analysis.scoring import decode_pairs, normalise_scored  # noqa: E402
from rff.decoders.nn_dp import (  # noqa: E402
    NEG,
    NNDecoder,
    PairScoreModel,
    logit_score,
    mask_allowed,
)

LONG_RANGE = 128


def read_jsonl(path: Path) -> list[dict]:
    with open(path) as fh:
        return [json.loads(line) for line in fh if line.strip()]


def dense_scores(scored: dict[tuple[int, int], float], seq: str) -> np.ndarray:
    """Log-odds matrix, ``NEG`` where the dump has no candidate or pairing is illegal."""
    n = len(seq)
    prob = np.zeros((n, n))
    present = np.zeros((n, n), dtype=bool)
    for (i, j), p in scored.items():
        if i < j:
            prob[i, j] = prob[j, i] = p
            present[i, j] = present[j, i] = True
    score = logit_score(prob)
    allowed = mask_allowed(seq)
    return np.where(present & allowed & np.triu(np.ones((n, n), dtype=bool), 1), score, NEG)


def make_model(seq: str, scored, args) -> tuple[PairScoreModel, np.ndarray]:
    S = dense_scores(scored, seq)
    model = PairScoreModel(
        pair_score=S,
        hairpin=np.array([args.hp] * (len(seq) + 2)) if args.hp else None,
        interior=np.full((args.max_loop + 1, args.max_loop + 2), args.interior)
        if args.interior
        else None,
        ml_closing=args.mlc,
        ml_branch=args.mlb,
        ml_unpaired=args.mlu,
        max_loop=args.max_loop,
    )
    return model, S


def counts(pred: set[tuple[int, int]], truth: set[tuple[int, int]]) -> tuple[int, int]:
    tp = len(pred & truth)
    lrfp = sum(1 for i, j in pred - truth if abs(i - j) > LONG_RANGE)
    return tp, lrfp


def bootstrap(delta: dict[str, dict[str, float]], n_boot: int, seed: int) -> dict:
    fams = sorted(f for f in delta if delta[f]["n"] > 0)
    if len(fams) < 3:
        return {"n_families": len(fams), "contrast": None, "ci": None}

    def stat(sel):
        d_tp = np.mean([delta[f]["d_tp"] for f in sel if delta[f]["n"]])
        d_lr = np.mean([delta[f]["d_lrfp"] for f in sel if delta[f]["n"]])
        return float(-d_lr - -d_tp)  # drop_lrfp - drop_tp

    point = stat(fams)
    rng = np.random.default_rng(seed)
    draws = np.array([stat([fams[k] for k in rng.integers(0, len(fams), len(fams))])
                      for _ in range(n_boot)])
    return {
        "n_families": len(fams),
        "contrast": round(point, 4),
        "ci": [round(float(np.percentile(draws, 2.5)), 4),
               round(float(np.percentile(draws, 97.5)), 4)],
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump", required=True, help="per-sequence dump jsonl")
    ap.add_argument("--split", required=True, help="jsonl with id -> seq")
    ap.add_argument("--name", default="")
    ap.add_argument("--threshold", type=float, default=0.5)
    ap.add_argument("--max-loop", type=int, default=30)
    ap.add_argument("--hp", type=float, default=0.0, help="hairpin penalty (score units)")
    ap.add_argument("--interior", type=float, default=0.0, help="interior/bulge penalty")
    ap.add_argument("--mlc", type=float, default=0.0)
    ap.add_argument("--mlb", type=float, default=0.0)
    ap.add_argument("--mlu", type=float, default=0.0)
    ap.add_argument("--n-boot", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument(
        "--energy",
        choices=["struct", "vienna"],
        default="struct",
        help="struct = A-struct (log-odds, zero penalties); vienna = A-energy (Turner + sc_add_bp)",
    )
    ap.add_argument("--energy-scale", type=float, default=1.0,
                    help="sc_add_bp energy per score unit (vienna only; calibrate on the server)")
    ap.add_argument("--recall-tol", type=float, default=0.002)
    ap.add_argument("--threshold-lo", type=float, default=0.05)
    ap.add_argument("--threshold-hi", type=float, default=0.99)
    ap.add_argument("--threshold-grid", type=int, default=24)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    vienna = None
    if args.energy == "vienna":
        try:
            import RNA  # noqa: F401

            from rff.decoders.vienna_soft import ViennaSoftDecoder

            vienna = ViennaSoftDecoder(scale=args.energy_scale)
        except ImportError:
            vienna = None

    seqs = {r["id"]: r["seq"] for r in read_jsonl(ROOT / args.split)}
    rows = read_jsonl(Path(args.dump))
    records = []
    n_skipped = 0
    for row in rows:
        seq = seqs.get(row["id"])
        if seq is None:
            n_skipped += 1
            continue
        truth = {tuple(p) for p in row["truth"]}
        scored = normalise_scored(row["scored"])
        model, S = make_model(seq, scored, args)
        if args.energy == "vienna":
            if vienna is None:
                raise SystemExit(
                    "ViennaRNA is not importable: the A-energy variant cannot run here. "
                    "Install ViennaRNA (or use --energy struct for the exploratory variant)."
                )
            dp_pairs, _ = vienna.decode(seq, S)
        else:
            dp_pairs, _ = NNDecoder(model).decode(seq, S)
        records.append(
            {
                "family": row["family"],
                "truth": truth,
                "scored": scored,
                "length": len(seq),
                "dp": set(dp_pairs),
            }
        )

    n_truth = sum(len(r["truth"]) for r in records)
    n_tp_dp = sum(len(r["dp"] & r["truth"]) for r in records)
    recall_dp = n_tp_dp / n_truth if n_truth else float("nan")

    def evaluate(th: float) -> dict:
        """All decoders at one greedy threshold; DP does not depend on ``th``."""
        per_family = defaultdict(lambda: {"d_tp": 0.0, "d_lrfp": 0.0, "n": 0.0})
        totals = {
            "tp_bare": 0, "lr_bare": 0,
            "tp_greedy": 0, "lr_greedy": 0,
            "tp_dp": 0, "lr_dp": 0,
        }
        per_family_totals = defaultdict(lambda: defaultdict(float))
        n_reduced = 0
        truth_total = 0
        for rec in records:
            fam, truth = rec["family"], rec["truth"]
            bare = {k for k, v in rec["scored"].items() if v >= th}
            greedy = set(decode_pairs(rec["scored"], th, rec["length"]))
            dp = rec["dp"]
            for key, pred in (("bare", bare), ("greedy", greedy), ("dp", dp)):
                tp, lr = counts(pred, truth)
                totals[f"tp_{key}"] += tp
                totals[f"lr_{key}"] += lr
                per_family_totals[fam][f"tp_{key}"] += tp
                per_family_totals[fam][f"lr_{key}"] += lr
            tp_g, lr_g = counts(greedy, truth)
            tp_d, lr_d = counts(dp, truth)
            per_family[fam]["d_tp"] += tp_d - tp_g
            per_family[fam]["d_lrfp"] += lr_d - lr_g
            per_family[fam]["n"] += 1
            truth_total += len(truth)
            if lr_d - lr_g < 0:
                n_reduced += 1

        def family_macro(key: str) -> float:
            vals = [
                per_family_totals[f][key] / per_family[f]["n"]
                for f in per_family_totals
                if per_family[f]["n"]
            ]
            return round(float(np.mean(vals)), 4) if vals else float("nan")

        contrast = bootstrap(per_family, args.n_boot, args.seed)
        ci = contrast["ci"]
        point = contrast["contrast"]
        return {
            "threshold": th,
            "greedy_recall": round(totals["tp_greedy"] / truth_total, 4) if truth_total else None,
            "dp_recall": round(recall_dp, 4),
            "family_macro_per_sequence": {
                "bare": {"tp": family_macro("tp_bare"), "long_range_fp": family_macro("lr_bare")},
                "greedy": {"tp": family_macro("tp_greedy"), "long_range_fp": family_macro("lr_greedy")},
                "dp": {"tp": family_macro("tp_dp"), "long_range_fp": family_macro("lr_dp")},
            },
            "mean_d_tp": round(float(np.mean([
                per_family[f]["d_tp"] / per_family[f]["n"] for f in per_family if per_family[f]["n"]
            ])), 4),
            "mean_d_long_range_fp": round(float(np.mean([
                per_family[f]["d_lrfp"] / per_family[f]["n"] for f in per_family if per_family[f]["n"]
            ])), 4),
            "contrast_drop_lrfp_minus_drop_tp": contrast,
            "n_sequences_lrfp_reduced": n_reduced,
            "passes": bool(point is not None and point > 0 and ci and ci[0] > 0),
        }

    def greedy_recall_at(th: float) -> float:
        tp = 0
        for rec in records:
            pred = set(decode_pairs(rec["scored"], th, rec["length"]))
            tp += len(pred & rec["truth"])
        return tp / n_truth if n_truth else float("nan")

    # Greedy nested decoding is **not** monotone in the threshold: dropping one
    # high-probability false pair can free two true ones, so recall can rise as
    # the threshold rises. A bisection would silently pick the wrong side of
    # such a step, so scan a grid and take the closest achievable recall.
    grid = np.unique(
        np.concatenate(
            [
                np.linspace(args.threshold_lo, args.threshold_hi, args.threshold_grid),
                [args.threshold],
            ]
        )
    )
    best = None
    for th in grid:
        gap = abs(greedy_recall_at(float(th)) - recall_dp)
        if best is None or gap < best[0]:
            best = (gap, float(th))
    gap, t_g = best
    matched = bool(gap <= args.recall_tol)

    b = evaluate(args.threshold)
    report = {
        "name": args.name or Path(args.dump).stem,
        "n_sequences": len(records),
        "n_sequences_without_seq": n_skipped,
        "long_range_definition": f"|i-j| > {LONG_RANGE}",
        "dp_config": {
            "hairpin_penalty": args.hp,
            "interior_penalty": args.interior,
            "ml_closing": args.mlc,
            "ml_branch": args.mlb,
            "ml_unpaired": args.mlu,
            "max_loop": args.max_loop,
            "score": "log-odds of the dumped cells; absent cells disallowed",
            "variant": (
                "A-energy (ViennaRNA Turner + sc_add_bp)"
                if args.energy == "vienna"
                else "A-struct (zero loop penalties) — exploratory, NOT the H2.4 test"
            ),
            "energy_scale": args.energy_scale if args.energy == "vienna" else None,
        },
        "at_0_5": b,
        "matched_recall": {
            **evaluate(t_g),
            "target_recall_dp": round(recall_dp, 4),
            "recall_matched": matched,
            "recall_tolerance": args.recall_tol,
            "recall_gap": round(float(gap), 4),
            "threshold_grid_size": int(len(grid)),
            "note": (
                "confirmatory readout per Amendment 3 A3.1; a single global "
                "greedy threshold chosen by grid search on pooled micro recall "
                "(greedy decoding is not monotone in the threshold)"
            ),
        },
        "H2_4": {
            "primary": "matched_recall",
            "matched_recall_threshold": round(t_g, 4),
            "recall_matched": matched,
            "status": None,
            "passes": None,
        },
    }
    if matched:
        report["H2_4"]["passes"] = report["matched_recall"]["passes"]
        report["H2_4"]["status"] = (
            "confirmed" if report["H2_4"]["passes"] else "not_confirmed"
        )
    else:
        # Greedy nested decoding cannot reach decoder A's recall at any
        # threshold here: the pairs that block its true positives outrank them,
        # so removing the blockers also removes the true pairs. The contrast is
        # then undefined rather than negative.
        report["H2_4"]["status"] = "inconclusive_no_recall_match"

    out = ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({report["name"]: report["H2_4"]}, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
