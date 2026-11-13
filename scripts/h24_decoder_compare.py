#!/usr/bin/env python3
"""H2.4 — does exact DP remove long-range false positives faster than true pairs?

Pre-registered wording, **restated in Amendment 3 §A3.1** because the production
decoder already enforced one partner and non-crossing:

    H2.4  Decoder A — exact nearest-neighbour DP over loop energies plus
          learned pair corrections — reduces long-range false positives
          relative to greedy nested matching AT MATCHED RECALL, and the
          reduction is larger than any accompanying loss of true positives.

So there are three readouts, in the order frozen by Amendment 4 §A4.1:

* ``recall_matched`` — **primary**. A single global greedy threshold chosen by
  grid search so that pooled micro recall equals decoder A's within
  ``--recall-tol``; if no exact match exists, the closest point at or above A's
  recall. Verdict: A must have fewer long-range false positives **and** at least
  as many true positives.
* ``pair_budget`` — **fallback, only when no threshold reaches A's recall**. The
  greedy threshold whose total pair count matches A's (within
  ``--pair-tol-frac``). Verdict: A must have **both** more true positives and
  fewer long-range false positives at equal count.
* ``greedy_recall_ceiling`` — descriptive only: the best greedy recall
  reachable over the grid. Reported always.

``at_0_5`` is also written (both decoders at the fixed 0.5 threshold) but it is
exploratory: it mostly reflects calibration, not consistency.

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
decreasing probability, so a pair that blocks true positives can outrank them;
removing it by raising the threshold removes the true pairs as well. If neither
the recall match nor the pair budget is feasible, the script reports
``H2_4.status = "inconclusive_no_recall_match"`` and no further fallback is
added (A4.1 item 4).
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

from rff.analysis.scoring import decode_pairs, length_band, normalise_scored  # noqa: E402
from rff.analysis.operating_point import Scan, select_operating_point  # noqa: E402
from rff.analysis.embargo import assert_allowed  # noqa: E402
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
                    help="sc_add_bp energy per score unit (vienna only; fit by 5-fold CV per A4.2)")
    ap.add_argument("--energy-offset", type=float, default=0.0,
                    help="constant term inside the pseudo-energy, -(scale*score + offset); "
                         "absorbs the pos_weight-20 shift (A4.2)")
    ap.add_argument("--recall-tol", type=float, default=0.002)
    ap.add_argument("--pair-tol-frac", type=float, default=0.01,
                    help="pair-budget fallback counts as matched within this fraction of "
                         "decoder A's predicted pair count (A4.1 item 2)")
    ap.add_argument("--threshold-lo", type=float, default=0.05)
    ap.add_argument("--threshold-hi", type=float, default=0.99)
    ap.add_argument("--threshold-grid", type=int, default=24)
    ap.add_argument("--label-version", choices=("v1", "v3"), default="v1")
    ap.add_argument("--min-true-pairs", type=int, default=0)
    ap.add_argument("--out", required=True)
    ap.add_argument(
        "--confirmatory",
        action="store_true",
        help="required to read a test_novel path; also verifies the amendment hashes (A6.5)",
    )
    args = ap.parse_args()

    assert_allowed(
        [args.dump, args.split],
        confirmatory=args.confirmatory,
        root=ROOT,
    )

    vienna = None
    if args.energy == "vienna":
        try:
            import RNA  # noqa: F401

            from rff.decoders.vienna_soft import ViennaSoftDecoder

            vienna = ViennaSoftDecoder(scale=args.energy_scale, offset=args.energy_offset)
        except ImportError:
            vienna = None

    seqs = {r["id"]: r["seq"] for r in read_jsonl(ROOT / args.split)}
    rows = read_jsonl(Path(args.dump))
    records = []
    n_skipped = 0
    n_below_min = 0
    n_parse_failed = 0
    truth_key = "truth" if args.label_version == "v1" else "truth_v3"
    for row in rows:
        seq = seqs.get(row["id"])
        if seq is None:
            n_skipped += 1
            continue
        raw_truth = row[truth_key]
        if raw_truth is None:
            n_parse_failed += 1
            continue
        if len(raw_truth) < args.min_true_pairs:
            n_below_min += 1
            continue
        truth = {tuple(p) for p in raw_truth}
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
        per_band = defaultdict(lambda: defaultdict(float))
        per_band_n = defaultdict(int)
        totals = {
            "tp_bare": 0, "lr_bare": 0,
            "tp_greedy": 0, "lr_greedy": 0,
            "tp_dp": 0, "lr_dp": 0,
            "n_bare_pairs": 0, "n_greedy_pairs": 0, "n_dp_pairs": 0,
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
                totals[f"n_{key}_pairs"] += len(pred)
                per_family_totals[fam][f"tp_{key}"] += tp
                per_family_totals[fam][f"lr_{key}"] += lr
                band = length_band(rec["length"])
                per_band[band][f"tp_{key}"] += tp
                per_band[band][f"lr_{key}"] += lr
                per_band[band][f"n_{key}_pairs"] += len(pred)
            per_band_n[length_band(rec["length"])] += 1
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
        by_band = {
            band: {
                "n": per_band_n[band],
                **{
                    key: {
                        "tp": per_band[band][f"tp_{key}"],
                        "long_range_fp": per_band[band][f"lr_{key}"],
                        "pairs": per_band[band][f"n_{key}_pairs"],
                    }
                    for key in ("bare", "greedy", "dp")
                },
            }
            for band in sorted(per_band)
        }
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
            "pooled": {
                key: {
                    "tp": totals[f"tp_{key}"],
                    "long_range_fp": totals[f"lr_{key}"],
                    "pairs": totals[f"n_{key}_pairs"],
                }
                for key in ("bare", "greedy", "dp")
            },
            "by_band": by_band,
            "contrast_drop_lrfp_minus_drop_tp": contrast,
            "n_sequences_lrfp_reduced": n_reduced,
            "passes": bool(point is not None and point > 0 and ci and ci[0] > 0),
        }

    def greedy_stats(th: float) -> tuple[float, int]:
        tp = 0
        n_pairs = 0
        for rec in records:
            pred = set(decode_pairs(rec["scored"], th, rec["length"]))
            tp += len(pred & rec["truth"])
            n_pairs += len(pred)
        return (tp / n_truth if n_truth else float("nan")), n_pairs

    # Greedy nested decoding need not be monotone in the threshold: dropping one
    # high-probability false pair can free two true ones, so recall can rise as
    # the threshold rises. A bisection would silently pick the wrong side of
    # such a step, so scan a grid and take the closest achievable point.
    grid = np.unique(
        np.concatenate(
            [
                np.linspace(args.threshold_lo, args.threshold_hi, args.threshold_grid),
                [args.threshold],
            ]
        )
    )
    dp_pairs_total = sum(len(r["dp"]) for r in records)
    scan = [Scan(float(th), *greedy_stats(float(th))) for th in grid]
    ceiling = max(s.recall for s in scan)
    t_ceiling = max(scan, key=lambda s: s.recall).threshold

    picked = select_operating_point(
        scan,
        dp_recall=recall_dp,
        dp_pairs=dp_pairs_total,
        recall_tol=args.recall_tol,
        pair_tol_frac=args.pair_tol_frac,
    )
    primary = picked["recall_matched"]
    pb = picked["pair_budget"]
    pb_th, pb_recall, pb_pairs = pb["threshold"], pb["recall"], pb["pairs"]
    pair_tol = pb["tolerance_pairs"]
    pair_gap = pb["gap_pairs"]
    pair_feasible = pb["feasible"]

    b = evaluate(args.threshold)
    verdicts = {}
    if primary is not None:
        m = evaluate(primary["threshold"])
        g, d = m["pooled"]["greedy"], m["pooled"]["dp"]
        verdicts["recall_matched"] = {
            "pooled_tp_greedy": g["tp"],
            "pooled_tp_dp": d["tp"],
            "pooled_long_range_fp_greedy": g["long_range_fp"],
            "pooled_long_range_fp_dp": d["long_range_fp"],
            "fewer_long_range_fp": bool(d["long_range_fp"] < g["long_range_fp"]),
            "no_fewer_true_positives": bool(d["tp"] >= g["tp"]),
            "passes": bool(d["long_range_fp"] < g["long_range_fp"] and d["tp"] >= g["tp"]),
            "by_band": m["by_band"],
        }
    if pair_feasible:
        p = evaluate(pb_th)
        g, d = p["pooled"]["greedy"], p["pooled"]["dp"]
        verdicts["pair_budget"] = {
            "pooled_pairs_greedy": g["pairs"],
            "pooled_pairs_dp": d["pairs"],
            "pooled_tp_greedy": g["tp"],
            "pooled_tp_dp": d["tp"],
            "pooled_long_range_fp_greedy": g["long_range_fp"],
            "pooled_long_range_fp_dp": d["long_range_fp"],
            "more_true_positives": bool(d["tp"] > g["tp"]),
            "fewer_long_range_fp": bool(d["long_range_fp"] < g["long_range_fp"]),
            "passes": bool(d["tp"] > g["tp"] and d["long_range_fp"] < g["long_range_fp"]),
            "by_band": p["by_band"],
        }

    if primary is not None:
        basis = "recall_matched"
    elif pair_feasible:
        basis = "pair_budget"
    else:
        basis = None

    report = {
        "name": args.name or Path(args.dump).stem,
        "n_sequences": len(records),
        "n_sequences_without_seq": n_skipped,
        "label_policy": {
            "version": args.label_version,
            "min_true_pairs": args.min_true_pairs,
            "n_excluded_below_min": n_below_min,
            "n_excluded_label_parse_failure": n_parse_failed,
        },
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
            "energy_offset": args.energy_offset if args.energy == "vienna" else None,
            "calibration_status": (
                "uncalibrated placeholder (A4.2 requires 5-fold CV fit on the frozen train split)"
                if args.energy == "vienna"
                and args.energy_scale == 1.0
                and args.energy_offset == 0.0
                else "explicit values supplied on the command line; record the fitting job"
                if args.energy == "vienna"
                else None
            ),
        },
        "at_0_5": b,
        "greedy_recall_ceiling": {
            "recall": round(float(ceiling), 4),
            "threshold": round(float(t_ceiling), 4),
            "note": "descriptive only (Amendment 4 A4.1 item 3): the best greedy recall "
                    "reachable on this split, over the whole threshold grid",
        },
        "H2_4": {
            "order": ["recall_matched", "pair_budget"],
            "basis": basis,
            "status": None,
            "passes": None,
            "verdicts": verdicts,
            "recall_matched_threshold": round(primary["threshold"], 4) if primary else None,
            "recall_matched_gap": round(float(primary["gap"]), 4) if primary else None,
            "recall_tolerance": args.recall_tol,
            "pair_budget_threshold": round(float(pb_th), 4),
            "pair_budget_gap_pairs": int(pair_gap),
            "pair_budget_tolerance_pairs": float(pair_tol),
            "pair_budget_feasible": pair_feasible,
            "threshold_grid_size": int(len(grid)),
            "note": (
                "Amendment 4 A4.1: recall-matched primary; pair-budget fallback only "
                "when no threshold reaches decoder A's recall; greedy recall ceiling "
                "always reported, descriptively"
            ),
        },
    }
    if basis is not None:
        report["H2_4"]["passes"] = verdicts[basis]["passes"]
        report["H2_4"]["status"] = "confirmed" if report["H2_4"]["passes"] else "not_confirmed"
    else:
        # Neither a recall-matched nor a pair-budget-matched greedy operating
        # point exists on this split. A4.1 item 4: report it, do not invent a
        # third fallback after seeing data.
        report["H2_4"]["status"] = "inconclusive_no_recall_match"

    out = ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({report["name"]: report["H2_4"]}, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
