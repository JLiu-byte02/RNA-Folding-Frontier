#!/usr/bin/env python3
"""H2.4 — does exact DP remove long-range false positives faster than true pairs?

Pre-registered wording (Amendment 2 §A2.4):

    H2.4  Global consistency removes them — constrained vs bare-threshold
          decoding (later decoder A) — passes if long-range FPs fall by more
          than true positives at fixed 0.5.

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
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    seqs = {r["id"]: r["seq"] for r in read_jsonl(ROOT / args.split)}
    rows = read_jsonl(Path(args.dump))
    totals = {k: defaultdict(float) for k in ("tp_bare", "lr_bare", "tp_greedy", "lr_greedy", "tp_dp", "lr_dp")}
    delta = defaultdict(lambda: {"d_tp": 0.0, "d_lrfp": 0.0, "n": 0.0})
    n_dp_better = 0
    n_skipped = 0
    for row in rows:
        seq = seqs.get(row["id"])
        if seq is None:
            n_skipped += 1
            continue
        truth = {tuple(p) for p in row["truth"]}
        scored = normalise_scored(row["scored"])
        fam = row["family"]

        bare = {k for k, v in scored.items() if v >= args.threshold}
        greedy = set(decode_pairs(scored, args.threshold, len(seq)))
        model, S = make_model(seq, scored, args)
        dp_pairs, _ = NNDecoder(model).decode(seq, S)
        dp = set(dp_pairs)

        for key, pred in (("bare", bare), ("greedy", greedy), ("dp", dp)):
            tp, lr = counts(pred, truth)
            totals[f"tp_{key}"][fam] += tp
            totals[f"lr_{key}"][fam] += lr

        tp_g, lr_g = counts(greedy, truth)
        tp_d, lr_d = counts(dp, truth)
        delta[fam]["d_tp"] += tp_d - tp_g
        delta[fam]["d_lrfp"] += lr_d - lr_g
        delta[fam]["n"] += 1
        if lr_d - lr_g < 0:
            n_dp_better += 1

    def family_macro(key: str) -> float:
        vals = [v / delta[f]["n"] for f, v in totals[key].items() if delta[f]["n"]]
        return round(float(np.mean(vals)), 4) if vals else float("nan")

    n_seq = sum(int(v["n"]) for v in delta.values())
    report = {
        "name": args.name or Path(args.dump).stem,
        "n_sequences": n_seq,
        "n_sequences_without_seq": n_skipped,
        "threshold": args.threshold,
        "long_range_definition": f"|i-j| > {LONG_RANGE}",
        "dp_config": {
            "hairpin_penalty": args.hp,
            "interior_penalty": args.interior,
            "ml_closing": args.mlc,
            "ml_branch": args.mlb,
            "ml_unpaired": args.mlu,
            "max_loop": args.max_loop,
            "score": "log-odds of the 0.5-thresholded dump cells; absent cells disallowed",
        },
        "family_macro_per_sequence": {
            "bare": {"tp": family_macro("tp_bare"), "long_range_fp": family_macro("lr_bare")},
            "greedy": {"tp": family_macro("tp_greedy"), "long_range_fp": family_macro("lr_greedy")},
            "dp": {"tp": family_macro("tp_dp"), "long_range_fp": family_macro("lr_dp")},
        },
        "H2_4": {
            "mean_d_tp": round(float(np.mean([delta[f]["d_tp"] / delta[f]["n"] for f in delta if delta[f]["n"]])), 4),
            "mean_d_long_range_fp": round(float(np.mean([delta[f]["d_lrfp"] / delta[f]["n"] for f in delta if delta[f]["n"]])), 4),
            "contrast_drop_lrfp_minus_drop_tp": bootstrap(delta, args.n_boot, args.seed),
            "n_sequences_lrfp_reduced": n_dp_better,
            "passes": None,
        },
    }
    ci = report["H2_4"]["contrast_drop_lrfp_minus_drop_tp"]["ci"]
    point = report["H2_4"]["contrast_drop_lrfp_minus_drop_tp"]["contrast"]
    report["H2_4"]["passes"] = bool(point is not None and point > 0 and ci and ci[0] > 0)

    out = ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({report["name"]: report["H2_4"]}, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
