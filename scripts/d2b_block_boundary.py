#!/usr/bin/env python3
"""D2b — H2.3: do extra false positives sit on a block boundary?

Pre-registered wording (Amendment 2 §A2.4):

    H2.3  (c) adds false positives at the window boundary — in P2-block, extra
          FPs under (c) vs (a) that link a block's native window to rewritten
          sequence — passes if >= 60% of the extra FPs.

P2-block embeds each block ``[s, e)`` with native sequence over
``[s - flank, e + flank)`` and donor content elsewhere (``build_block_cache.py``).
For a position ``i`` the *native window* is therefore

    W(i) = [max(0, s_i - flank), min(L, e_i + flank)),   s_i = floor(i / block) * block

A pair ``(i, j)`` is **mutually native** when ``j in W(i)`` and ``i in W(j)``:
neither endpoint's embedding saw rewritten content where its partner sits. Every
other pair is *context-dependent* — at least one endpoint's evidence for the
pair comes from content that was replaced.

**Known degeneracy, stated up front.** ``|W(i)| <= block + 2 * flank``, so any
pair with ``|i - j| > block + 2 * flank`` is context-dependent by construction.
The 60% rule above is therefore not evidence of a *boundary* artefact on its
own: it can pass simply because the extra FPs are long-range, which is H2.1's
claim, not H2.3's. The script reports the pre-registered statistic unchanged and
then stratifies by distance, giving the informative ratio only where the
distinction exists (``d <= block + 2 * flank``) against a distance-matched null
of random canonical, non-true pairs under the same block geometry. Read the
stratified number when judging the mechanism; read the 60% rule when judging the
pre-registration.

Inputs are P2-block dumps from ``e1prime_dump.py`` for conditions (a) and (c).
The script is pure numpy and runs on CPU.
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

CANONICAL = set("AU UA GC CG GU UG".split())


def read_jsonl(path: Path) -> list[dict]:
    with open(path) as fh:
        return [json.loads(line) for line in fh if line.strip()]


def block_window(i: int, length: int, block: int, flank: int) -> tuple[int, int]:
    s = (i // block) * block
    e = min(s + block, length)
    return max(0, s - flank), min(length, e + flank)


def mutually_native(i: int, j: int, length: int, block: int, flank: int) -> bool:
    lo_i, hi_i = block_window(i, length, block, flank)
    lo_j, hi_j = block_window(j, length, block, flank)
    return (lo_i <= j < hi_i) and (lo_j <= i < hi_j)


def distance_band(d: int, span: int) -> str:
    if d <= span // 2:
        return "d<=half_span"
    if d <= span:
        return "half_span<d<=span"
    return "d>span"


def null_pair(
    seq: str,
    d: int,
    truth: set[tuple[int, int]],
    rng: np.random.Generator,
    tries: int = 20,
) -> tuple[int, int] | None:
    if d <= 0 or d >= len(seq):
        return None
    for _ in range(tries):
        i = int(rng.integers(0, len(seq) - d))
        j = i + d
        if (i, j) in truth or (j, i) in truth:
            continue
        if seq[i] + seq[j] in CANONICAL:
            return i, j
    return None


def ratio(sum_a: float, n_a: float, sum_b: float, n_b: float) -> float:
    if n_a == 0 or n_b == 0 or sum_b == 0:
        return float("nan")
    return (sum_a / n_a) / (sum_b / n_b)


def bootstrap_counts(
    per_family: dict[str, dict[str, float]], n_boot: int, seed: int
) -> dict:
    """Family-block bootstrap of observed/null and observed/shared ratios."""
    fams = sorted(f for f in per_family if per_family[f]["extra_n"] > 0)
    if len(fams) < 3:
        return {
            "n_families": len(fams),
            "observed_frac": None,
            "null_frac": None,
            "ratio_obs_null": None,
            "ratio_ci": None,
        }

    def stats(sel):
        agg = {k: 0.0 for k in ("extra_n", "extra_b", "null_n", "null_b")}
        for f in sel:
            row = per_family[f]
            for k in agg:
                agg[k] += row[k]
        obs = agg["extra_b"] / agg["extra_n"] if agg["extra_n"] else float("nan")
        nul = agg["null_b"] / agg["null_n"] if agg["null_n"] else float("nan")
        return obs, nul, ratio(agg["extra_b"], agg["extra_n"], agg["null_b"], agg["null_n"])

    obs, nul, point = stats(fams)
    rng = np.random.default_rng(seed)
    draws = []
    for _ in range(n_boot):
        idx = rng.integers(0, len(fams), len(fams))
        draws.append(stats([fams[k] for k in idx])[2])
    draws = np.array([d for d in draws if np.isfinite(d)])
    return {
        "n_families": len(fams),
        "observed_frac": round(float(obs), 4) if np.isfinite(obs) else None,
        "null_frac": round(float(nul), 4) if np.isfinite(nul) else None,
        "ratio_obs_null": round(float(point), 4) if np.isfinite(point) else None,
        "ratio_ci": [round(float(np.percentile(draws, 2.5)), 4),
                     round(float(np.percentile(draws, 97.5)), 4)]
        if draws.size
        else None,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--a-dump", required=True, help="condition (a) P2-block dump jsonl")
    ap.add_argument("--c-dump", required=True, help="condition (c) P2-block dump jsonl")
    ap.add_argument("--split", required=True, help="jsonl with id -> seq")
    ap.add_argument("--name", default="p2_block_ac")
    ap.add_argument("--block", type=int, default=16)
    ap.add_argument("--flank", type=int, default=32)
    ap.add_argument("--threshold", type=float, default=0.5)
    ap.add_argument("--decode", choices=["bare", "constrained"], default="bare")
    ap.add_argument("--n-boot", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--min-extra", type=int, default=30)
    ap.add_argument("--label-version", choices=("v1", "v3"), default="v1")
    ap.add_argument("--min-true-pairs", type=int, default=0)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    seqs = {r["id"]: r["seq"] for r in read_jsonl(ROOT / args.split)}
    a_rows = {r["id"]: r for r in read_jsonl(Path(args.a_dump))}
    c_rows = {r["id"]: r for r in read_jsonl(Path(args.c_dump))}
    common = sorted(set(a_rows) & set(c_rows))
    if not common:
        raise SystemExit("no shared sequence ids between the two dumps")

    span = args.block + 2 * args.flank
    rng = np.random.default_rng(args.seed)
    per_family: dict[str, dict[str, float]] = defaultdict(
        lambda: {
            "extra_n": 0.0,
            "extra_b": 0.0,
            "null_n": 0.0,
            "null_b": 0.0,
            "shared_n": 0.0,
            "shared_b": 0.0,
            "tp_n": 0.0,
            "tp_b": 0.0,
        }
    )
    bands: dict[str, dict[str, float]] = defaultdict(
        lambda: {"extra_n": 0.0, "extra_b": 0.0, "null_n": 0.0, "null_b": 0.0}
    )
    n_skipped = 0
    n_extra_total = 0
    n_excluded_below_min = 0
    n_excluded_parse = 0
    truth_key = "truth" if args.label_version == "v1" else "truth_v3"

    for rid in common:
        seq = seqs.get(rid)
        if seq is None:
            n_skipped += 1
            continue
        raw_truth = c_rows[rid][truth_key]
        if raw_truth is None:
            n_excluded_parse += 1
            continue
        if len(raw_truth) < args.min_true_pairs:
            n_excluded_below_min += 1
            continue
        truth = {tuple(p) for p in raw_truth}
        fam = c_rows[rid]["family"]

        def pairs_of(row, length: int) -> set[tuple[int, int]]:
            scored = normalise_scored(row["scored"])
            if args.decode == "bare":
                return {k for k, v in scored.items() if v >= args.threshold}
            return set(decode_pairs(scored, args.threshold, length))

        pred_a = pairs_of(a_rows[rid], len(seq))
        pred_c = pairs_of(c_rows[rid], len(seq))
        extra = (pred_c - pred_a) - truth
        shared_fp = (pred_a & pred_c) - truth
        tp_c = truth & pred_c
        n_extra_total += len(extra)

        for i, j in extra:
            d = abs(i - j)
            bad = not mutually_native(i, j, len(seq), args.block, args.flank)
            per_family[fam]["extra_n"] += 1
            per_family[fam]["extra_b"] += float(bad)
            band = distance_band(d, span)
            bands[band]["extra_n"] += 1
            bands[band]["extra_b"] += float(bad)
            pair = null_pair(seq, d, truth, rng)
            if pair is not None:
                nbad = not mutually_native(*pair, len(seq), args.block, args.flank)
                per_family[fam]["null_n"] += 1
                per_family[fam]["null_b"] += float(nbad)
                bands[band]["null_n"] += 1
                bands[band]["null_b"] += float(nbad)

        for i, j in shared_fp:
            per_family[fam]["shared_n"] += 1
            per_family[fam]["shared_b"] += float(
                not mutually_native(i, j, len(seq), args.block, args.flank)
            )
        for i, j in tp_c:
            per_family[fam]["tp_n"] += 1
            per_family[fam]["tp_b"] += float(
                not mutually_native(i, j, len(seq), args.block, args.flank)
            )

    def pooled(key_n: str, key_b: str) -> tuple[float, int]:
        n = sum(v[key_n] for v in per_family.values())
        b = sum(v[key_b] for v in per_family.values())
        return (round(float(b / n), 4) if n else float("nan")), int(n)

    extra_frac, extra_n = pooled("extra_n", "extra_b")
    shared_frac, shared_n = pooled("shared_n", "shared_b")
    tp_frac, tp_n = pooled("tp_n", "tp_b")
    boot = bootstrap_counts(per_family, args.n_boot, args.seed)

    band_out = {}
    for band, row in sorted(bands.items()):
        n, b = row["extra_n"], row["extra_b"]
        nn, nb = row["null_n"], row["null_b"]
        band_out[band] = {
            "extra_n": int(n),
            "extra_boundary_frac": round(float(b / n), 4) if n else None,
            "null_n": int(nn),
            "null_boundary_frac": round(float(nb / nn), 4) if nn else None,
            "ratio_obs_null": round(ratio(b, n, nb, nn), 4)
            if np.isfinite(ratio(b, n, nb, nn))
            else None,
        }

    report = {
        "name": args.name,
        "definition": {
            "block": args.block,
            "flank": args.flank,
            "native_span": span,
            "decoder": args.decode,
            "threshold": args.threshold,
            "mutually_native": "j in W(i) and i in W(j); W(i) = block(i) +/- flank",
        },
        "n_sequences": len(common)
        - n_skipped
        - n_excluded_below_min
        - n_excluded_parse,
        "n_sequences_without_seq": n_skipped,
        "label_policy": {
            "version": args.label_version,
            "min_true_pairs": args.min_true_pairs,
            "n_excluded_below_min": n_excluded_below_min,
            "n_excluded_label_parse_failure": n_excluded_parse,
        },
        "n_extra_fp": extra_n,
        "H2_3": {
            "extra_boundary_frac": extra_frac if extra_n else None,
            "passes_60pct_rule": bool(extra_n and extra_frac >= 0.60),
            "underpowered": bool(extra_n < args.min_extra),
            "null_frac": boot["null_frac"],
            "ratio_obs_null": boot["ratio_obs_null"],
            "ratio_ci": boot["ratio_ci"],
            "n_families": boot["n_families"],
        },
        "reference_rates": {
            "shared_fp": {"frac": shared_frac if shared_n else None, "n": shared_n},
            "true_positive_c": {"frac": tp_frac if tp_n else None, "n": tp_n},
        },
        "by_distance_band": band_out,
        "note": (
            "Pairs with |i-j| > native_span are boundary by construction; the "
            "60% rule can therefore pass for purely long-range extra FPs. Judge "
            "the mechanism on by_distance_band, where null_frac is the "
            "distance-matched expectation."
        ),
    }

    out = ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({args.name: {"extra_n": extra_n,
                                  "frac": report["H2_3"]["extra_boundary_frac"],
                                  "passes": report["H2_3"]["passes_60pct_rule"]}}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
