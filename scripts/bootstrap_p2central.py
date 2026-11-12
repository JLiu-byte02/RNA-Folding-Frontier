#!/usr/bin/env python3
"""Family-bootstrap CIs for the corrected P2-central retention ratios.

The aggregate table in `docs/p2central_fixed_2026-11-12.md` has point estimates
only. Ratios of differences are unstable, so this script recomputes the
per-family values from the p2fix dumps and bootstraps the c/d retention ratios.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402

from rff.analysis.block_eval import (  # noqa: E402
    SeqScore,
    bootstrap_ratio,
    common_ids,
    family_means,
    seed_average_family_means,
)
from rff.analysis.scoring import decode_pairs, normalise_scored  # noqa: E402


def rnd(value: float, nd: int = 6) -> float:
    return round(float(value), nd)


def score_rows(path: Path, kept_only: bool, min_true: int) -> dict[str, SeqScore]:
    out = {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if not line.strip():
                continue
            row = json.loads(line)
            raw = row.get("truth_v3")
            if raw is None or len(raw) < min_true:
                continue
            truth = {tuple(int(x) for x in pair) for pair in raw}
            pred = set(
                decode_pairs(
                    normalise_scored(row["scored"]),
                    0.5,
                    int(row["length"]),
                    partner=True,
                    nesting=True,
                )
            )
            if kept_only:
                lo, hi = row["kept"]
                truth = {p for p in truth if lo <= p[0] < hi and lo <= p[1] < hi}
                pred = {p for p in pred if lo <= p[0] < hi and lo <= p[1] < hi}
                if not truth:
                    continue
            out[row["id"]] = SeqScore(
                id=row["id"],
                family=row["family"],
                tp=len(truth & pred),
                fp=len(pred - truth),
                fn=len(truth - pred),
            )
    return out


def ratio_for(
    per_seed: dict[int, dict[str, dict[str, SeqScore]]],
    numerator: str,
    kept_only: bool,
) -> dict:
    ids_by_seed = {
        seed: common_ids(tables["a"], tables[numerator], tables["e"])
        for seed, tables in per_seed.items()
    }
    ids = set.intersection(*ids_by_seed.values())
    values = {
        cond: seed_average_family_means(
            [per_seed[seed][cond] for seed in per_seed], ids
        )
        for cond in ("a", numerator, "e")
    }
    ratio = bootstrap_ratio(
        values, numerator=numerator, denominator_a="a", denominator_e="e"
    )
    return {
        "n_sequences": len(ids),
        "n_families": ratio["n_families"],
        "R": {k: rnd(v) if isinstance(v, float) else v for k, v in ratio.items()},
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump-root", default="experiments/002_e1prime_runs")
    ap.add_argument("--seeds", nargs="+", type=int, default=[17, 23, 41])
    ap.add_argument("--splits", nargs="+", default=["id_test", "test_withheld"])
    ap.add_argument("--min-true-pairs", type=int, default=4)
    ap.add_argument(
        "--out", default="experiments/002_e1prime_runs/p2central_ratios_ci.json"
    )
    ap.add_argument("--md-out", default="docs/p2central_ratios_ci_2026-11-12.md")
    args = ap.parse_args()

    result = {"threshold": 0.5, "min_true_pairs": args.min_true_pairs, "splits": {}}
    for split in args.splits:
        per_seed = {}
        for seed in args.seeds:
            folder = ROOT / args.dump_root / f"dump_p2fix_seed{seed}"
            per_seed[seed] = {
                cond: score_rows(
                    folder / f"{cond}_{split}.jsonl",
                    kept_only=False,
                    min_true=args.min_true_pairs,
                )
                for cond in ("a", "c", "d", "e")
            }
        result["splits"][split] = {
            "all_pairs": {
                "c": ratio_for(per_seed, "c", kept_only=False),
                "d": ratio_for(per_seed, "d", kept_only=False),
            },
            "kept_kept": {
                "c": ratio_for(
                    {
                        seed: {
                            cond: score_rows(
                                ROOT
                                / args.dump_root
                                / f"dump_p2fix_seed{seed}"
                                / f"{cond}_{split}.jsonl",
                                kept_only=True,
                                min_true=args.min_true_pairs,
                            )
                            for cond in ("a", "c", "d", "e")
                        }
                        for seed in args.seeds
                    },
                    "c",
                    kept_only=True,
                ),
                "d": ratio_for(
                    {
                        seed: {
                            cond: score_rows(
                                ROOT
                                / args.dump_root
                                / f"dump_p2fix_seed{seed}"
                                / f"{cond}_{split}.jsonl",
                                kept_only=True,
                                min_true=args.min_true_pairs,
                            )
                            for cond in ("a", "c", "d", "e")
                        }
                        for seed in args.seeds
                    },
                    "d",
                    kept_only=True,
                ),
            },
        }

    out = ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    md = render_markdown(result)
    md_path = ROOT / args.md_out
    md_path.parent.mkdir(parents=True, exist_ok=True)
    md_path.write_text(md, encoding="utf-8")
    print(json.dumps({"out": args.out, "md": args.md_out}))
    return 0


def render_markdown(result: dict) -> str:
    lines = [
        "# P2-central retention ratios with family-bootstrap CIs",
        "",
        "Computed by `scripts/bootstrap_p2central.py` from the corrected p2fix",
        "dumps; fixed threshold 0.5. The ratio is `(X - e) / (a - e)` on",
        "family-macro F1; CIs are 10,000 family-bootstrap resamples.",
        "",
        "| split | readout | ratio | n | families | point | 95% CI | undefined |",
        "|---|---|---|---:|---:|---:|---|---:|",
    ]
    for split, blocks in result["splits"].items():
        for readout in ("kept_kept", "all_pairs"):
            for cond in ("c", "d"):
                row = blocks[readout][cond]
                r = row["R"]
                lines.append(
                    f"| {split} | {readout} | R_{cond} | {row['n_sequences']} | "
                    f"{row['n_families']} | {r['ratio']:.4f} | "
                    f"[{r['lo']:.4f}, {r['hi']:.4f}] | "
                    f"{r['undefined_fraction']:.4f} |"
                )
    lines += [
        "",
        "The ratio is undefined when `a - e <= 0`; the undefined fraction counts",
        "bootstrap draws in that case. These CIs are diagnostic and do not",
        "replace the P3 D1-core verdict.",
    ]
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    raise SystemExit(main())
