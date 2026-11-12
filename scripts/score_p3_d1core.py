#!/usr/bin/env python3
"""Score the P3 train-and-test arms and apply the Amendment 10 D1-core rule."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402

from rff.analysis.block_eval import (  # noqa: E402
    bootstrap_mean,
    common_ids,
    family_macro,
    family_means,
    micro,
    score_dump,
    seed_average_family_means,
)

ARMS = ("a", "e", "b")


def rnd(value: float, nd: int = 6) -> float:
    return round(float(value), nd)


def load_arm(
    root: Path, seed: int, arm: str, split: str, nested: bool
) -> dict:
    if arm == "a":
        path = root / f"dump_p2fix_seed{seed}" / f"a_{split}.jsonl"
    else:
        path = root / f"dump_p3_{arm}_seed{seed}" / f"{split}.jsonl"
    if not path.exists():
        raise SystemExit(f"missing dump: {path}")
    return score_dump(path, nested=nested)


def arm_summary(per_seed: dict[int, dict[str, dict]], arm: str) -> dict:
    ids = common_ids(*[per_seed[s][arm] for s in per_seed])
    means = seed_average_family_means([per_seed[s][arm] for s in per_seed], ids)
    boot = bootstrap_mean(means)
    rows = []
    for seed in sorted(per_seed):
        table = per_seed[seed][arm]
        rows.append(
            {
                "seed": seed,
                "family_macro_f1": rnd(family_macro(table, ids)),
                "micro_f1": rnd(micro(table, ids)),
            }
        )
    return {
        "n_sequences": len(ids),
        "n_families": len(means),
        "seed_averaged_family_macro_f1": rnd(boot["mean"]),
        "family_macro_ci95": [rnd(boot["lo"]), rnd(boot["hi"])],
        "per_seed": rows,
        "mean_micro_f1": rnd(np.mean([r["micro_f1"] for r in rows])),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", nargs="+", type=int, default=[17, 23, 41])
    ap.add_argument("--split", default="test_withheld")
    ap.add_argument("--dump-root", default="experiments/002_e1prime_runs")
    ap.add_argument("--nested-report", action="store_true")
    ap.add_argument(
        "--out", default="experiments/002_e1prime_runs/d1core_test_withheld.json"
    )
    ap.add_argument("--md-out", default="docs/d1core_2026-11-12.md")
    args = ap.parse_args()
    dump_root = ROOT / args.dump_root

    per_seed_allpairs = {}
    for seed in args.seeds:
        per_seed_allpairs[seed] = {
            arm: load_arm(dump_root, seed, arm, args.split, nested=False)
            for arm in ARMS
        }

    result: dict = {
        "amendment": "AMENDMENT_10.md",
        "split": args.split,
        "truth": "v3 all-pairs, >=4 true pairs, threshold 0.5",
        "seeds": args.seeds,
        "arms": {
            "a": "native full context (p2fix)",
            "e": "e-block trained/tested (P3)",
            "b": "b-block trained/tested (P3)",
        },
    }
    result["unpaired_allpairs"] = {
        arm: arm_summary(per_seed_allpairs, arm) for arm in ARMS
    }

    # Paired family deltas: each seed's e/b F1 minus the same seed's a F1.
    per_seed_delta: dict[str, dict[int, dict[str, float]]] = {arm: {} for arm in ("e", "b")}
    for seed in args.seeds:
        for arm in ("e", "b"):
            ids = common_ids(
                per_seed_allpairs[seed][arm], per_seed_allpairs[seed]["a"]
            )
            f_block = family_means(per_seed_allpairs[seed][arm], ids)
            f_a = family_means(per_seed_allpairs[seed]["a"], ids)
            per_seed_delta[arm][seed] = {
                family: f_block[family] - f_a[family]
                for family in set(f_block) & set(f_a)
            }

    d1: dict = {}
    for arm in ("e", "b"):
        families = set(per_seed_delta[arm][args.seeds[0]])
        for seed in args.seeds[1:]:
            families &= set(per_seed_delta[arm][seed])
        families = sorted(families)
        averaged = {
            family: float(
                np.mean([per_seed_delta[arm][seed][family] for seed in args.seeds])
            )
            for family in families
        }
        boot = bootstrap_mean(averaged)
        per_seed_rows = []
        for seed in args.seeds:
            vals = [per_seed_delta[arm][seed][f] for f in families]
            per_seed_rows.append(
                {
                    "seed": seed,
                    "mean_delta": rnd(np.mean(vals)),
                    "positive_families": int(sum(v > 0 for v in vals)),
                    "n_families": len(vals),
                }
            )
        d1[arm] = {
            "n_families": len(families),
            "seed_averaged_delta": rnd(boot["mean"]),
            "delta_ci95": [rnd(boot["lo"]), rnd(boot["hi"])],
            "per_seed": per_seed_rows,
            "passes_rule": bool(boot["mean"] >= 0.02 and boot["lo"] > 0),
        }
    result["d1_core"] = {
        "rule": (
            "supported iff at least one block arm has seed-averaged "
            "Delta >= .02 and a 95% family-bootstrap CI lower bound > 0"
        ),
        "arms": d1,
        "supported": bool(d1["e"]["passes_rule"] or d1["b"]["passes_rule"]),
        "falsified": bool(not (d1["e"]["passes_rule"] or d1["b"]["passes_rule"])),
    }

    if args.nested_report:
        nested = {}
        for seed in args.seeds:
            nested[seed] = {
                arm: load_arm(dump_root, seed, arm, args.split, nested=True)
                for arm in ARMS
            }
        result["nested_primary_secondary"] = {
            arm: arm_summary(nested, arm) for arm in ARMS
        }

    result["invariants"] = {}
    gaps = []
    for arm in ARMS:
        for row in result["unpaired_allpairs"][arm]["per_seed"]:
            gaps.append(abs(row["family_macro_f1"] - row["micro_f1"]))
    result["invariants"]["max_macro_micro_gap"] = rnd(max(gaps))
    result["invariants"]["macro_micro_gate_0.15"] = bool(max(gaps) <= 0.15)

    out = ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    md_path = ROOT / args.md_out
    md_path.parent.mkdir(parents=True, exist_ok=True)
    md_path.write_text(render_markdown(result), encoding="utf-8")
    print(
        json.dumps(
            {
                "out": args.out,
                "md": args.md_out,
                "d1_core_falsified": result["d1_core"]["falsified"],
                "d1_core_supported": result["d1_core"]["supported"],
                "e": result["d1_core"]["arms"]["e"],
                "b": result["d1_core"]["arms"]["b"],
            }
        )
    )
    return 0


def render_markdown(result: dict) -> str:
    d1 = result["d1_core"]
    lines = [
        "# P3 train-and-test and the D1-core verdict",
        "",
        "Generated by `scripts/score_p3_d1core.py` from the dumps.",
        "",
        f"Split: `{result['split']}`; {result['truth']}.",
        "",
        "## Unpaired all-pairs family-macro F1",
        "",
        "| arm | n | families | family-macro F1 | 95% CI | mean micro F1 |",
        "|---|---:|---:|---:|---|---:|",
    ]
    for arm in ARMS:
        row = result["unpaired_allpairs"][arm]
        lines.append(
            f"| {arm} ({result['arms'][arm]}) | {row['n_sequences']} | "
            f"{row['n_families']} | {row['seed_averaged_family_macro_f1']:.4f} | "
            f"[{row['family_macro_ci95'][0]:.4f}, {row['family_macro_ci95'][1]:.4f}] | "
            f"{row['mean_micro_f1']:.4f} |"
        )
    lines += [
        "",
        "## D1-core paired deltas",
        "",
        "| arm | seed-averaged Delta | 95% CI | passes | per-seed mean Delta |",
        "|---|---:|---|---|---|",
    ]
    for arm in ("e", "b"):
        row = d1["arms"][arm]
        per_seed = "; ".join(
            f"{x['seed']}: {x['mean_delta']:+.4f}" for x in row["per_seed"]
        )
        lines.append(
            f"| {arm} | {row['seed_averaged_delta']:+.4f} | "
            f"[{row['delta_ci95'][0]:+.4f}, {row['delta_ci95'][1]:+.4f}] | "
            f"{row['passes_rule']} | {per_seed} |"
        )
    lines += [
        "",
        f"**D1-core falsified: {d1['falsified']}**; supported: {d1['supported']}.",
        "",
        "The rule is Amendment 10 §A10.5: at least one block arm must have a "
        "seed-averaged Delta ≥ .02 with a 95% CI lower bound above zero.",
    ]
    if "nested_primary_secondary" in result:
        lines += [
            "",
            "## Nested-subset primary (secondary readout)",
            "",
            "| arm | family-macro F1 | 95% CI | mean micro F1 |",
            "|---|---:|---|---:|",
        ]
        for arm in ARMS:
            row = result["nested_primary_secondary"][arm]
            lines.append(
                f"| {arm} | {row['seed_averaged_family_macro_f1']:.4f} | "
                f"[{row['family_macro_ci95'][0]:.4f}, {row['family_macro_ci95'][1]:.4f}] | "
                f"{row['mean_micro_f1']:.4f} |"
            )
    lines += [
        "",
        "## Reporting invariants",
        "",
        f"- Max macro/micro gap: {result['invariants']['max_macro_micro_gap']:.4f} "
        f"(gate ≤ .15: {result['invariants']['macro_micro_gate_0.15']})",
        "",
        "Round 4 v1-label numbers are not used in this table.",
    ]
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    raise SystemExit(main())
