#!/usr/bin/env python3
"""Score the P2-block test-time arm and evaluate H3 (Amendment 10).

Reads `dump_p2block_seed{seed}/{condition}_{split}.jsonl` and writes the
family-macro/micro table, the H3 retention ratios with family-bootstrap CIs,
and the exploratory homolog-identity dose-response.
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
    bootstrap_mean,
    bootstrap_ratio,
    common_ids,
    family_macro,
    family_means,
    micro,
    score_dump,
    seed_average_family_means,
)

CONDITIONS = ("a", "b", "c", "d", "e")
LABEL = {
    "a": "native",
    "b": "dinucleotide shuffle",
    "c": "same-family homolog",
    "d": "other-clan",
    "e": "N mask",
}


def rnd(value: float, nd: int = 6) -> float:
    return round(float(value), nd)


def load_seed_tables(
    dump_root: Path, seed: int, split: str, threshold: float, min_true: int
) -> dict[str, dict]:
    folder = dump_root / f"dump_p2block_seed{seed}"
    tables = {}
    for cond in CONDITIONS:
        path = folder / f"{cond}_{split}.jsonl"
        if path.exists():
            tables[cond] = score_dump(
                path, threshold=threshold, min_true_pairs=min_true
            )
    return tables


def condition_family_values(
    per_seed: dict[int, dict[str, dict]], cond: str, ids: set[str]
) -> dict[str, float]:
    return seed_average_family_means([per_seed[s][cond] for s in per_seed], ids)


def summarize_condition(
    per_seed: dict[int, dict[str, dict]], cond: str
) -> dict:
    ids = common_ids(*[per_seed[s][cond] for s in per_seed])
    means = condition_family_values(per_seed, cond, ids)
    boot = bootstrap_mean(means)
    per_seed_rows = []
    for seed in sorted(per_seed):
        table = per_seed[seed][cond]
        per_seed_rows.append(
            {
                "seed": seed,
                "n": len(ids),
                "family_macro_f1": rnd(family_macro(table, ids)),
                "micro_f1": rnd(micro(table, ids)),
            }
        )
    return {
        "n_sequences": len(ids),
        "n_families": len(means),
        "seed_averaged_family_macro_f1": rnd(boot["mean"]),
        "family_macro_ci95": [rnd(boot["lo"]), rnd(boot["hi"])],
        "per_seed": per_seed_rows,
        "mean_micro_f1": rnd(np.mean([r["micro_f1"] for r in per_seed_rows])),
    }


def identity_for(record_id: str, homolog: dict[str, dict]) -> float | None:
    for key in (record_id, record_id.split(":", 1)[-1]):
        hit = homolog.get(key)
        if hit and hit.get("identity") is not None:
            return float(hit["identity"])
    return None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--dump-root", default="experiments/002_e1prime_runs"
    )
    ap.add_argument("--seeds", nargs="+", type=int, default=[17, 23, 41])
    ap.add_argument("--split", default="test_withheld")
    ap.add_argument("--threshold", type=float, default=0.5)
    ap.add_argument("--min-true-pairs", type=int, default=4)
    ap.add_argument(
        "--homolog-map", default="data/processed/homologs/test_withheld.json"
    )
    ap.add_argument(
        "--out",
        default="experiments/002_e1prime_runs/p2block_h3_test_withheld.json",
    )
    ap.add_argument(
        "--md-out", default="docs/p2block_h3_2026-11-12.md"
    )
    args = ap.parse_args()

    per_seed = {}
    for seed in args.seeds:
        tables = load_seed_tables(
            ROOT / args.dump_root,
            seed,
            args.split,
            args.threshold,
            args.min_true_pairs,
        )
        missing = [c for c in CONDITIONS if c not in tables]
        if missing:
            raise SystemExit(f"seed {seed}: missing conditions {missing}")
        per_seed[seed] = tables

    result: dict = {
        "amendment": "AMENDMENT_10.md",
        "split": args.split,
        "threshold": args.threshold,
        "min_true_pairs": args.min_true_pairs,
        "seeds": args.seeds,
        "conditions": {c: LABEL[c] for c in CONDITIONS},
    }

    # Full-split descriptive table. Condition c is restricted by donor coverage
    # and is reported both on its own population and on the common H3 subset.
    result["full_split"] = {
        cond: summarize_condition(per_seed, cond) for cond in CONDITIONS
    }

    # H3 population: records present for a, c and e in every seed.
    h3_tables = {
        s: common_ids(per_seed[s]["a"], per_seed[s]["c"], per_seed[s]["e"])
        for s in per_seed
    }
    h3_ids = set.intersection(*h3_tables.values())
    result["h3_common_c"] = {
        "n_sequences": len(h3_ids),
        "n_families": len(condition_family_values(per_seed, "a", h3_ids)),
        "per_condition": {
            cond: summarize_condition_on_ids(per_seed, cond, h3_ids)
            for cond in CONDITIONS
        },
    }
    h3_values = {
        cond: condition_family_values(per_seed, cond, h3_ids)
        for cond in CONDITIONS
    }
    rc = bootstrap_ratio(
        h3_values, numerator="c", denominator_a="a", denominator_e="e"
    )
    rd = bootstrap_ratio(
        h3_values, numerator="d", denominator_a="a", denominator_e="e"
    )
    if h3_values["a"] and h3_values["e"]:
        point_den = float(np.mean(list(h3_values["a"].values()))) - float(
            np.mean(list(h3_values["e"].values()))
        )
    else:
        point_den = None
    denominator_ok = bool(point_den is not None and point_den > 0)
    h3_supported = bool(
        denominator_ok
        and rc["undefined_fraction"] <= 0.01
        and rd["undefined_fraction"] <= 0.01
        and rc["lo"] >= 0.90
        and rd["hi"] <= 0.30
    )
    result["h3"] = {
        "estimand": "(condition - e) / (a - e), family-macro all-pairs F1 at 0.5",
        "a_minus_e": rnd(point_den) if point_den is not None else None,
        "R_c": {k: rnd(v) if isinstance(v, float) else v for k, v in rc.items()},
        "R_d": {k: rnd(v) if isinstance(v, float) else v for k, v in rd.items()},
        "supported": h3_supported,
        "criteria": {
            "R_c_lo_ge_0.90": bool(rc["lo"] >= 0.90),
            "R_d_hi_le_0.30": bool(rd["hi"] <= 0.30),
            "denominator_positive": denominator_ok,
            "defined_fraction_ge_0.99": bool(
                rc["undefined_fraction"] <= 0.01
                and rd["undefined_fraction"] <= 0.01
            ),
        },
    }

    # Exploratory dose-response on the same common-c population.
    homolog = json.loads((ROOT / args.homolog_map).read_text())
    identity = {rid: identity_for(rid, homolog) for rid in h3_ids}
    bins = {
        "0.40-0.50": lambda x: x is not None and 0.40 <= x < 0.50,
        "0.50-0.60": lambda x: x is not None and 0.50 <= x < 0.60,
    }
    dose = {}
    for name, keep in bins.items():
        ids = {rid for rid, value in identity.items() if keep(value)}
        if not ids:
            dose[name] = {"n_sequences": 0, "n_families": 0}
            continue
        values = {
            cond: condition_family_values(per_seed, cond, ids)
            for cond in ("a", "c", "e")
        }
        ratio = bootstrap_ratio(
            values, numerator="c", denominator_a="a", denominator_e="e"
        )
        dose[name] = {
            "n_sequences": len(ids),
            "n_families": ratio["n_families"],
            "R_c": {k: rnd(v) if isinstance(v, float) else v for k, v in ratio.items()},
        }
    result["dose_response"] = dose

    # Reporting invariants (Amendment 9 §A9.6).
    gaps = []
    for cond in CONDITIONS:
        for row in result["full_split"][cond]["per_seed"]:
            gaps.append(
                {
                    "condition": cond,
                    "seed": row["seed"],
                    "gap": abs(row["family_macro_f1"] - row["micro_f1"]),
                }
            )
    result["invariants"] = {
        "max_macro_micro_gap": rnd(max(g["gap"] for g in gaps)),
        "max_macro_micro_gap_row": max(gaps, key=lambda g: g["gap"]),
        "macro_micro_gate_0.15": bool(max(g["gap"] for g in gaps) <= 0.15),
        "shuffle_floor_max": rnd(
            max(
                row["family_macro_f1"]
                for row in result["full_split"]["b"]["per_seed"]
            )
        ),
        "shuffle_floor_gate_0.05": bool(
            max(
                row["family_macro_f1"]
                for row in result["full_split"]["b"]["per_seed"]
            )
            <= 0.05
        ),
    }

    out = ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    md = render_markdown(result)
    md_path = ROOT / args.md_out
    md_path.parent.mkdir(parents=True, exist_ok=True)
    md_path.write_text(md, encoding="utf-8")
    print(
        json.dumps(
            {
                "out": args.out,
                "md": args.md_out,
                "h3_supported": h3_supported,
                "R_c": result["h3"]["R_c"],
                "R_d": result["h3"]["R_d"],
            }
        )
    )
    return 0


def summarize_condition_on_ids(
    per_seed: dict[int, dict[str, dict]], cond: str, ids: set[str]
) -> dict:
    means = condition_family_values(per_seed, cond, ids)
    boot = bootstrap_mean(means)
    rows = []
    for seed in sorted(per_seed):
        table = per_seed[seed][cond]
        rows.append(
            {
                "seed": seed,
                "family_macro_f1": rnd(family_macro(table, ids)),
                "micro_f1": rnd(micro(table, ids)),
            }
        )
    return {
        "n_families": len(means),
        "seed_averaged_family_macro_f1": rnd(boot["mean"]),
        "family_macro_ci95": [rnd(boot["lo"]), rnd(boot["hi"])],
        "per_seed": rows,
        "mean_micro_f1": rnd(np.mean([r["micro_f1"] for r in rows])),
    }


def render_markdown(result: dict) -> str:
    lines = [
        "# P2-block test-time arm and H3",
        "",
        "Generated by `scripts/score_p2_block.py` from the dumps; no GPU is used.",
        "",
        f"Split: `{result['split']}`; fixed threshold {result['threshold']}; "
        f"v3 all-pairs truth; sequences with fewer than {result['min_true_pairs']} "
        "true pairs excluded.",
        "",
        "## Full-split descriptive table",
        "",
        "| condition | n | families | family-macro F1 | 95% CI | mean micro F1 |",
        "|---|---:|---:|---:|---|---:|",
    ]
    for cond in CONDITIONS:
        row = result["full_split"][cond]
        lines.append(
            f"| {cond} ({LABEL[cond]}) | {row['n_sequences']} | {row['n_families']} | "
            f"{row['seed_averaged_family_macro_f1']:.4f} | "
            f"[{row['family_macro_ci95'][0]:.4f}, {row['family_macro_ci95'][1]:.4f}] | "
            f"{row['mean_micro_f1']:.4f} |"
        )
    h3 = result["h3"]
    lines += [
        "",
        "Condition c is restricted to records with a verified homolog; the H3 "
        "comparison below uses the common a/c/e population.",
        "",
        "## H3: evolutionarily consistent context transfers",
        "",
        f"Common a/c/e population: {result['h3_common_c']['n_sequences']} sequences, "
        f"{result['h3_common_c']['n_families']} families.",
        "",
        "| quantity | point | 95% CI | undefined draws |",
        "|---|---:|---|---:|",
        f"| R_c | {h3['R_c']['ratio']:.4f} | "
        f"[{h3['R_c']['lo']:.4f}, {h3['R_c']['hi']:.4f}] | "
        f"{h3['R_c']['undefined_fraction']:.4f} |",
        f"| R_d | {h3['R_d']['ratio']:.4f} | "
        f"[{h3['R_d']['lo']:.4f}, {h3['R_d']['hi']:.4f}] | "
        f"{h3['R_d']['undefined_fraction']:.4f} |",
        "",
        f"**H3 supported: {h3['supported']}** "
        f"(R_c lower bound ≥ .90, R_d upper bound ≤ .30, positive denominator).",
        "",
        "## Exploratory dose-response",
        "",
        "| homolog identity | n | families | R_c | 95% CI |",
        "|---|---:|---:|---:|---|",
    ]
    for name, row in result["dose_response"].items():
        if row.get("n_families"):
            rc = row["R_c"]
            lines.append(
                f"| {name} | {row['n_sequences']} | {row['n_families']} | "
                f"{rc['ratio']:.4f} | [{rc['lo']:.4f}, {rc['hi']:.4f}] |"
            )
        else:
            lines.append(f"| {name} | 0 | 0 | — | — |")
    inv = result["invariants"]
    lines += [
        "",
        "## Reporting invariants",
        "",
        f"- Max macro/micro gap: {inv['max_macro_micro_gap']:.4f} "
        f"(gate ≤ .15: {inv['macro_micro_gate_0.15']})",
        f"- Shuffle-floor max family-macro F1: {inv['shuffle_floor_max']:.4f} "
        f"(gate ≤ .05: {inv['shuffle_floor_gate_0.05']})",
        "",
        "Round 4 v1-label numbers are not used in this table.",
    ]
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    raise SystemExit(main())
