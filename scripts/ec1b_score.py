#!/usr/bin/env python3
"""E-C1b: decomposition and no-pretraining one-hot specificity control."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402

from rff.analysis.block_eval import score_dump  # noqa: E402


def rankdata(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(len(values), dtype=float)
    sorted_vals = values[order]
    start = 0
    while start < len(values):
        end = start + 1
        while end < len(values) and sorted_vals[end] == sorted_vals[start]:
            end += 1
        ranks[order[start:end]] = 0.5 * (start + end - 1) + 1
        start = end
    return ranks


def residualize(y: np.ndarray, covariates: np.ndarray) -> np.ndarray:
    design = np.column_stack([np.ones(len(y)), covariates])
    beta, *_ = np.linalg.lstsq(design, y, rcond=None)
    return y - design @ beta


def partial_spearman(x: np.ndarray, y: np.ndarray, cov: np.ndarray) -> float:
    rx, ry = rankdata(x), rankdata(y)
    rc = np.column_stack([rankdata(cov[:, k]) for k in range(cov.shape[1])])
    ex, ey = residualize(rx, rc), residualize(ry, rc)
    if np.std(ex) == 0 or np.std(ey) == 0:
        return float("nan")
    return float(np.corrcoef(ex, ey)[0, 1])


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--ec1", default="experiments/002_e1prime_runs/ec1_p1_p2.json"
    )
    ap.add_argument(
        "--dump-root", default="experiments/002_e1prime_runs"
    )
    ap.add_argument("--seeds", nargs="+", type=int, default=[17, 23, 41])
    ap.add_argument(
        "--out", default="experiments/002_e1prime_runs/ec1b_specificity.json"
    )
    ap.add_argument("--md-out", default="docs/ec1b_2026-11-13.md")
    args = ap.parse_args()

    ec1 = json.loads((ROOT / args.ec1).read_text(encoding="utf-8"))
    families = {row["family"]: row for row in ec1["families"]}
    if not families or "fm_f1" not in next(iter(families.values())):
        raise SystemExit("E-C1 artifact lacks fm_f1/rnafold_f1; rerun ec1_score.py")

    onehot_by_seed = []
    for seed in args.seeds:
        path = ROOT / args.dump_root / f"dump_onehot_seed{seed}" / "test_withheld.jsonl"
        if not path.exists():
            raise SystemExit(f"missing one-hot dump: {path}")
        onehot_by_seed.append(score_dump(path, nested=True))

    common = set(onehot_by_seed[0])
    for table in onehot_by_seed[1:]:
        common &= set(table)
    # Keep only families represented in every one-hot seed.
    fam_ids = {f: [] for f in families}
    for rid in common:
        fam = onehot_by_seed[0][rid].family
        if fam in families:
            fam_ids[fam].append(rid)
    fam_ids = {f: ids for f, ids in fam_ids.items() if ids}

    rows = []
    for family, ids in sorted(fam_ids.items()):
        onehot = float(
            np.mean(
                [
                    np.mean(
                        [
                            next(t[rid].f1 for rid in ids if rid in t)
                            for t in onehot_by_seed
                        ]
                    )
                ]
            )
        )
        base = families[family]
        rows.append(
            {
                "family": family,
                "fm_f1": base["fm_f1"],
                "rnafold_f1": base["rnafold_f1"],
                "onehot_f1": onehot,
                "p1": base["p1"],
                "delta_fm": base["fm_f1"] - base["rnafold_f1"],
                "delta_onehot": onehot - base["rnafold_f1"],
                "cov": [
                    base["median_length"],
                    base["gc"],
                    base["pair_density"],
                    base["multiloop"],
                    base["label_source"],
                ],
            }
        )

    y_fm = np.asarray([r["delta_fm"] for r in rows])
    y_one = np.asarray([r["delta_onehot"] for r in rows])
    x = np.asarray([r["p1"] for r in rows])
    cov = np.asarray([r["cov"] for r in rows])
    rho_fm = partial_spearman(x, y_fm, cov)
    rho_one = partial_spearman(x, y_one, cov)
    rho_rna = partial_spearman(
        x, np.asarray([r["rnafold_f1"] for r in rows]), cov
    )
    delta = rho_fm - rho_one

    rng = np.random.default_rng(20261113)
    vals = []
    for _ in range(10_000):
        idx = rng.integers(0, len(rows), len(rows))
        vals.append(
            partial_spearman(x[idx], y_fm[idx], cov[idx])
            - partial_spearman(x[idx], y_one[idx], cov[idx])
        )
    finite = np.asarray([v for v in vals if np.isfinite(v)])
    lo = float(np.percentile(finite, 2.5)) if len(finite) else float("nan")
    hi = float(np.percentile(finite, 97.5)) if len(finite) else float("nan")
    passes = bool(delta >= 0.15 and lo > 0)

    result = {
        "amendment": "AMENDMENT_13.md",
        "n_families": len(rows),
        "rho_fm": round(rho_fm, 6),
        "rho_rnafold": round(rho_rna, 6),
        "rho_onehot": round(rho_one, 6),
        "delta_rho": round(delta, 6),
        "delta_ci95": [round(lo, 6), round(hi, 6)],
        "passes": passes,
        "families": rows,
    }
    out = ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    md = [
        "# E-C1b specificity control — 2026-11-13",
        "",
        "| quantity | partial ρ |",
        "|---|---:|",
        f"| FM head F1 vs P1 | {rho_fm:.4f} |",
        f"| RNAfold F1 vs P1 | {rho_rna:.4f} |",
        f"| one-hot head ΔF1 vs P1 | {rho_one:.4f} |",
        "",
        f"Δρ = ρ_FM − ρ_onehot = **{delta:.4f}**, 95% CI "
        f"[{lo:.4f}, {hi:.4f}].  ",
        f"Passes (Δρ ≥ .15 and CI excludes 0): **{passes}**.",
        "",
        "The decomposition tells whether the E-C1 signal comes from the FM head, "
        "from RNAfold, or from a generic learned-model advantage.",
    ]
    (ROOT / args.md_out).write_text("\n".join(md) + "\n", encoding="utf-8")
    print(json.dumps(result, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
