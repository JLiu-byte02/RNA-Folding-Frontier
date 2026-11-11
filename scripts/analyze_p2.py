#!/usr/bin/env python3
"""Build the P2 table and test the mechanism against its nuisance control.

Two things are computed here, in the order the pre-registration asks for them.

1. **The P2 table.** Family-macro base-pair F1 per distal-context condition, with
   a paired bootstrap over families (the family is the resampling unit because
   sequences inside a family are not independent). "Paired" means one bootstrap
   draw is applied to every condition, so the interval is on the *difference*,
   which is what P2 is about.

2. **Distortion versus clan decodability (Amendment 1.3).** A swap can lower F1
   for two reasons: it removed family identity (the hypothesis) or it pushed the
   input off the pretraining distribution (a nuisance). *Distortion* measures the
   nuisance; the clan probe's probability on the true clan measures the
   hypothesis. F1 is regressed on both per sequence, and the mechanism is only
   credited with the part of the effect that survives conditioning on distortion.

Eval names carry the convention ``<condition>_<split>``, e.g. ``b_test_withheld``.
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np


def read_jsonl(path: Path) -> list[dict]:
    with open(path) as fh:
        return [json.loads(line) for line in fh if line.strip()]


def bootstrap_delta(
    base: dict[str, float], other: dict[str, float], n_boot: int, seed: int
) -> dict:
    """Paired family bootstrap of mean(other) - mean(base)."""
    fams = sorted(set(base) & set(other))
    if len(fams) < 3:
        return {"delta": None, "lo": None, "hi": None, "n_families": len(fams)}
    a = np.array([base[f] for f in fams])
    b = np.array([other[f] for f in fams])
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(fams), size=(n_boot, len(fams)))
    d = b[idx].mean(1) - a[idx].mean(1)
    return {
        "delta": round(float(b.mean() - a.mean()), 4),
        "lo": round(float(np.percentile(d, 2.5)), 4),
        "hi": round(float(np.percentile(d, 97.5)), 4),
        "n_families": len(fams),
        "frac_le_0": round(float((d <= 0).mean()), 4),
    }


def ols(X: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, float]:
    X1 = np.column_stack([np.ones(len(X)), X])
    beta, *_ = np.linalg.lstsq(X1, y, rcond=None)
    resid = y - X1 @ beta
    ss_res = float(resid @ resid)
    ss_tot = float(((y - y.mean()) ** 2).sum())
    return beta, (1.0 - ss_res / ss_tot) if ss_tot else 0.0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", nargs="+", required=True, help="run dir or results.json")
    ap.add_argument("--diagnostics", nargs="*", default=[], help="context_diagnostics *.json")
    ap.add_argument("--per-seq-diag", nargs="*", default=[], help="context_diagnostics *.jsonl")
    ap.add_argument("--n-boot", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=17)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    run_dirs = [
        Path(p).parent if Path(p).name == "results.json" else Path(p) for p in args.runs
    ]
    runs = [json.loads((d / "results.json").read_text()) for d in run_dirs]

    names = sorted(set.intersection(*[set(r["evals"]) for r in runs]))
    table: dict[str, dict] = {}
    per_family: dict[str, dict[str, float]] = {}
    for name in names:
        fams: dict[str, list[float]] = defaultdict(list)
        for run in runs:
            for fam, f1 in run["evals"][name]["per_family"].items():
                fams[fam].append(f1)
        per_family[name] = {f: float(np.mean(v)) for f, v in fams.items()}
        table[name] = {
            "micro_f1_mean": round(
                float(np.mean([r["evals"][name]["micro_f1"] for r in runs])), 4
            ),
            "family_macro_f1_mean": round(float(np.mean(list(per_family[name].values()))), 4),
            "family_macro_f1_per_seed": [r["evals"][name]["family_macro_f1"] for r in runs],
            "n": runs[0]["evals"][name]["n"],
            "by_band": {
                band: round(
                    float(np.mean([r["evals"][name]["by_band"].get(band, np.nan) for r in runs])),
                    4,
                )
                for band in runs[0]["evals"][name]["by_band"]
            },
        }

    for name in names:
        base_name = "a_" + name.split("_", 1)[1]
        if name == base_name or base_name not in per_family:
            continue
        table[name]["vs_native"] = bootstrap_delta(
            per_family[base_name], per_family[name], args.n_boot, args.seed
        )

    report = {"runs": [str(d) for d in run_dirs], "table": table}

    if args.per_seq_diag:
        diag = {}
        for p in args.per_seq_diag:
            for row in read_jsonl(Path(p)):
                diag[(row["id"], row["condition"])] = row
        reg_rows = []
        for d in run_dirs:
            for f in sorted((d / "per_seq").glob("*.jsonl")):
                name = f.stem
                cond = name.split("_", 1)[0]
                for row in read_jsonl(f):
                    key = (row["id"], cond)
                    if key in diag:
                        dr = diag[key]
                        if dr.get("p_true") is None or dr.get("distortion_all") is None:
                            continue
                        reg_rows.append(
                            {
                                "split": name.split("_", 1)[1],
                                "condition": cond,
                                "f1": row["f1"],
                                "distortion": dr["distortion_all"],
                                "p_true": dr["p_true"],
                            }
                        )
        reg: dict[str, dict] = {}
        for key in sorted({(r["split"], r["condition"]) for r in reg_rows}):
            sel = [r for r in reg_rows if (r["split"], r["condition"]) == key]
            if len(sel) < 50:
                continue
            y = np.array([r["f1"] for r in sel])
            dist = np.array([r["distortion"] for r in sel])
            p_true = np.array([r["p_true"] for r in sel])
            _b_dist, r2_dist = ols(dist[:, None], y)
            beta, r2_both = ols(np.column_stack([dist, p_true]), y)
            reg[f"{key[0]}|{key[1]}"] = {
                "n": len(sel),
                "r2_distortion": round(r2_dist, 4),
                "r2_distortion_plus_decodability": round(r2_both, 4),
                "incremental_r2": round(r2_both - r2_dist, 4),
                "beta_intercept_distortion_p_true": [round(float(v), 4) for v in beta],
            }
        report["mechanism_regression"] = reg
        report["regression_rows"] = len(reg_rows)

    if args.diagnostics:
        report["diagnostics"] = [json.loads(Path(p).read_text()) for p in args.diagnostics]

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"table": table}, indent=2))
    if "mechanism_regression" in report:
        print(json.dumps({"mechanism_regression": report["mechanism_regression"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
