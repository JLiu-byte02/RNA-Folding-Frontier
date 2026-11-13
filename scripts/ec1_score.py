#!/usr/bin/env python3
"""E-C1: coverage proxies versus family-level FM advantage (Amendment 12)."""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402

from rff.analysis.labels import max_nested_subset  # noqa: E402
from rff.analysis.scoring import decode_pairs, normalise_scored  # noqa: E402
from rff.data.io import pairs_from_dbn  # noqa: E402


def read_jsonl(path: Path) -> list[dict]:
    with open(path, encoding="utf-8") as fh:
        return [json.loads(line) for line in fh if line.strip()]


def f1(truth: set[tuple[int, int]], pred: set[tuple[int, int]]) -> float:
    tp = len(truth & pred)
    fp = len(pred - truth)
    fn = len(truth - pred)
    if tp + fp + fn == 0:
        return 1.0
    if tp == 0:
        return 0.0
    return 2 * tp / (2 * tp + fp + fn)


def has_multiloop(pairs: set[tuple[int, int]]) -> bool:
    ordered = sorted(pairs, key=lambda p: (p[0], -p[1]))
    children: dict[tuple[int, int], int] = {}
    stack: list[tuple[int, int]] = []
    for i, j in ordered:
        while stack and not (stack[-1][0] < i < j < stack[-1][1]):
            stack.pop()
        if stack:
            children[stack[-1]] = children.get(stack[-1], 0) + 1
        stack.append((i, j))
    return any(count > 1 for count in children.values())


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


def partial_spearman(
    x: np.ndarray, y: np.ndarray, covariates: np.ndarray
) -> float:
    rx = rankdata(x)
    ry = rankdata(y)
    cov_ranks = np.column_stack([rankdata(covariates[:, k]) for k in range(covariates.shape[1])])
    ex = residualize(rx, cov_ranks)
    ey = residualize(ry, cov_ranks)
    if np.std(ex) == 0 or np.std(ey) == 0:
        return float("nan")
    return float(np.corrcoef(ex, ey)[0, 1])


def bootstrap_partial(
    x: np.ndarray,
    y: np.ndarray,
    covariates: np.ndarray,
    *,
    n_boot: int = 10_000,
    seed: int = 20261113,
) -> dict:
    point = partial_spearman(x, y, covariates)
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(n_boot):
        idx = rng.integers(0, len(x), size=len(x))
        vals.append(partial_spearman(x[idx], y[idx], covariates[idx]))
    finite = np.asarray([v for v in vals if np.isfinite(v)])
    if not np.isfinite(point) or not len(finite):
        return {"rho": float("nan"), "lo": float("nan"), "hi": float("nan")}
    return {
        "rho": round(point, 6),
        "lo": round(float(np.percentile(finite, 2.5)), 6),
        "hi": round(float(np.percentile(finite, 97.5)), 6),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--families", default="experiments/002_e1prime_runs/ec1_families.json"
    )
    ap.add_argument("--split", default="splits/s1prime_filtered/test_withheld.jsonl")
    ap.add_argument(
        "--fm-dump",
        default="experiments/002_e1prime_runs/dump_p2fix_seed17/a_test_withheld.jsonl",
    )
    ap.add_argument("--rfam-family", default="data/raw/rfam/family.txt")
    ap.add_argument(
        "--p2", default="experiments/002_e1prime_runs/ec1_p2_pppl.json"
    )
    ap.add_argument("--threshold", type=float, default=0.5)
    ap.add_argument(
        "--out", default="experiments/002_e1prime_runs/ec1_p1_p2.json"
    )
    ap.add_argument("--md-out", default="docs/ec1_2026-11-13.md")
    args = ap.parse_args()

    from RNA import fold

    family_data = json.loads((ROOT / args.families).read_text(encoding="utf-8"))
    wanted = {
        rid for entry in family_data["families"].values() for rid in entry["ids"]
    }
    seqs = {
        rec["id"]: rec["seq"]
        for rec in read_jsonl(ROOT / args.split)
        if rec["id"] in wanted
    }
    fm_rows = {
        row["id"]: row for row in read_jsonl(ROOT / args.fm_dump) if row["id"] in wanted
    }

    full_member: dict[str, int] = {}
    with open(ROOT / args.rfam_family, encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 16 and parts[15].isdigit():
                full_member[parts[0]] = int(parts[15])

    p2_values = {}
    if (ROOT / args.p2).exists():
        blob = json.loads((ROOT / args.p2).read_text(encoding="utf-8"))
        for rid, row in blob["values"].items():
            p2_values.setdefault(row["family"], []).append(row["pppl"])

    per_family = []
    for family, entry in sorted(family_data["families"].items()):
        ids = entry["ids"]
        if family not in full_member:
            raise SystemExit(f"family missing from family.txt: {family}")
        fm_f1 = []
        rnafold_f1 = []
        lengths = []
        gc = []
        density = []
        multiloop = False
        for rid in ids:
            seq = seqs[rid]
            row = fm_rows[rid]
            truth_raw = row.get("truth_v3")
            if truth_raw is None:
                continue
            truth = max_nested_subset({tuple(p) for p in truth_raw})
            pred_fm = set(
                decode_pairs(
                    normalise_scored(row["scored"]),
                    args.threshold,
                    len(seq),
                    partner=True,
                    nesting=True,
                )
            )
            dbn, _ = fold(seq)
            pred_rna = set(pairs_from_dbn(dbn))
            fm_f1.append(f1(truth, pred_fm))
            rnafold_f1.append(f1(truth, pred_rna))
            lengths.append(len(seq))
            gc.append((seq.count("G") + seq.count("C")) / len(seq))
            density.append(len(truth) / len(seq))
            multiloop = multiloop or has_multiloop(truth)
        if not fm_f1:
            continue
        entry_out = {
            "family": family,
            "n": len(fm_f1),
            "fm_f1": float(np.mean(fm_f1)),
            "rnafold_f1": float(np.mean(rnafold_f1)),
            "delta_f1": float(np.mean(fm_f1) - np.mean(rnafold_f1)),
            "p1": math.log1p(full_member[family]),
            "p2": float(np.median(p2_values[family])) if family in p2_values else None,
            "median_length": float(np.median(lengths)),
            "gc": float(np.mean(gc)),
            "pair_density": float(np.mean(density)),
            "multiloop": float(multiloop),
            "label_source": 1.0 if "rfam" in entry.get("sources", []) else 0.0,
        }
        per_family.append(entry_out)

    x = np.asarray([r["p1"] for r in per_family])
    y = np.asarray([r["delta_f1"] for r in per_family])
    cov = np.asarray(
        [
            [
                r["median_length"],
                r["gc"],
                r["pair_density"],
                r["multiloop"],
                r["label_source"],
            ]
            for r in per_family
        ]
    )
    p1 = bootstrap_partial(x, y, cov)
    p2_rows = [r for r in per_family if r["p2"] is not None]
    p2 = (
        bootstrap_partial(
            np.asarray([r["p2"] for r in p2_rows]),
            np.asarray([r["delta_f1"] for r in p2_rows]),
            np.asarray(
                [
                    [
                        r["median_length"],
                        r["gc"],
                        r["pair_density"],
                        r["multiloop"],
                        r["label_source"],
                    ]
                    for r in p2_rows
                ]
            ),
        )
        if p2_rows
        else None
    )

    supported = bool(
        p1["rho"] >= 0.30 and p1["lo"] > 0 and p2 is not None and p2["rho"] <= -0.20
    )
    falsified = bool(p1["rho"] < 0.10 or p1["lo"] <= 0 <= p1["hi"])
    result = {
        "amendment": "AMENDMENT_12.md",
        "n_families": len(per_family),
        "p1": p1,
        "p2": p2,
        "supported": supported,
        "falsified": falsified,
        "inconclusive": bool(not supported and not falsified),
        "families": per_family,
    }
    out = ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    md_path = ROOT / args.md_out
    md_path.parent.mkdir(parents=True, exist_ok=True)
    md_path.write_text(render_markdown(result), encoding="utf-8")
    print(json.dumps({k: result[k] for k in ("n_families", "p1", "p2", "supported", "falsified", "inconclusive")}))
    return 0


def render_markdown(result: dict) -> str:
    lines = [
        "# E-C1 coverage–advantage test — 2026-11-13",
        "",
        "Amendment 12; 409-family common withheld population; outcome is the",
        "family-level FM advantage over RNAfold on nested v3 truth.",
        "",
        "| proxy | partial rho | 95% CI |",
        "|---|---:|---|",
        f"| P1 Rfam full-region membership | {result['p1']['rho']:.4f} | "
        f"[{result['p1']['lo']:.4f}, {result['p1']['hi']:.4f}] |",
    ]
    if result["p2"]:
        lines.append(
            f"| P2 RiNALMo pseudo-perplexity | {result['p2']['rho']:.4f} | "
            f"[{result['p2']['lo']:.4f}, {result['p2']['hi']:.4f}] |"
        )
    lines += [
        "",
        f"Supported: **{result['supported']}**; falsified: "
        f"**{result['falsified']}**; inconclusive: **{result['inconclusive']}**.",
        "",
        "The pre-registered success rule was P1 rho ≥ .30 with a CI excluding 0,",
        "and P2 in the negative direction (rho ≤ −.20).",
    ]
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    raise SystemExit(main())
