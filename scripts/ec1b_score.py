#!/usr/bin/env python3
"""E-C1b: decomposition and no-pretraining one-hot specificity control.

Both sides of the contrast are recomputed on ONE common per-family id set:
ids that are present in the FM dump (with >= min_true_pairs nested v3 truth)
and in every one-hot seed dump.  The earlier version mixed an E-C1 family mean
(FM side, its own id set) with a one-hot family mean over a different id set,
which can bias the paired delta.  RNAfold F1 is recomputed on the same ids and
cached, because the E-C1 artifact stores only family-level deltas.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402

from rff.analysis.block_eval import score_dump  # noqa: E402
from rff.analysis.labels import max_nested_subset  # noqa: E402
from rff.data.io import pairs_from_dbn  # noqa: E402


def read_jsonl(path: Path) -> list[dict]:
    with open(path, encoding="utf-8") as fh:
        return [json.loads(line) for line in fh if line.strip()]


def f1_from_counts(tp: int, fp: int, fn: int) -> float:
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


def partial_spearman(x: np.ndarray, y: np.ndarray, cov: np.ndarray) -> float:
    rx, ry = rankdata(x), rankdata(y)
    rc = np.column_stack([rankdata(cov[:, k]) for k in range(cov.shape[1])])
    ex, ey = residualize(rx, rc), residualize(ry, rc)
    if np.std(ex) == 0 or np.std(ey) == 0:
        return float("nan")
    return float(np.corrcoef(ex, ey)[0, 1])


def load_rnafold_cache(path: Path) -> dict[str, dict]:
    cache: dict[str, dict] = {}
    if not path.exists():
        return cache
    for row in read_jsonl(path):
        cache[row["id"]] = row
    return cache


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ec1", default="experiments/002_e1prime_runs/ec1_p1_p2.json")
    ap.add_argument(
        "--families", default="experiments/002_e1prime_runs/ec1_families.json"
    )
    ap.add_argument("--split", default="splits/s1prime_filtered/test_withheld.jsonl")
    ap.add_argument(
        "--fm-dump",
        default="experiments/002_e1prime_runs/dump_p2fix_seed17/a_test_withheld.jsonl",
    )
    ap.add_argument("--dump-root", default="experiments/002_e1prime_runs")
    ap.add_argument("--seeds", nargs="+", type=int, default=[17, 23, 41])
    ap.add_argument("--min-true-pairs", type=int, default=4)
    ap.add_argument(
        "--rnafold-cache",
        default="experiments/002_e1prime_runs/ec1b_rnafold_counts.jsonl",
    )
    ap.add_argument(
        "--out", default="experiments/002_e1prime_runs/ec1b_specificity.json"
    )
    ap.add_argument("--md-out", default="docs/ec1b_2026-11-13.md")
    ap.add_argument("--n-boot", type=int, default=10_000)
    args = ap.parse_args()

    from RNA import fold

    ec1 = json.loads((ROOT / args.ec1).read_text(encoding="utf-8"))
    stored = {row["family"]: row for row in ec1["families"]}
    fam_data = json.loads((ROOT / args.families).read_text(encoding="utf-8"))
    wanted: dict[str, set[str]] = {
        fam: set(entry["ids"]) for fam, entry in fam_data["families"].items()
    }
    sources = {
        fam: set(entry.get("sources", []))
        for fam, entry in fam_data["families"].items()
    }

    split_rows = {
        row["id"]: row for row in read_jsonl(ROOT / args.split)
    }

    # FM side: nested v3 truth with at least min_true_pairs pairs.
    fm_table = score_dump(
        ROOT / args.fm_dump,
        nested=True,
        min_true_pairs=args.min_true_pairs,
    )
    truth_by_id: dict[str, set[tuple[int, int]]] = {}
    length_by_id: dict[str, int] = {}
    for row in read_jsonl(ROOT / args.fm_dump):
        rid = row["id"]
        if rid not in fm_table:
            continue
        truth_by_id[rid] = max_nested_subset(
            {tuple(int(x) for x in pair) for pair in row["truth_v3"]}
        )
        length_by_id[rid] = int(row["length"])

    onehot_tables = []
    for seed in args.seeds:
        path = ROOT / args.dump_root / f"dump_onehot_seed{seed}" / "test_withheld.jsonl"
        if not path.exists():
            raise SystemExit(f"missing one-hot dump: {path}")
        onehot_tables.append(score_dump(path, nested=True))

    common = set(fm_table)
    for table in onehot_tables:
        common &= set(table)

    # Restrict to the E-C1 family universe and to the ids that family listed.
    fam_ids: dict[str, list[str]] = {}
    skipped_seq_mismatch = 0
    for rid in sorted(common):
        fam = fm_table[rid].family
        if fam not in wanted or rid not in wanted[fam]:
            continue
        if rid not in split_rows:
            continue
        seq = split_rows[rid]["seq"]
        if len(seq) != length_by_id[rid]:
            skipped_seq_mismatch += 1
            continue
        fam_ids.setdefault(fam, []).append(rid)

    cache_path = ROOT / args.rnafold_cache
    cache = load_rnafold_cache(cache_path)
    new_cache_rows = []
    for fam, ids in sorted(fam_ids.items()):
        for rid in ids:
            if rid in cache:
                continue
            dbn, _ = fold(split_rows[rid]["seq"])
            pred = set(pairs_from_dbn(dbn))
            truth = truth_by_id[rid]
            row = {
                "id": rid,
                "tp": len(truth & pred),
                "fp": len(pred - truth),
                "fn": len(truth - pred),
                "length": length_by_id[rid],
            }
            cache[rid] = row
            new_cache_rows.append(row)
    if new_cache_rows:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        with open(cache_path, "w", encoding="utf-8") as fh:
            for rid in sorted(cache):
                fh.write(json.dumps(cache[rid], sort_keys=True) + "\n")

    rows = []
    for fam, ids in sorted(fam_ids.items()):
        fm_f1 = float(np.mean([fm_table[rid].f1 for rid in ids]))
        rna_f1 = float(
            np.mean(
                [
                    f1_from_counts(
                        cache[rid]["tp"], cache[rid]["fp"], cache[rid]["fn"]
                    )
                    for rid in ids
                ]
            )
        )
        onehot_f1 = float(
            np.mean(
                [
                    np.mean([table[rid].f1 for rid in ids])
                    for table in onehot_tables
                ]
            )
        )
        lengths = [length_by_id[rid] for rid in ids]
        gcs = [
            (split_rows[rid]["seq"].count("G") + split_rows[rid]["seq"].count("C"))
            / len(split_rows[rid]["seq"])
            for rid in ids
        ]
        densities = [len(truth_by_id[rid]) / length_by_id[rid] for rid in ids]
        multi = float(
            any(has_multiloop(truth_by_id[rid]) for rid in ids)
        )
        base = stored.get(fam, {})
        rows.append(
            {
                "family": fam,
                "n_common": len(ids),
                "n_ec1": int(base.get("n", 0)),
                "fm_f1": fm_f1,
                "rnafold_f1": rna_f1,
                "onehot_f1": onehot_f1,
                "p1": float(base.get("p1", float("nan"))),
                "delta_fm": fm_f1 - rna_f1,
                "delta_onehot": onehot_f1 - rna_f1,
                "delta_f1_ec1": float(base.get("delta_f1", float("nan"))),
                "cov": [
                    float(np.median(lengths)),
                    float(np.mean(gcs)),
                    float(np.mean(densities)),
                    multi,
                    1.0 if "rfam" in sources.get(fam, set()) else 0.0,
                ],
            }
        )

    rows = [r for r in rows if np.isfinite(r["p1"]) and r["n_common"] > 0]
    x = np.asarray([r["p1"] for r in rows])
    y_fm = np.asarray([r["delta_fm"] for r in rows])
    y_one = np.asarray([r["delta_onehot"] for r in rows])
    cov = np.asarray([r["cov"] for r in rows])

    rho_fm = partial_spearman(x, y_fm, cov)
    rho_one = partial_spearman(x, y_one, cov)
    rho_rna = partial_spearman(x, np.asarray([r["rnafold_f1"] for r in rows]), cov)
    delta = rho_fm - rho_one

    rng = np.random.default_rng(20261113)
    vals = []
    for _ in range(args.n_boot):
        idx = rng.integers(0, len(rows), len(rows))
        vals.append(
            partial_spearman(x[idx], y_fm[idx], cov[idx])
            - partial_spearman(x[idx], y_one[idx], cov[idx])
        )
    finite = np.asarray([v for v in vals if np.isfinite(v)])
    lo = float(np.percentile(finite, 2.5)) if len(finite) else float("nan")
    hi = float(np.percentile(finite, 97.5)) if len(finite) else float("nan")
    passes = bool(delta >= 0.15 and lo > 0)

    ec1_delta = np.asarray([r["delta_f1_ec1"] for r in rows if np.isfinite(r["delta_f1_ec1"])])
    result = {
        "amendment": "AMENDMENT_13.md",
        "population": "common ids: FM dump with >= %d nested v3 pairs AND all %d one-hot seed dumps"
        % (args.min_true_pairs, len(args.seeds)),
        "n_families": len(rows),
        "n_sequences": int(sum(r["n_common"] for r in rows)),
        "n_sequences_ec1": int(sum(r["n_ec1"] for r in rows)),
        "skipped_seq_length_mismatch": skipped_seq_mismatch,
        "seeds": args.seeds,
        "rho_fm": round(rho_fm, 6),
        "rho_rnafold": round(rho_rna, 6),
        "rho_onehot": round(rho_one, 6),
        "delta_rho": round(delta, 6),
        "delta_ci95": [round(lo, 6), round(hi, 6)],
        "passes": passes,
        "mean_delta_fm": round(float(y_fm.mean()), 6),
        "mean_delta_onehot": round(float(y_one.mean()), 6),
        "mean_delta_f1_ec1_same_families": round(float(ec1_delta.mean()), 6)
        if ec1_delta.size
        else None,
        "families": rows,
    }
    out = ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    md = [
        "# E-C1b specificity control — 2026-11-13",
        "",
        f"Common population: {result['n_sequences']} sequences in "
        f"{result['n_families']} families (ids present in the FM dump with "
        f">= {args.min_true_pairs} nested v3 pairs and in all "
        f"{len(args.seeds)} one-hot seed dumps).",
        "",
        "| quantity | partial rho vs P1 |",
        "|---|---:|",
        f"| FM head F1 (level) | {rho_fm:.4f} |",
        f"| RNAfold F1 (level) | {rho_rna:.4f} |",
        f"| FM head dF1 vs RNAfold | {rho_fm:.4f} |",
        f"| one-hot head dF1 vs RNAfold | {rho_one:.4f} |",
        "",
        f"d_rho = rho_FM - rho_onehot = **{delta:.4f}**, 95% family-bootstrap CI "
        f"[{lo:.4f}, {hi:.4f}].  ",
        f"Passes (d_rho >= .15 and CI excludes 0): **{passes}**.",
        "",
        f"Mean dF1: FM {y_fm.mean():+.4f}, one-hot {y_one.mean():+.4f}. "
        f"For reference, the E-C1 (unfiltered) mean dF1 on these same families "
        f"is {result['mean_delta_f1_ec1_same_families']:+.4f}."
        if result["mean_delta_f1_ec1_same_families"] is not None
        else "",
        "",
        "Both arms are scored on identical ids, identical nested v3 truth, and "
        "identical decoder settings, so the paired difference is population-"
        "consistent.",
    ]
    (ROOT / args.md_out).write_text("\n".join(md) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                k: result[k]
                for k in (
                    "n_families",
                    "n_sequences",
                    "n_sequences_ec1",
                    "skipped_seq_length_mismatch",
                    "rho_fm",
                    "rho_rnafold",
                    "rho_onehot",
                    "delta_rho",
                    "delta_ci95",
                    "passes",
                )
            }
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
