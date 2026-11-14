"""Render the common-cohort tables directly from ``paired_<split>.json``.

The 2026-11-14 audits found two numbers in the baseline document that had been
transcribed by hand and drifted from the artifact.  This script makes that
class of error impossible for the method table and the paired-delta table: both
are emitted from the JSON the bootstrap actually wrote, in a fixed row order,
at a fixed precision.

Usage:
    python scripts/render_cohort_tables.py \
        --paired experiments/002_e1prime_runs/cohort/paired_test_withheld.json \
        [--paired ...] --out docs/cohort_tables.generated.md
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

ORDER = ["RNAfold", "EternaFold", "RNAstructure", "MXfold2", "head17", "head23", "head41", "head_mean"]
# The paired JSON stores lower-case tool names for the external baselines.
ALIAS = {name.lower(): name for name in ORDER}
LABEL = {
    "RNAfold": "RNAfold",
    "EternaFold": "EternaFold",
    "RNAstructure": "RNAstructure",
    "MXfold2": "MXfold2 (retrained)",
    "head17": "head, seed 17",
    "head23": "head, seed 23",
    "head41": "head, seed 41",
    "head_mean": "head, 3-seed mean",
}


def row(name: str, m: dict) -> str:
    lo, hi = m["family_macro_f1_ci95"]
    return (
        f"| {LABEL.get(name, name)} | {m['family_macro_f1']:.4f} "
        f"[{lo:.4f}, {hi:.4f}] | {m['micro_f1']:.4f} |"
    )


def delta_row(p: dict) -> str:
    lo, hi = p["ci95"]
    return (
        f"| {p['a']} - {p['b']} | {p['delta']:+.4f} [{lo:+.4f}, {hi:+.4f}] | "
        f"{p['p_two_sided']:.4f} | {p['frac_families_a_better']:.3f} |"
    )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--paired", action="append", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    blocks = []
    for path in args.paired:
        d = json.loads(Path(path).read_text(encoding="utf-8"))
        raw = list(d["methods"])
        names = sorted(raw, key=lambda n: ORDER.index(ALIAS.get(n.lower(), n)) if ALIAS.get(n.lower(), n) in ORDER else len(ORDER))
        names = [ALIAS.get(n.lower(), n) for n in names]
        lookup = {ALIAS.get(n.lower(), n): d["methods"][n] for n in raw}
        blocks.append(f"### {d['split']} — {d['n_sequences']:,} sequences / {d['n_families']} families\n")
        blocks.append("| method | family-macro F1 [95% CI] | micro F1 |")
        blocks.append("|---|---:|---:|")
        blocks += [row(n, lookup[n]) for n in names]
        blocks.append("")
        blocks.append("| comparison | delta [95% CI] | p | frac families better |")
        blocks.append("|---|---:|---:|---:|")
        blocks += [delta_row(p) for p in d["paired"]]
        blocks.append("")
        blocks.append(
            f"Bootstrap: {d['n_boot']:,} family resamples, seed {d['bootstrap_seed']}; "
            f"min true pairs {d['min_true_pairs']}; labels `{d['label_set']}`."
        )
        blocks.append("")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(blocks).rstrip() + "\n", encoding="utf-8", newline="\n")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
