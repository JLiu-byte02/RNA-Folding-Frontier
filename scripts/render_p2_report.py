#!/usr/bin/env python3
"""Render the P2 evidence as a committed markdown table.

Pulls the numbers from the run directories and the diagnostics files so the
report cannot drift from the artifacts: every cell is read from a file that is
in the repo or referenced by hash.
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]

COND_LABEL = {
    "a": "a native",
    "b": "b dinuc-shuffled",
    "c": "c homolog <60% id",
    "d": "d other clan",
    "e": "e N-masked",
}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", nargs="+", required=True)
    ap.add_argument("--diagnostics", nargs="*", default=[])
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    run_dirs = [Path(p) for p in args.runs]
    runs = [json.loads((d / "results.json").read_text()) for d in run_dirs]
    names = sorted(set.intersection(*[set(r["evals"]) for r in runs]))

    rows: dict[str, dict[str, float]] = defaultdict(dict)
    for name in names:
        rows[name]["micro"] = float(np.mean([r["evals"][name]["micro_f1"] for r in runs]))
        fams = defaultdict(list)
        for r in runs:
            for fam, f1 in r["evals"][name]["per_family"].items():
                fams[fam].append(f1)
        rows[name]["macro"] = float(np.mean([np.mean(v) for v in fams.values()]))
        rows[name]["n"] = runs[0]["evals"][name]["n"]
        rows[name]["bands"] = {
            b: float(
                np.mean(
                    [r["evals"][name]["by_band"].get(b, np.nan) for r in runs]
                )
            )
            for b in runs[0]["evals"][name]["by_band"]
        }

    lines = ["# P2 — test-time distal-context intervention", ""]
    lines.append(
        "One set of decoder-B weights per seed, trained on native inputs only "
        "(condition a). At evaluation the distal region is rewritten and nothing "
        "else changes. Numbers are means over seeds "
        f"({', '.join(d.name for d in run_dirs)})."
    )
    lines.append("")

    for split in sorted({n.split("_", 1)[1] for n in names}):
        lines.append(f"## {split}")
        lines.append("")
        lines.append("| condition | n | micro F1 | family-macro F1 |")
        lines.append("|---|---|---|---|")
        for cond in "abcde":
            name = f"{cond}_{split}"
            if name not in rows:
                continue
            lines.append(
                f"| {COND_LABEL[cond]} | {rows[name]['n']} | "
                f"{rows[name]['micro']:.4f} | {rows[name]['macro']:.4f} |"
            )
        lines.append("")
        bands = sorted({b for n in names if n.endswith(split) for b in rows[n]["bands"]})
        if bands:
            header = "| condition | " + " | ".join(bands) + " |"
            lines.append(header)
            lines.append("|" + "---|" * (len(bands) + 1))
            for cond in "abcde":
                name = f"{cond}_{split}"
                if name not in rows:
                    continue
                cells = " | ".join(f"{rows[name]['bands'].get(b, float('nan')):.3f}" for b in bands)
                lines.append(f"| {cond} | {cells} |")
            lines.append("")

    if args.diagnostics:
        lines.append("## What the swap did to the representation")
        lines.append("")
        lines.append(
            "`distal` restricts the cosine distance to the rewritten positions; "
            "`local` measures positions whose tokens were untouched, i.e. the "
            "indirect effect of global attention."
        )
        lines.append("")
        lines.append("| split | condition | distortion (all) | distal | local | probe bal-acc |")
        lines.append("|---|---|---|---|---|---|")
        for p in args.diagnostics:
            d = json.loads(Path(p).read_text())
            split = Path(d["split"]).stem
            for cond in "abcde":
                c = d["conditions"].get(cond)
                if not c:
                    continue
                pr = c.get("probe") or {}
                ba = pr.get("balanced_acc")
                lines.append(
                    f"| {split} | {cond} | "
                    f"{c.get('distortion_all') if c.get('distortion_all') is not None else '—'} | "
                    f"{c.get('distortion_distal') if c.get('distortion_distal') is not None else '—'} | "
                    f"{c.get('distortion_local') if c.get('distortion_local') is not None else '—'} | "
                    f"{ba if ba is not None else 'n/a'} |"
                )
        lines.append("")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
