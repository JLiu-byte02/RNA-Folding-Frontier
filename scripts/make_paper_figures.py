"""Deterministic paper figures rendered from committed artifacts.

Every number in these figures is read from a JSON artifact; nothing is typed in
by hand, so a figure cannot drift from the table. Rendering uses Agg, a fixed
palette, and fixed sizes, so re-running produces the same files.

    python scripts/make_paper_figures.py --out-dir figures
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
COHORT = ROOT / "experiments" / "002_e1prime_runs" / "cohort"

PALETTE = {
    "rnafold": "#1f77b4",
    "eternafold": "#2ca02c",
    "rnastructure": "#9467bd",
    "mxfold2": "#7f7f7f",
    "head_mean": "#d62728",
}
ORDER = ["rnafold", "eternafold", "rnastructure", "mxfold2", "head_mean"]
LABELS = {
    "rnafold": "RNAfold",
    "eternafold": "EternaFold",
    "rnastructure": "RNAstructure",
    "mxfold2": "MXfold2 (retrained)",
    "head_mean": "head (3-seed mean)",
}


def load_paired(split: str) -> dict:
    return json.loads((COHORT / f"paired_{split}.json").read_text(encoding="utf-8"))


def method(entry: dict, key: str) -> dict:
    for name, value in entry["methods"].items():
        if name.lower() == key:
            return value
    raise KeyError(key)


def fig7(out_dir: Path) -> None:
    splits = [("test_withheld", "Withheld"), ("id_test", "id_test")]
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.4), sharey=True)
    for ax, (split, title) in zip(axes, splits):
        d = load_paired(split)
        ys = list(range(len(ORDER)))[::-1]
        for y, key in zip(ys, ORDER):
            m = method(d, key)
            lo, hi = m["family_macro_f1_ci95"]
            ax.plot([lo, hi], [y, y], color=PALETTE[key], lw=2, solid_capstyle="butt")
            ax.plot(m["family_macro_f1"], y, "o", color=PALETTE[key], ms=7)
        ax.axvline(method(d, "eternafold")["family_macro_f1"], color="#bbbbbb", ls=":", lw=1, zorder=0)
        ax.set_yticks(ys)
        ax.set_yticklabels([LABELS[k] for k in ORDER])
        ax.set_xlabel("family-macro F1 (95% CI)")
        ax.set_title(f"{title} — {d['n_sequences']:,} seqs / {d['n_families']} families", fontsize=10)
        ax.grid(axis="x", color="#eeeeee", lw=0.8)
        ax.set_axisbelow(True)
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)
    fig.suptitle("Frozen common cohort: no head advantage over the thermodynamic tools", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    save(fig, out_dir / "fig7_four_method_cohort")


def fig2(out_dir: Path) -> None:
    v3 = json.loads((ROOT / "experiments/011_labels_v3/labels_v3_summary.json").read_text(encoding="utf-8"))
    n = v3["n_sequences"]
    parse = [
        1 - v3["n_v1_parse_fail"] / n,
        1 - v3["n_v2_parse_fail"] / n,
        1 - v3["n_v3_parse_fail"] / n,
    ]
    pairs = [v3["pairs_effective_v1"], v3["pairs_effective_v2"], v3["pairs_v3"]]

    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.2))
    names = ["v1 string", "v2 char-gap", "v3 pair-level"]
    colors = ["#999999", "#4c72b0", "#c44e52"]
    # Dots rather than bars: the axis is truncated at 0.75 to show the
    # differences, and a truncated bar chart would misrepresent their size.
    for x, (value, colour) in enumerate(zip(parse, colors)):
        axes[0].plot(x, value, "o", ms=10, color=colour)
        axes[0].text(x, value + 0.008, f"{value:.3f}", ha="center", fontsize=9)
    axes[0].set_xticks(range(len(names)))
    axes[0].set_xticklabels(names)
    axes[0].set_ylim(0.75, 1.02)
    axes[0].set_ylabel("fraction of rows that parse")
    axes[0].set_title(f"Parse rate ({n:,} seed rows; axis starts at 0.75)", fontsize=10)

    bars2 = axes[1].bar(names, [x / 1e6 for x in pairs], color=colors)
    axes[1].set_ylabel("total pairs (millions)")
    for b, p in zip(bars2, pairs):
        axes[1].text(b.get_x() + b.get_width() / 2, p / 1e6 + 0.03, f"{p/1e6:.2f}M", ha="center", fontsize=9)
    axes[1].set_title("Recovered base pairs", fontsize=10)

    for ax in axes:
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)
        ax.grid(axis="y", color="#eeeeee", lw=0.8)
        ax.set_axisbelow(True)
    fig.suptitle("String-level SS_cons projection loses labels; pair-level projection does not", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    save(fig, out_dir / "fig2_label_projection")


def save(fig, stem: Path) -> None:
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)
    print("wrote", stem.with_suffix(".png"), "and .svg")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default="figures")
    args = ap.parse_args()
    out = Path(args.out_dir)
    if not out.is_absolute():
        out = ROOT / out
    fig7(out)
    fig2(out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
