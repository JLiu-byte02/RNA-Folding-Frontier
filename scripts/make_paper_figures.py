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



def fig1(out_dir: Path) -> None:
    e0 = ROOT / "experiments" / "000_e0_leakage"
    ov = json.loads((e0 / "union_overlap/summary.json").read_text(encoding="utf-8"))
    rs = json.loads((e0 / "label_discordance_rnastralign/summary.json").read_text(encoding="utf-8"))
    bp = json.loads((e0 / "label_discordance_bpfold_archiveII/summary.json").read_text(encoding="utf-8"))

    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.2))
    labels = ["exact duplicate", "near-duplicate (k=15)", "internal duplicates"]
    values = [
        ov["exact_union_rate"],
        ov["near_union_rate"],
        ov["within_archiveii_duplicate_sequences"] / ov["n_archiveii"],
    ]
    colors = ["#4c72b0", "#55a868", "#c44e52"]
    bars = axes[0].bar(labels, [v * 100 for v in values], color=colors)
    axes[0].set_ylim(0, 100)
    axes[0].set_ylabel("% of ArchiveII")
    for b, v in zip(bars, values):
        axes[0].text(b.get_x() + b.get_width() / 2, v * 100 + 2, f"{v*100:.1f}%", ha="center", fontsize=9)
    axes[0].set_title(f"ArchiveII overlap with RNAStrAlign + TR0 (n={ov['n_archiveii']:,})", fontsize=10)
    axes[0].tick_params(axis="x", labelsize=8)

    names = ["vs RNAStrAlign", "vs bpRNA TR0"]
    f1s = [rs["mean_pair_f1"], bp["mean_pair_f1"]]
    # Dots, not bars: the axis starts at 0.85 and a truncated bar would
    # exaggerate the gap between the two corpora.
    for x, (value, colour) in enumerate(zip(f1s, ["#4c72b0", "#8172b3"])):
        axes[1].plot(x, value, "o", ms=10, color=colour)
        axes[1].text(x + 0.12, value, f"{value:.4f}", ha="left", va="center", fontsize=9)
    axes[1].set_xticks(range(len(names)))
    axes[1].set_xticklabels(names)
    axes[1].set_xlim(-0.4, 1.6)
    axes[1].set_ylim(0.85, 1.01)
    axes[1].set_ylabel("mean pair F1")
    axes[1].set_title(f"Label agreement of duplicates (n={rs['n_archiveii']:,}; y from 0.85)", fontsize=10)
    for ax in axes:
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)
        ax.grid(axis="y", color="#eeeeee", lw=0.8)
        ax.set_axisbelow(True)
    fig.suptitle("Benchmark contamination: overlap plus label identity", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    save(fig, out_dir / "fig1_overlap_and_labels")


def fig6(out_dir: Path) -> None:
    ec1 = json.loads((ROOT / "experiments/002_e1prime_runs/ec1_p1_p2.json").read_text(encoding="utf-8"))
    spec = json.loads((ROOT / "experiments/002_e1prime_runs/ec1b_specificity.json").read_text(encoding="utf-8"))
    fm = [(row["p1"], row["delta_f1"], row["n"]) for row in ec1["families"]]
    oh = [(row["p1"], row["delta_onehot"], row.get("n_common", 0)) for row in spec["families"]]

    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.4), sharex=True)
    for ax, pts, colour, rho, title in (
        (axes[0], fm, "#d62728", ec1["p1"]["rho"], "frozen FM head"),
        (axes[1], oh, "#7f7f7f", spec["rho_onehot"], "no-pretraining one-hot control"),
    ):
        xs = [p[0] for p in pts]
        ys = [p[1] for p in pts]
        sizes = [max(6, min(40, (p[2] or 1) ** 0.5 / 3)) for p in pts]
        ax.scatter(xs, ys, s=sizes, color=colour, alpha=0.55, linewidths=0)
        ax.axhline(0, color="#999999", lw=0.8, ls=":")
        ax.set_xscale("symlog", linthresh=1.0)
        ax.set_xlabel("Rfam full-region coverage proxy (p1, symlog; 0 included)")
        ax.set_title(f"{title} — rho = {rho:.3f}", fontsize=10)
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)
        ax.grid(color="#eeeeee", lw=0.8)
        ax.set_axisbelow(True)
    axes[0].set_ylabel("per-family delta F1 (head - RNAfold)")
    fig.suptitle(
        f"Coverage predicts accuracy for both heads: delta-rho = {spec['delta_rho']:.3f} < .15",
        fontsize=11,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    save(fig, out_dir / "fig6_coverage_and_control")

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
    fig1(out)
    fig6(out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
