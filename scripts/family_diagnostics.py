"""Family-level diagnostics derived from the paired cohort artifacts.

Three analyses that the paper quotes but that were originally produced by
throw-away scripts:

* ``family_delta_structure.json`` -- correlation of the head-tool delta with
  log family size and mean length, with family-bootstrap CIs;
* ``family_delta_concentration.json`` -- how concentrated the gains and losses
  are across families;
* ``seed_variation.json`` -- per-seed deltas and a t interval across the three
  training seeds.

Everything is read from ``paired_<split>.json``; no predictions or labels are
re-read, so the numbers are exactly those behind Table T7.

    python scripts/family_diagnostics.py --cohort experiments/002_e1prime_runs/cohort
"""

from __future__ import annotations

import argparse
import json
import math
import statistics
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SPLITS = ("test_withheld", "id_test")
TOOLS = ("rnafold", "eternafold", "rnastructure")
SEEDS = ("head17", "head23", "head41")


def load(coh: Path, split: str) -> dict:
    return json.loads((coh / f"paired_{split}.json").read_text(encoding="utf-8"))


def method(entry: dict, key: str) -> dict:
    for name, value in entry["methods"].items():
        if name.lower() == key.lower():
            return value
    raise KeyError(key)


def family_lengths(coh: Path, split: str) -> tuple[dict[str, float], dict[str, int]]:
    fam_len: dict[str, list[int]] = {}
    for line in (coh / f"{split}.jsonl").read_text(encoding="utf-8").splitlines():
        if line.strip():
            row = json.loads(line)
            fam_len.setdefault(row["family"], []).append(int(row["length"]))
    return ({f: sum(v) / len(v) for f, v in fam_len.items()},
            {f: len(v) for f, v in fam_len.items()})


def structure(coh: Path, out: Path) -> None:
    rng = np.random.default_rng(20261114)
    result = {}
    for split in SPLITS:
        d = load(coh, split)
        mean_len, count = family_lengths(coh, split)
        head = method(d, "head_mean")["per_family_mean_f1"]
        fams = sorted(head)
        ln = np.log(np.array([count[f] for f in fams], dtype=float))
        L = np.log(np.array([mean_len[f] for f in fams]))
        rows = {}
        for tool in TOOLS:
            t = method(d, tool)["per_family_mean_f1"]
            delta = np.array([head[f] - t[f] for f in fams])
            idx = rng.integers(0, len(fams), size=(5000, len(fams)))
            rhos = np.empty(len(idx))
            meds = np.empty(len(idx))
            for k, ii in enumerate(idx):
                rhos[k] = np.corrcoef(ln[ii], delta[ii])[0, 1]
                meds[k] = np.median(delta[ii])
            rows[tool] = {
                "rho_log_n": round(float(np.corrcoef(ln, delta)[0, 1]), 4),
                "rho_log_n_ci95": [round(float(np.percentile(rhos, 2.5)), 4),
                                   round(float(np.percentile(rhos, 97.5)), 4)],
                "rho_log_len": round(float(np.corrcoef(L, delta)[0, 1]), 4),
                "median_delta": round(float(np.median(delta)), 4),
                "median_delta_ci95": [round(float(np.percentile(meds, 2.5)), 4),
                                      round(float(np.percentile(meds, 97.5)), 4)],
                "mean_delta": round(float(delta.mean()), 4),
                "frac_pos": round(float((delta > 0).mean()), 4),
            }
        result[split] = {"n_families": len(fams), "n_boot": 5000, "seed": 20261114, "vs": rows}
    (out / "family_delta_structure.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print("wrote family_delta_structure.json")


def concentration(coh: Path, out: Path) -> None:
    result = {}
    for split in SPLITS:
        d = load(coh, split)
        head = method(d, "head_mean")["per_family_mean_f1"]
        fams = sorted(head)
        rows = {}
        for tool in ("rnafold", "eternafold"):
            t = method(d, tool)["per_family_mean_f1"]
            delta = np.array([head[f] - t[f] for f in fams])
            neg, pos = delta[delta < 0], delta[delta > 0]
            order = np.argsort(delta)
            rows[tool] = {
                "mean_delta": round(float(delta.mean()), 4),
                "median_delta": round(float(np.median(delta)), 4),
                "n_negative": int(len(neg)),
                "n_positive": int(len(pos)),
                "negative_mass": round(float(-neg.sum()), 3),
                "positive_mass": round(float(pos.sum()), 3),
                "worst10_mass": round(float(-delta[order[:10]].sum()), 3),
                "worst10_share_of_loss": round(float(-delta[order[:10]].sum() / -neg.sum()), 4),
                "best10_mass": round(float(delta[order[-10:]].sum()), 3),
                "best10_share_of_gain": round(float(delta[order[-10:]].sum() / pos.sum()), 4),
            }
        result[split] = {"n_families": len(fams), "vs": rows}
    (out / "family_delta_concentration.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print("wrote family_delta_concentration.json")


def seed_variation(coh: Path, out: Path) -> None:
    result = {}
    for split in SPLITS:
        d = load(coh, split)
        head = method(d, "head_mean")["per_family_mean_f1"]
        fams = sorted(head)
        rows = {}
        for tool in TOOLS:
            t = method(d, tool)["per_family_mean_f1"]
            per_seed = {}
            for seed in SEEDS:
                h = method(d, seed)["per_family_mean_f1"]
                per_seed[seed] = sum(h[f] - t[f] for f in fams) / len(fams)
            vals = list(per_seed.values())
            mean = statistics.fmean(vals)
            sd = statistics.stdev(vals)
            half = 4.303 * sd / math.sqrt(len(vals))
            rows[tool] = {
                "per_seed": {k: round(v, 4) for k, v in per_seed.items()},
                "seed_mean": round(mean, 4),
                "seed_sd": round(sd, 4),
                "t95": [round(mean - half, 4), round(mean + half, 4)],
            }
        result[split] = rows
    (out / "seed_variation.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print("wrote seed_variation.json")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", default="experiments/002_e1prime_runs/cohort")
    ap.add_argument("--out-dir", default=None, help="defaults to the cohort directory")
    args = ap.parse_args()
    coh = Path(args.cohort)
    if not coh.is_absolute():
        coh = ROOT / coh
    out = Path(args.out_dir) if args.out_dir else coh
    if not out.is_absolute():
        out = ROOT / out
    out.mkdir(parents=True, exist_ok=True)
    structure(coh, out)
    concentration(coh, out)
    seed_variation(coh, out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
