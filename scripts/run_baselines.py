#!/usr/bin/env python3
"""Run and score an external secondary-structure predictor on one split.

The queue item is "BPfold, RNAstructure, EternaFold, MXfold2 on id_test and
test_withheld". Their command lines differ and some write files rather than
stdout, so the runner is deliberately tool-agnostic: it takes a command
template, collects stdout plus any text file the command leaves in its output
directory, and hands the whole blob to ``rff.analysis.dotbracket``. Predictions
that cannot be matched to a sequence are reported as missing, never guessed.

Scoring uses the pre-registered metric: family-macro base-pair F1 with a
family-block bootstrap interval, plus sequence-micro F1 as a secondary number.

Examples (fill the paths on the server, they are not installed locally):

    python scripts/run_baselines.py --split splits/s1prime_filtered/id_test.jsonl \
        --tool rnafold --out results/baselines/id_test \
        --cmd "RNAfold --noPS -d0 < {fasta} > rnafold.txt"

    python scripts/run_baselines.py --split ... --tool mxfold2 \
        --cmd "mxfold2 predict -f {fasta} -o mxfold2"
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rff.analysis.dotbracket import (  # noqa: E402
    ct_to_pairs,
    extract_predictions,
    parse_dotbracket,
)
from rff.analysis.embargo import assert_allowed  # noqa: E402
from rff.analysis.scoring import prf  # noqa: E402

TEXT_SUFFIXES = {".txt", ".out", ".db", ".dbn", ".dotbracket", ".ct", ".log", ".ps"}


def read_jsonl(path: Path) -> list[dict]:
    with open(path) as fh:
        return [json.loads(line) for line in fh if line.strip()]


def write_fasta(rows: list[dict], path: Path) -> None:
    with open(path, "w") as fh:
        for row in rows:
            fh.write(f">{row['id']}\n{row['seq']}\n")


def collect_text(stdout: str, outdir: Path | None, budget: int = 5_000_000) -> str:
    parts = [stdout]
    if outdir and outdir.exists():
        for path in sorted(outdir.rglob("*")):
            if not path.is_file() or path.suffix.lower() not in TEXT_SUFFIXES:
                continue
            if path.stat().st_size > budget:
                continue
            try:
                parts.append(path.read_text(errors="ignore"))
            except OSError:
                continue
    return "\n".join(parts)


def ct_predictions(outdir: Path | None, ids: set[str]) -> dict[str, set[tuple[int, int]]]:
    out: dict[str, set[tuple[int, int]]] = {}
    if not outdir or not outdir.exists():
        return out
    for path in outdir.rglob("*.ct"):
        stem = path.stem
        if stem in ids:
            out[stem] = ct_to_pairs(path.read_text(errors="ignore"))
    return out


def bootstrap_f1(per_family: dict[str, float], n_boot: int, seed: int) -> dict:
    fams = sorted(per_family)
    values = np.array([per_family[f] for f in fams]) if fams else np.array([])
    point = round(float(values.mean()), 4) if values.size else None
    if len(fams) < 3:
        return {"n_families": len(fams), "mean": point, "ci": None}
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(fams), (n_boot, len(fams)))
    means = values[idx].mean(axis=1)
    return {
        "n_families": len(fams),
        "mean": point,
        "ci": [round(float(np.percentile(means, 2.5)), 4),
               round(float(np.percentile(means, 97.5)), 4)],
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--split", required=True, help="jsonl with id/seq/family/pairs")
    ap.add_argument("--tool", required=True)
    ap.add_argument("--cmd", default="", help="shell command; {fasta} and {outdir} are substituted")
    ap.add_argument("--predictions", default="", help="parse this file instead of running --cmd")
    ap.add_argument("--outdir", default="", help="working/output directory for the tool")
    ap.add_argument("--out", required=True)
    ap.add_argument("--n-boot", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--keep", action="store_true", help="keep the temporary worksheet")
    ap.add_argument(
        "--confirmatory",
        action="store_true",
        help="required to read a test_novel path; also verifies the amendment hashes (A6.5)",
    )
    args = ap.parse_args()

    assert_allowed(
        [args.split, args.predictions] if args.predictions else [args.split],
        confirmatory=args.confirmatory,
        root=ROOT,
    )

    rows = read_jsonl(ROOT / args.split)
    seq_by_id = {r["id"]: r["seq"] for r in rows}
    workdir = Path(args.outdir) if args.outdir else Path(tempfile.mkdtemp(prefix="baseline_"))
    workdir.mkdir(parents=True, exist_ok=True)
    fasta = workdir / "input.fasta"
    write_fasta(rows, fasta)

    run_info = {"command": args.cmd, "returncode": None}
    if args.predictions:
        text = (ROOT / args.predictions).read_text(errors="ignore")
    else:
        if not args.cmd:
            raise SystemExit("either --cmd or --predictions is required")
        cmd = args.cmd.replace("{fasta}", str(fasta)).replace("{outdir}", str(workdir))
        proc = subprocess.run(cmd, shell=True, cwd=workdir, capture_output=True, text=True)
        run_info["returncode"] = proc.returncode
        run_info["stderr_tail"] = (proc.stderr or "")[-2000:]
        text = collect_text(proc.stdout or "", workdir)

    found = extract_predictions(text, seq_by_id)
    found_ct = ct_predictions(workdir, set(seq_by_id))

    per_family: dict[str, list[float]] = defaultdict(list)
    totals = {"tp": 0, "fp": 0, "fn": 0}
    by_band: dict[str, dict[str, int]] = defaultdict(lambda: {"tp": 0, "fp": 0, "fn": 0})
    missing: list[str] = []
    for row in rows:
        rid, length = row["id"], len(row["seq"])
        truth = {tuple(p) for p in row.get("pairs", [])}
        if rid in found_ct:
            pred = found_ct[rid]
        elif rid in found:
            try:
                pred = parse_dotbracket(found[rid], length)
            except ValueError:
                missing.append(rid)
                continue
        else:
            missing.append(rid)
            continue
        tp = len(pred & truth)
        fp = len(pred - truth)
        fn = len(truth - pred)
        per_family[row["family"]].append(prf(tp, fp, fn)["f1"])
        for k, v in (("tp", tp), ("fp", fp), ("fn", fn)):
            totals[k] += v
        band = "<80" if length < 80 else ("80-149" if length < 150 else ("150-299" if length < 300 else ">=300"))
        for k, v in (("tp", tp), ("fp", fp), ("fn", fn)):
            by_band[band][k] += v

    fam_mean = {f: float(np.mean(v)) for f, v in per_family.items() if v}
    micro = prf(totals["tp"], totals["fp"], totals["fn"])
    report = {
        "tool": args.tool,
        "split": args.split,
        "run": run_info,
        "n_sequences": len(rows),
        "n_scored": sum(len(v) for v in per_family.values()),
        "n_missing": len(missing),
        "missing_examples": missing[:20],
        "family_macro_f1": bootstrap_f1(fam_mean, args.n_boot, args.seed),
        "micro": micro,
        "by_length_band": {
            band: prf(v["tp"], v["fp"], v["fn"]) for band, v in sorted(by_band.items())
        },
    }

    out = ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"tool": args.tool, "scored": report["n_scored"],
                      "missing": report["n_missing"],
                      "family_macro_f1": report["family_macro_f1"]}))

    if not args.keep and not args.outdir:
        shutil.rmtree(workdir, ignore_errors=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
