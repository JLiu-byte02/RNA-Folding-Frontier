#!/usr/bin/env python3
"""End-to-end validation of the decoder-A training pipeline on synthetic labels.

PI Round 4c asked for a toy run "using random head outputs to check gradient
direction, then a tiny real run". This is that toy run, and it needs no GPU, no
backbone and no cluster data:

* labels are `RNAfold` MFE structures on random sequences, so a correct pipeline
  *can* reach F1 ≈ 1 — the labels are representable by the decoder;
* the head is simulated by a noisy score matrix, which is what a badly trained
  head looks like;
* training is the structured max-margin step in `rff.training.maxmargin`, and
  decoding is `ViennaSoftDecoder` (Turner + `sc_add_bp`).

What it proves: scores → soft constraints → MFE → subgradient → better structures
is wired in the right direction. What it does not prove: anything about real
data, real backbones, or the calibration of `--energy-scale` on log-odds.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import RNA  # noqa: E402

from rff.analysis.scoring import prf  # noqa: E402
from rff.decoders.nn_dp import NEG, mask_allowed  # noqa: E402
from rff.decoders.vienna_soft import ViennaSoftDecoder  # noqa: E402
from rff.training.maxmargin import indicator, structured_step  # noqa: E402


def mfe_pairs(seq: str) -> set[tuple[int, int]]:
    structure, _ = RNA.fold(seq)
    table = RNA.ptable(structure)
    return {(i - 1, table[i] - 1) for i in range(1, len(seq) + 1) if table[i] > i}


def decode_pairs(decoder):
    def decode(seq, scores):
        return decoder.decode(seq, scores)

    return decode


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-sequences", type=int, default=40)
    ap.add_argument("--epochs", type=int, default=12)
    ap.add_argument("--margin", type=float, default=2.0)
    ap.add_argument("--lr", type=float, default=0.15)
    ap.add_argument("--noise", type=float, default=0.7)
    ap.add_argument("--scale", type=float, default=1.0)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--min-pairs", type=int, default=3)
    ap.add_argument("--out", default="experiments/003_decoder_a/toy_maxmargin.json")
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    seqs = []
    while len(seqs) < args.n_sequences:
        length = int(rng.integers(20, 41))
        seq = "".join(rng.choice(list("ACGU"), length))
        truth = mfe_pairs(seq)
        if len(truth) >= args.min_pairs:
            seqs.append((seq, sorted(truth)))

    theta = {}
    for seq, _truth in seqs:
        mask = mask_allowed(seq)
        t = np.full((len(seq), len(seq)), NEG)
        t[mask] = rng.normal(0.0, args.noise, mask.sum())
        theta[seq] = t

    decoder = ViennaSoftDecoder(scale=args.scale)
    decode = decode_pairs(decoder)

    def evaluate() -> dict:
        tp = fp = fn = 0
        for seq, truth in seqs:
            pred, _ = decoder.decode(seq, theta[seq])
            p, t = set(pred), set(truth)
            tp += len(p & t)
            fp += len(p - t)
            fn += len(t - p)
        return prf(tp, fp, fn)

    curve = [{"epoch": 0, **{k: round(v, 4) for k, v in evaluate().items()}}]
    for epoch in range(1, args.epochs + 1):
        lr = args.lr / (1.0 + 0.25 * (epoch - 1))
        for seq, truth in seqs:
            structured_step(theta[seq], seq, truth, mask_allowed(seq), decode,
                            margin=args.margin, lr=lr)
        if epoch % 2 == 0 or epoch == args.epochs:
            curve.append({"epoch": epoch, **{k: round(v, 4) for k, v in evaluate().items()}})

    report = {
        "settings": {
            "n_sequences": len(seqs),
            "epochs": args.epochs,
            "margin": args.margin,
            "lr": args.lr,
            "noise": args.noise,
            "energy_scale": args.scale,
            "vienna_version": getattr(RNA, "__version__", "unknown"),
            "labels": "RNAfold MFE",
        },
        "curve": curve,
        "f1_initial": curve[0]["f1"],
        "f1_final": curve[-1]["f1"],
    }
    report["improved"] = bool(report["f1_final"] > report["f1_initial"] + 0.2)

    out = ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: report[k] for k in ("f1_initial", "f1_final", "improved")}))
    return 0 if report["improved"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
