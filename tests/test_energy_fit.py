"""Metric contract for the (s,b) calibration fit."""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from fit_energy_calibration import family_macro_f1  # noqa: E402


def test_family_macro_uses_per_sequence_f1_not_pooled_counts():
    items = [
        {"family": "f1", "truth": {(0, 10), (1, 9)}},
        {"family": "f1", "truth": {(0, 10)}},
        {"family": "f2", "truth": {(0, 10), (1, 9)}},
    ]
    preds = {
        0: {(0, 10), (1, 9)},  # f1 sequence 1: 1.0
        1: {(2, 8)},            # f1 sequence 2: 0.0
        2: {(0, 10), (1, 9)},  # f2: 1.0
    }
    # Per-family means: f1 = 0.5, f2 = 1.0 -> macro = 0.75.
    # Pooled counts would give f1 = 8/9 and f2 = 1.0 -> 17/18.
    assert family_macro_f1(preds, items, [0, 1, 2]) == 0.75
