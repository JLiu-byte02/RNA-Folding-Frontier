"""Contract tests for the D1/D2 scoring primitives."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rff.analysis.scoring import (  # noqa: E402
    Region,
    average_precision,
    best_threshold,
    distance_bin,
    pair_scores,
    prf,
    summarise,
)
from rff.analysis.labels import pair_prf  # noqa: E402


def test_prf_edges():
    assert prf(0, 0, 0)["f1"] == 1.0
    assert prf(1, 0, 0)["f1"] == 1.0
    assert prf(0, 3, 3)["f1"] == 0.0


def test_pair_prf_zero_overlap_is_not_vacuous():
    assert pair_prf([], [])["f1"] == 1.0
    assert pair_prf([(0, 5)], [])["f1"] == 0.0
    assert pair_prf([], [(0, 5)])["f1"] == 0.0
    assert pair_prf([(0, 5)], [(1, 6)])["f1"] == 0.0
    assert pair_prf([(0, 5)], [(0, 5)])["f1"] == 1.0


def test_old_zero_overlap_convention_is_rejected():
    # The pre-fix implementation returned 1.0 here because precision+recall==0.
    result = pair_prf([(0, 5)], [(1, 6)])
    assert result["tp"] == 0
    assert result["f1"] == 0.0


def test_region_partition_covers_all_pairs():
    reg = Region(100, 165)
    assert reg.where(110, 150) == "both_kept"
    assert reg.where(110, 300) == "one_kept"
    assert reg.where(10, 300) == "neither_kept"
    assert reg.where(10, 99) == "neither_kept"


def test_distance_bins_are_ordered():
    assert distance_bin(0, 10) == "d<32"
    assert distance_bin(0, 32) == "d32-128"
    assert distance_bin(0, 129) == "d>128"


def test_average_precision_perfect_and_worst():
    y = np.array([1, 1, 0, 0])
    s = np.array([0.9, 0.8, 0.2, 0.1])
    assert average_precision(y, s) == 1.0
    s_bad = np.array([0.1, 0.2, 0.8, 0.9])
    assert average_precision(y, s_bad) < 0.6


def test_oracle_threshold_beats_fixed_cut_when_scale_shifted():
    # Labels are separable, but all scores are pushed up: a 0.5 cut predicts
    # everything, while the oracle cut recovers the perfect split.
    y = np.array([1, 1, 0, 0])
    s = np.array([0.99, 0.98, 0.97, 0.96])
    assert prf(int(((s >= 0.5) & (y == 1)).sum()), int(((s >= 0.5) & (y == 0)).sum()), 0)["f1"] < 0.7
    cut, f1 = best_threshold(y, s)
    assert f1 == 1.0
    assert cut > 0.5


def test_summarise_splits_by_region_and_distance():
    rows = [
        {
            "truth": [[100, 150], [170, 200]],
            "scored": {"100,150": 0.9, "170,200": 0.9, "120,190": 0.7},
            "length": 210,
        }
    ]
    out = summarise(rows, kept_lo=95, kept_hi=160)
    assert out["tp"] == 2
    assert out["fp"] == 1
    assert out["by_region"]["both_kept"]["tp"] == 1
    assert out["by_region"]["neither_kept"]["tp"] == 1
    assert out["by_region"]["one_kept"]["fp"] == 1
    assert out["mean_predicted_pairs"] == 3.0
    assert out["auprc"] > 0.8


def _load_ec2_module():
    import importlib.util

    path = ROOT / "scripts" / "score_ec2_gate.py"
    spec = importlib.util.spec_from_file_location("score_ec2_gate_under_test", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_ec2_gate_decodes_each_sequence_once():
    """The scoring loop must decode one scale per sequence, not one per arm.

    Job 15627 burned 56 CPU-h because ``_decode_item`` looped over per-arm
    scales, re-decoding every sequence ``len(scales)`` times.  ``_init`` now
    binds a single per-sequence scale and reports the planned cost so the
    caller (and the sbatch timing gate) can assert linearity.
    """
    import pytest

    module = _load_ec2_module()
    items = [{"seq": "ACGU", "truth": set(), "scored": [], "z": 0.0, "log_len": 1.0}] * 5
    plan = module._init(items, [0.1] * 5)
    assert plan == {"n_items": 5, "decodes_per_item": 1, "total_decodes": 5}
    assert module._ARMS == [0.1] * 5

    with pytest.raises(SystemExit) as excinfo:
        module._init(items, [0.1] * 3)
    assert "3 entries for 5 items" in str(excinfo.value)
