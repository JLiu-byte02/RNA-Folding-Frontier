import json

import numpy as np
import pytest

from rff.analysis.block_eval import (
    SeqScore,
    bootstrap_delta,
    bootstrap_mean,
    bootstrap_ratio,
    common_ids,
    family_macro,
    family_means,
    micro,
    score_dump,
    seed_average_family_means,
)


def dump_row(rid, family, truth, scored, length=30):
    return {
        "id": rid,
        "family": family,
        "length": length,
        "truth_v3": truth,
        "scored": {f"{i},{j}": p for (i, j), p in scored.items()},
    }


def write_dump(path, rows):
    path.write_text("".join(json.dumps(r) + "\n" for r in rows))


def test_zero_overlap_is_zero_not_vacuous():
    assert SeqScore("x", "f", tp=0, fp=5, fn=4).f1 == 0.0
    assert SeqScore("x", "f", tp=4, fp=0, fn=0).f1 == 1.0


def test_score_dump_uses_v3_truth_and_min_pairs(tmp_path):
    rows = [
        dump_row(
            "keep",
            "f1",
            [[0, 10], [1, 9], [2, 8], [3, 7]],
            {(0, 10): 0.9, (1, 9): 0.9, (4, 8): 0.9},
        ),
        dump_row("drop", "f1", [[0, 10], [1, 9], [2, 8]], {(0, 10): 0.9}),
    ]
    path = tmp_path / "a.jsonl"
    write_dump(path, rows)
    table = score_dump(path)
    assert set(table) == {"keep"}
    assert table["keep"].tp == 2
    assert table["keep"].fp == 1
    assert table["keep"].fn == 2


def test_family_macro_and_micro():
    table = {
        "a": SeqScore("a", "f1", tp=3, fp=0, fn=0),
        "b": SeqScore("b", "f1", tp=0, fp=1, fn=0),
        "c": SeqScore("c", "f2", tp=1, fp=1, fn=1),
    }
    assert family_macro(table) == pytest.approx(0.5)
    assert micro(table) == pytest.approx(8 / 11)


def test_seed_average_requires_common_ids():
    t1 = {
        "a": SeqScore("a", "f1", tp=4, fp=0, fn=0),
        "b": SeqScore("b", "f2", tp=0, fp=0, fn=4),
    }
    t2 = {
        "a": SeqScore("a", "f1", tp=0, fp=4, fn=0),
        "b": SeqScore("b", "f2", tp=4, fp=0, fn=0),
    }
    ids = common_ids(t1, t2)
    means = seed_average_family_means([t1, t2], ids)
    assert means["f1"] == pytest.approx(0.5)
    assert means["f2"] == pytest.approx(0.5)


def test_bootstrap_ratio_recovers_known_ratio():
    values = {
        "a": {"f1": 0.8, "f2": 0.7, "f3": 0.6, "f4": 0.5},
        "c": {"f1": 0.7, "f2": 0.6, "f3": 0.5, "f4": 0.4},
        "e": {"f1": 0.2, "f2": 0.2, "f3": 0.2, "f4": 0.2},
    }
    result = bootstrap_ratio(values, numerator="c", denominator_a="a", denominator_e="e")
    # (0.55 - 0.2) / (0.65 - 0.2) = 0.777...
    assert result["ratio"] == pytest.approx(7 / 9)
    assert result["lo"] <= result["ratio"] <= result["hi"]
    assert result["undefined_fraction"] == 0.0


def test_nonpositive_denominator_is_undefined():
    values = {
        "a": {"f1": 0.2, "f2": 0.2},
        "c": {"f1": 0.4, "f2": 0.4},
        "e": {"f1": 0.3, "f2": 0.3},
    }
    result = bootstrap_ratio(values, numerator="c", denominator_a="a", denominator_e="e")
    assert np.isnan(result["ratio"])
    assert result["undefined_fraction"] == 1.0
