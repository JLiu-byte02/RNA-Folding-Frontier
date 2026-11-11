"""Contract tests for H2.3 (block-boundary false positives).

Two things must hold: the window geometry matches ``build_block_cache.py``, and
the statistic actually separates a synthetic boundary artefact from a synthetic
long-range-only difference. The second is the one that keeps the number honest:
long-range extra FPs are boundary by construction, which is why the script
reports distance bands and a distance-matched null.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from d2b_block_boundary import (  # noqa: E402
    block_window,
    distance_band,
    mutually_native,
)


def test_block_window_matches_block_cache_flank_rule():
    # block=16, flank=32: middle block [32, 48) -> window [0, 80)
    assert block_window(40, 300, 16, 32) == (0, 80)
    # first block clips at the sequence start
    assert block_window(0, 300, 16, 32) == (0, 48)
    # last partial block clips at the sequence end
    assert block_window(295, 300, 16, 32) == (256, 300)


def test_mutual_nativeness_boundary_is_symmetric():
    length, block, flank = 300, 16, 32
    assert mutually_native(40, 70, length, block, flank)
    assert not mutually_native(40, 90, length, block, flank)
    assert mutually_native(90, 40, length, block, flank) == mutually_native(
        40, 90, length, block, flank
    )
    # beyond the native span nothing can be mutually native
    span = block + 2 * flank
    assert not mutually_native(10, 10 + span + 1, length, block, flank)


def test_distance_bands_split_at_the_native_span():
    span = 80
    assert distance_band(10, span) == "d<=half_span"
    assert distance_band(60, span) == "half_span<d<=span"
    assert distance_band(81, span) == "d>span"


def _write_dump(path: Path, rows: list[dict]) -> None:
    with open(path, "w") as fh:
        for row in rows:
            fh.write(json.dumps(row) + "\n")


def _run(tmp_path, a_pairs, c_pairs, seq, block=16, flank=32):
    """Run the CLI on one sequence; returns the parsed report."""
    import subprocess

    split = tmp_path / "split.jsonl"
    _write_dump(split, [{"id": "s1", "seq": seq, "pairs": []}])
    a = tmp_path / "a.jsonl"
    c = tmp_path / "c.jsonl"
    common = {"id": "s1", "family": "f1", "length": len(seq), "truth": []}
    _write_dump(a, [{**common, "scored": {f"{i},{j}": 0.9 for i, j in a_pairs}}])
    _write_dump(c, [{**common, "scored": {f"{i},{j}": 0.9 for i, j in c_pairs}}])
    out = tmp_path / "out.json"
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts" / "d2b_block_boundary.py"),
            "--a-dump", str(a),
            "--c-dump", str(c),
            "--split", str(split),
            "--block", str(block),
            "--flank", str(flank),
            "--out", str(out),
            "--n-boot", "200",
        ],
        check=True,
        capture_output=True,
    )
    return json.loads(out.read_text())


def test_local_extra_pair_is_not_counted_as_boundary(tmp_path):
    seq = ("ACGU" * 100)[:300]
    # a single extra FP that is inside one block's native window
    report = _run(tmp_path, [], [(40, 60)], seq)
    assert report["n_extra_fp"] == 1
    assert report["H2_3"]["extra_boundary_frac"] == 0.0
    assert report["H2_3"]["underpowered"] is True


def test_long_range_extra_pair_is_boundary_by_construction(tmp_path):
    seq = ("ACGU" * 100)[:300]
    report = _run(tmp_path, [], [(20, 220)], seq)
    assert report["n_extra_fp"] == 1
    assert report["H2_3"]["extra_boundary_frac"] == 1.0
    band = report["by_distance_band"]["d>span"]
    assert band["extra_n"] == 1
    # the null is distance-matched, so it is boundary too: no separation here
    assert band["ratio_obs_null"] in (None, 1.0)


def test_pairs_present_in_both_conditions_are_not_extra(tmp_path):
    seq = ("ACGU" * 100)[:300]
    report = _run(tmp_path, [(20, 220)], [(20, 220)], seq)
    assert report["n_extra_fp"] == 0
    assert report["H2_3"]["extra_boundary_frac"] is None


def test_true_pairs_are_excluded_from_extra(tmp_path):
    seq = ("ACGU" * 100)[:300]
    import subprocess

    split = tmp_path / "split.jsonl"
    _write_dump(split, [{"id": "s1", "seq": seq, "pairs": []}])
    a = tmp_path / "a.jsonl"
    c = tmp_path / "c.jsonl"
    _write_dump(a, [{"id": "s1", "family": "f1", "length": len(seq), "truth": [[20, 220]], "scored": {}}])
    _write_dump(c, [{"id": "s1", "family": "f1", "length": len(seq), "truth": [[20, 220]],
                     "scored": {"20,220": 0.9}}])
    out = tmp_path / "out.json"
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts" / "d2b_block_boundary.py"),
            "--a-dump", str(a), "--c-dump", str(c), "--split", str(split),
            "--out", str(out), "--n-boot", "100",
        ],
        check=True,
        capture_output=True,
    )
    report = json.loads(out.read_text())
    assert report["n_extra_fp"] == 0
    assert report["reference_rates"]["true_positive_c"]["n"] == 1
    assert report["reference_rates"]["true_positive_c"]["frac"] == 1.0
