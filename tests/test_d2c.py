"""Checks for the exploratory FP-helix energy helper."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
pytest.importorskip("RNA")

from d2c_helix_energy import isolated_energy, stack_pairs  # noqa: E402


def test_stack_pairs_and_isolated_energy():
    seq = "GGGAAACCC"
    pairs = stack_pairs(seq, 0, 8)
    assert pairs == [(0, 8), (1, 7), (2, 6)]
    assert isolated_energy(seq, pairs) < 0


def test_noncanonical_pair_has_no_stack():
    assert stack_pairs("AAAAAAAAA", 0, 8) == []
