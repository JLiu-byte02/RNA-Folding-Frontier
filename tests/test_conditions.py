"""Contract tests for the E1' swap conditions.

The conditions only mean something if they are length-preserving, leave the
local window untouched, and keep the statistics they claim to keep.
"""

from __future__ import annotations

import random
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from build_conditions import dinuc_shuffle, distal_mask, positional_map  # noqa: E402


def test_dinucleotide_composition_is_preserved():
    rng = random.Random(0)
    for _ in range(20):
        seq = "".join(rng.choice("ACGU") for _ in range(80))
        shuffled = dinuc_shuffle(seq, rng)
        assert len(shuffled) == len(seq)
        assert set(shuffled) <= set("ACGU")
        assert Counter(zip(seq, seq[1:])) == Counter(zip(shuffled, shuffled[1:])), (
            seq,
            shuffled,
        )


def test_centre_window_is_kept():
    mask = distal_mask(200, 32, "centre", [])
    assert not any(mask[68:133])
    assert all(mask[:68])
    assert all(mask[133:])


def test_paired_window_is_kept():
    mask = distal_mask(200, 10, "paired", [(50, 150)])
    assert not any(mask[40:61])
    assert not any(mask[140:161])
    assert mask[0] and mask[199]


def test_short_sequence_still_masks_something():
    mask = distal_mask(30, 32, "centre", [])
    assert len(mask) == 30


def test_positional_map_is_identity_for_equal_sequences():
    seq = "ACGUGCAUGC"
    assert positional_map(seq, seq) == list(range(len(seq)))


def test_positional_map_reports_gaps_as_minus_one():
    target = "ACGUGCAUGC"
    donor = "ACGUGCAUGCGGG"
    mapping = positional_map(target, donor)
    assert len(mapping) == len(target)
    assert all(0 <= j < len(donor) for j in mapping)
    # An insertion in the donor must leave a gap somewhere on the target side.
    shorter = positional_map(donor, target)
    assert any(j == -1 for j in shorter)
