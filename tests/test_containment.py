"""The contamination screen must fire on a planted near-duplicate.

A containment statistic that silently returns zero would pass every split, so
these tests assert both ends: a 1-nt mutant of a training sequence is caught, and
an unrelated sequence is not.
"""

from __future__ import annotations

import sys
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from leakage_filter import build_index, containment, kmer_codes  # noqa: E402

_RNG = np.random.default_rng(7)
_RANDOM = "".join(_RNG.choice(list("ACGU"), size=120))

TRAIN = [
    {"id": "t0", "seq": "GGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUU"},
    {"id": "t1", "seq": "ACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGU"},
    {"id": "t2", "seq": "UUUUGGGGCCCCAAAAGGGGCCCCAAAAGGGGCCCCAAAAGGGGCCCC"},
    {"id": "t3", "seq": _RANDOM},
]


def test_kmer_codes_are_unique_per_window():
    codes = kmer_codes("ACGUACGUACGU", 5)
    assert codes.size == 4  # ACGUA, CGUAC, GUACG, UACGU
    assert len(set(codes.tolist())) == codes.size


def test_identical_sequence_is_fully_contained():
    index, df, cut = build_index(TRAIN, 15, 0.05)
    raw, info, n_info = containment(TRAIN[0]["seq"], 15, index, df, cut)
    assert raw > 0.99
    if n_info >= 20:
        assert info > 0.99


def test_one_nt_mutant_is_still_detected():
    seq = list(_RANDOM)
    seq[60] = "A" if seq[60] != "A" else "C"
    mutant = "".join(seq)
    index, df, cut = build_index(TRAIN, 15, 0.05)
    _raw, info, n_info = containment(mutant, 15, index, df, cut)
    assert n_info >= 20
    # one substitution breaks at most 15 of the 106 windows
    assert info >= 0.75


def test_unrelated_sequence_is_not_contained():
    rng = np.random.default_rng(0)
    seq = "".join(rng.choice(list("ACGU"), size=120))
    index, df, cut = build_index(TRAIN, 15, 0.05)
    _raw, info, n_info = containment(seq, 15, index, df, cut)
    if n_info >= 20:
        assert info < 0.2
