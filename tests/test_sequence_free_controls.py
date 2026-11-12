"""Controls for the template-prior inputs."""

from __future__ import annotations

import random
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from build_sequence_free_controls import dinucleotide_shuffle  # noqa: E402


def test_shuffle_preserves_dinucleotides_endpoints_and_length():
    seq = "GGGAAACCCUUGGAAC"
    out = dinucleotide_shuffle(seq, random.Random(17))
    assert len(out) == len(seq)
    assert out[0] == seq[0] and out[-1] == seq[-1]
    assert Counter(zip(seq, seq[1:])) == Counter(zip(out, out[1:]))


def test_shuffle_is_reproducible():
    seq = "ACGUACGUACGU"
    assert dinucleotide_shuffle(seq, random.Random(3)) == dinucleotide_shuffle(
        seq, random.Random(3)
    )
