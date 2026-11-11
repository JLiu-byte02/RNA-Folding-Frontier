"""Checks against the real ViennaRNA package, skipped when it is not installed.

The stub tests pin the wiring; these pin the physics. They are deliberately
small (a few short sequences) so they can run in CI as well as on the cluster.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

RNA = pytest.importorskip("RNA")

from rff.decoders.nn_dp import NEG, mask_allowed  # noqa: E402
from rff.decoders.vienna_soft import ViennaSoftDecoder  # noqa: E402

SEQS = ["GGGAAACCC", "GCGCUUCGCCGC", "AAAAAAAAAAAAAAAAAAAA", "GCCGGAUCCGGC"]


def test_zero_scores_reproduce_rnafold():
    dec = ViennaSoftDecoder()
    for seq in SEQS:
        ref_struct, ref_energy = RNA.fold(seq)
        pairs, energy = dec.decode(seq, np.zeros((len(seq), len(seq))))
        table = RNA.ptable(ref_struct)
        ref = {(i - 1, table[i] - 1) for i in range(1, len(seq) + 1) if table[i] > i}
        assert set(pairs) == ref
        assert energy == pytest.approx(ref_energy, abs=1e-6)


def test_bonus_attracts_the_asked_pair():
    dec = ViennaSoftDecoder(scale=1.0)
    seq = "GCGCUUCGCCGC"
    mask = mask_allowed(seq)
    ii, jj = np.where(mask)
    i, j = int(ii[-1]), int(jj[-1])
    scores = np.full((len(seq), len(seq)), NEG)
    scores[i, j] = 20.0
    pairs, _ = dec.decode(seq, scores)
    assert (i, j) in set(pairs)
