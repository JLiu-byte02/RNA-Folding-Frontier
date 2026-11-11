"""Contract tests for the ViennaRNA soft-constraint decoder.

ViennaRNA is not installed on the development machine, so these tests drive the
wrapper with a stub ``RNA`` module. They pin the three things that would
silently corrupt results on the server: index base, sign, and which cells get a
pseudo-energy.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rff.decoders.nn_dp import NEG  # noqa: E402
from rff.decoders.vienna_soft import (  # noqa: E402
    ViennaSoftDecoder,
    parse_dotbracket,
    zero_scores,
)


class FakeFoldCompound:
    def __init__(self, seq: str, structure: str, energy: float):
        self.seq = seq
        self.structure = structure
        self.energy = energy
        self.bp_calls: list[tuple[int, int, float]] = []

    def sc_add_bp(self, i: int, j: int, e: float) -> None:
        self.bp_calls.append((i, j, e))

    def mfe(self):
        return self.structure, self.energy


class FakeRNA:
    def __init__(self, structure="((..))", energy=-5.0):
        self.structure = structure
        self.energy = energy
        self.made: list[FakeFoldCompound] = []

    def fold_compound(self, seq: str) -> FakeFoldCompound:
        fc = FakeFoldCompound(seq, self.structure, self.energy)
        self.made.append(fc)
        return fc


def test_parse_dotbracket():
    assert parse_dotbracket("((..))") == {(0, 5), (1, 4)}
    assert parse_dotbracket("......") == set()


def test_scores_map_to_negative_pseudo_energy():
    rna = FakeRNA()
    dec = ViennaSoftDecoder(scale=2.0, rna_module=rna)
    seq = "GGGAAACCC" if len("GGGAAACCC") == 9 else "GGGAAACCC"
    scores = np.full((len(seq), len(seq)), NEG)
    scores[0, 8] = 3.0
    scores[1, 7] = -1.0
    pairs, energy = dec.decode(seq, scores)
    assert pairs == [(0, 5), (1, 4)]
    assert energy == -5.0
    calls = rna.made[0].bp_calls
    assert calls == [(0, 8, -6.0), (1, 7, 2.0)]


def test_disallowed_cells_are_not_constrained():
    rna = FakeRNA()
    dec = ViennaSoftDecoder(rna_module=rna)
    seq = "ACGUACGUAC"
    scores = np.full((len(seq), len(seq)), NEG)
    dec.decode(seq, scores)
    assert rna.made[0].bp_calls == []
    assert dec.last["n_constraints"] == 0


def test_zero_scores_are_the_plain_fold():
    rna = FakeRNA(structure="(((....)))", energy=-9.3)
    dec = ViennaSoftDecoder(rna_module=rna)
    seq = "GGGAAACCCGGG"
    dec.decode(seq, zero_scores(seq))
    # every allowed cell is added, but with energy 0: a no-op that keeps the
    # Turner model in charge, which is the parity requirement
    assert all(e == 0.0 for _i, _j, e in rna.made[0].bp_calls)
    structure, energy = dec.fold_plain(seq)
    assert structure == "(((....)))"
    assert energy == -9.3
