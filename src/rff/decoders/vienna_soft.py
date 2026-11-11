"""Decoder A via ViennaRNA soft constraints ("Turner plus learned corrections").

PI Round 4c rejected re-implementing Turner 2004 by hand: the parameter model is
weeks of parity work and nothing in our claims needs to own it. Instead the
learned scores enter ViennaRNA as per-pair pseudo-energies
(``fold_compound.sc_add_bp``), so with all corrections set to zero the decoder
reproduces ``RNAfold`` exactly and parity is structural rather than earned.

Conventions that must be checked against the installed version before any
number is quoted (marked [CHECK] in the Round 4c reply):

* ``sc_add_bp(i, j, e)`` takes **1-based** indices — the docstring of the
  matrix overload states ``1 <= i < j <= n``, and the single-pair overload
  follows the same convention (verified against ViennaRNA 2.8.0 by
  `scripts/check_vienna_parity.py`: a 0-based call silently constrained the
  wrong pair and 59 of 60 bonuses failed to attract the pair they were given);
  it adds ``e`` in kcal/mol every time the pair is formed and returns 0 when the
  constraint is rejected;
* the sign: a positive pseudo-energy destabilises a pair, so a score that is
  "higher is better" maps to ``e = -scale * score``;
* the useful magnitude of ``scale``: log-odds and kcal/mol are different units,
  and a scale that is too large makes pseudo-energies dominate the physical
  model. Calibrate on the server (one sequence, corrections on/off) rather than
  guessing here; the default 1.0 is a placeholder, not a calibrated value.

The class imports ``RNA`` lazily so the module and its tests run on machines
without ViennaRNA.
"""

from __future__ import annotations

import numpy as np

from rff.decoders.nn_dp import NEG, pairs_to_dotbracket


def parse_dotbracket(db: str) -> set[tuple[int, int]]:
    stack: list[int] = []
    pairs: set[tuple[int, int]] = set()
    for i, ch in enumerate(db):
        if ch == "(":
            stack.append(i)
        elif ch == ")":
            if not stack:
                raise ValueError(f"unbalanced ')' at {i}")
            pairs.add((stack.pop(), i))
    if stack:
        raise ValueError("unbalanced '(' in ViennaRNA output")
    return pairs


class ViennaSoftDecoder:
    """Wrap ``RNA.fold_compound`` + ``sc_add_bp`` as a score-matrix decoder.

    ``decode(seq, scores)`` returns ``(pairs, energy)`` where ``scores`` follows
    the same convention as :mod:`rff.decoders.nn_dp`: higher is better, ``NEG``
    marks disallowed or absent cells. Cells whose score is ``NEG`` are simply not
    given a pseudo-energy; ViennaRNA then decides them from the Turner model,
    which is why a zero matrix reproduces plain ``RNAfold``.
    """

    def __init__(self, scale: float = 1.0, rna_module=None):
        self.scale = float(scale)
        self._rna = rna_module

    @property
    def rna(self):
        if self._rna is None:
            import RNA  # imported here so CPU-only tests do not need it

            self._rna = RNA
        return self._rna

    def decode(self, seq: str, scores: np.ndarray) -> tuple[list[tuple[int, int]], float]:
        n = len(seq)
        if scores.shape != (n, n):
            raise ValueError(f"scores must be {n}x{n}, got {scores.shape}")
        fc = self.rna.fold_compound(seq)
        n_constraints = 0
        n_rejected = 0
        for i in range(n):
            for j in range(i + 1, n):
                s = float(scores[i, j])
                if s <= NEG:
                    continue
                ok = fc.sc_add_bp(i + 1, j + 1, -self.scale * s)
                if ok:
                    n_constraints += 1
                else:
                    n_rejected += 1
        structure, energy = fc.mfe()
        self.last = {
            "n_constraints": n_constraints,
            "n_rejected": n_rejected,
            "energy": float(energy),
        }
        return sorted(parse_dotbracket(structure)), float(energy)

    def fold_plain(self, seq: str) -> tuple[str, float]:
        """Turner-only structure, no constraints: the parity control."""
        fc = self.rna.fold_compound(seq)
        structure, energy = fc.mfe()
        return structure, float(energy)


def zero_scores(seq: str) -> np.ndarray:
    """Score matrix whose decode must equal plain RNAfold (parity check)."""
    n = len(seq)
    return np.zeros((n, n))


__all__ = ["ViennaSoftDecoder", "parse_dotbracket", "pairs_to_dotbracket", "zero_scores"]
