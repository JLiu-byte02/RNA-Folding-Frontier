"""WUSS consensus-structure parsing, including letter-encoded pseudoknots.

The v1 seed builder (``scripts/build_rfam_index.py``) keeps only ``()[]{}<>`` and
maps every other consensus character to ``.``. Rfam marks many pseudoknots with
letter pairs (``AAA...aaa``), so v1 silently drops them — measured in Amendment 6
§A6.2 and left in place there because rebuilding the splits would change frozen
hashes. PI Round 4g §2 asked instead for a separate ``labels_v2`` artifact, which
this module supplies.

WUSS conventions used here:

* ``(`` ``[`` ``{`` ``<`` open and their mirrors close, as before;
* an uppercase letter opens a pseudoknot helix and the matching lowercase letter
  closes it (``A``…``a``, ``B``…``b``, …);
* ``.`` ``,`` ``:`` ``_`` ``-`` are unpaired (they may be used to grade
  accessibility) and ``~`` is a gap;
* anything else is treated as unpaired rather than guessed at.

This module is separate from :func:`rff.data.io.pairs_from_dbn` on purpose: the
v1 parser is frozen inside an already-hashed pipeline, and quietly teaching it a
new notation would make old labels unreproducible.
"""

from __future__ import annotations

import gzip
import json
import string
from pathlib import Path

from .io import pairs_from_dbn

OPEN_TO_CLOSE = {"(": ")", "[": "]", "{": "}", "<": ">"}
CLOSE_TO_OPEN = {v: k for k, v in OPEN_TO_CLOSE.items()}
LETTERS = string.ascii_uppercase


def project_wuss(ss_cons: str, aligned: str) -> str:
    """Project a consensus structure onto one gapped sequence row.

    Columns where ``aligned`` is a gap are dropped, exactly as the v1 projection
    does, so the result lines up with the ungapped sequence.
    """
    out: list[str] = []
    for col, ch in enumerate(ss_cons):
        if col >= len(aligned):
            break
        if aligned[col] in "-._~":
            continue
        if ch in OPEN_TO_CLOSE or ch in CLOSE_TO_OPEN or ch in LETTERS or ch in LETTERS.lower():
            out.append(ch)
        else:
            out.append(".")
    return "".join(out)


def pairs_from_wuss(text: str) -> list[tuple[int, int]]:
    """Pairs from a WUSS string, honouring bracket levels and letter pairs.

    Raises ``ValueError`` on an unmatched closing symbol, which would otherwise
    turn a corrupt consensus row into a plausible-looking label. Unmatched
    *opening* symbols are tolerated and dropped, because they are what a
    truncated consensus row looks like and dropping them loses no pair that the
    row actually closed.
    """
    stacks: dict[str, list[int]] = {c: [] for c in OPEN_TO_CLOSE}
    letters: dict[str, list[int]] = {c: [] for c in LETTERS}
    pairs: list[tuple[int, int]] = []
    for idx, ch in enumerate(text):
        if ch in stacks:
            stacks[ch].append(idx)
        elif ch in CLOSE_TO_OPEN:
            opener = CLOSE_TO_OPEN[ch]
            if not stacks[opener]:
                raise ValueError(f"unmatched '{ch}' at position {idx}")
            pairs.append((stacks[opener].pop(), idx))
        elif ch in LETTERS:
            letters[ch].append(idx)
        elif ch in LETTERS.lower():
            key = ch.upper()
            if not letters[key]:
                raise ValueError(f"unmatched '{ch}' at position {idx}")
            pairs.append((letters[key].pop(), idx))
    return sorted(pairs)


def canonical_pairs(
    seq: str, pairs: list[tuple[int, int]]
) -> set[tuple[int, int]]:
    """Keep AU/UA/GC/CG/GU/UG pairs, the label policy frozen in Amendment 6 §A6.3."""
    allowed = {("A", "U"), ("U", "A"), ("G", "C"), ("C", "G"), ("G", "U"), ("U", "G")}
    return {
        (i, j)
        for i, j in pairs
        if 0 <= i < len(seq) and 0 <= j < len(seq) and (seq[i], seq[j]) in allowed
    }


def v1_pairs(seq: str, dbn_v1: str) -> set[tuple[int, int]]:
    """Pairs implied by the frozen v1 projection (brackets only)."""
    return canonical_pairs(seq, pairs_from_dbn(dbn_v1))


def load_labels_v2_delta(path: str | Path) -> dict[str, set[tuple[int, int]]]:
    """Read ``labels_v2_delta.json.gz`` into ``{seed id: canonical pair set}``.

    Keys are ``<RF accession>|<seed name>``; a split entry has id
    ``rfam:<RF accession>|<seed name>``, so strip the ``rfam:`` prefix when
    applying the delta.
    """
    with gzip.open(Path(path), "rt", encoding="utf-8") as fh:
        raw = json.load(fh)
    return {key: {(int(i), int(j)) for i, j in pairs} for key, pairs in raw.items()}


__all__ = [
    "canonical_pairs",
    "load_labels_v2_delta",
    "pairs_from_wuss",
    "project_wuss",
    "v1_pairs",
]
