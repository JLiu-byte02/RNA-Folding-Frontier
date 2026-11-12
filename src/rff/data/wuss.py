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
GAP_CHARS = frozenset("-._~")


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
    return sorted(pairs_from_wuss_columns(text))


def pairs_from_wuss_columns(text: str) -> list[tuple[int, int]]:
    """Pairs in *alignment-column* coordinates.

    This is the parser underlying :func:`pairs_from_wuss`; it is exposed because
    pair-level projection must map columns to sequence coordinates before any
    pair is discarded. Parsing a projected string instead can fabricate a pair
    when a gap removes one helix endpoint and two unrelated partners become
    adjacent in the projected string.
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


def project_pairs(ss_cons: str, aligned: str) -> list[tuple[int, int]]:
    """Project consensus pairs onto one gapped row, pair by pair.

    PI Round 4h §3: parse ``SS_cons`` in alignment-column coordinates, map both
    endpoints to ungapped sequence coordinates, and drop a pair if either
    alignment column is a gap. This differs from parsing the character-projected
    string: when a gap removes one endpoint, the latter can pair the surviving
    endpoint with a different closing column and silently invent a base pair.

    ``ss_cons`` and ``aligned`` must have equal length. A corrupt consensus row
    raises ``ValueError`` rather than yielding a partially guessed label.
    """
    if len(ss_cons) != len(aligned):
        raise ValueError(
            f"alignment length mismatch: SS_cons={len(ss_cons)} row={len(aligned)}"
        )
    col_to_seq: dict[int, int] = {}
    for col, ch in enumerate(aligned):
        if ch in GAP_CHARS:
            continue
        col_to_seq[col] = len(col_to_seq)

    projected: list[tuple[int, int]] = []
    for i, j in pairs_from_wuss_columns(ss_cons):
        if i in col_to_seq and j in col_to_seq:
            projected.append((col_to_seq[i], col_to_seq[j]))
    return sorted(projected)


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


def load_labels_v3(path: str | Path) -> dict[str, set[tuple[int, int]] | None]:
    """Read the full v3 replacement map, preserving explicit parse failures."""
    with gzip.open(Path(path), "rt", encoding="utf-8") as fh:
        raw = json.load(fh)
    return {
        key: None if pairs is None else {(int(i), int(j)) for i, j in pairs}
        for key, pairs in raw.items()
    }


__all__ = [
    "canonical_pairs",
    "GAP_CHARS",
    "load_labels_v2_delta",
    "load_labels_v3",
    "pairs_from_wuss_columns",
    "pairs_from_wuss",
    "project_pairs",
    "project_wuss",
    "v1_pairs",
]
