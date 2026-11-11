"""Dot-bracket and CT parsing, plus best-effort extraction from tool output.

Every external predictor (RNAfold, MXfold2, EternaFold, RNAstructure, BPfold)
prints its structure differently. Rather than one brittle parser per tool, this
module recognises the three shapes that cover all of them:

1. ``<id><TAB><dot-bracket>`` or ``<id> <dot-bracket>``;
2. the RNAfold layout — a sequence line followed by its dot-bracket line;
3. RNAstructure ``.ct`` files, parsed separately.

Anything that matches none of these is reported as missing rather than guessed
at: a silently mis-assigned structure would be far worse than a missing row.
"""

from __future__ import annotations

import re

BRACKETS = {"(": ")", "[": "]", "{": "}", "<": ">"}
OPEN = set(BRACKETS)
CLOSE = set(BRACKETS.values())
DB_CHARS = set(".()-[]{}<>:,_~") | set("ACGUacgu")


def parse_dotbracket(db: str, length: int | None = None) -> set[tuple[int, int]]:
    """Pairs from a dot-bracket string; supports the four nested bracket types."""
    if length is not None and len(db) != length:
        raise ValueError(f"dot-bracket length {len(db)} != sequence length {length}")
    stacks: dict[str, list[int]] = {c: [] for c in OPEN}
    pairs: set[tuple[int, int]] = set()
    for idx, ch in enumerate(db):
        if ch in OPEN:
            stacks[ch].append(idx)
        elif ch in CLOSE:
            opener = next(o for o, c in BRACKETS.items() if c == ch)
            if not stacks[opener]:
                raise ValueError(f"unmatched '{ch}' at {idx}")
            pairs.add((stacks[opener].pop(), idx))
        elif ch in ".:,_-~":
            continue
        else:
            raise ValueError(f"unexpected character {ch!r} at {idx}")
    leftover = {o: len(s) for o, s in stacks.items() if s}
    if leftover:
        raise ValueError(f"unmatched opening brackets: {leftover}")
    return pairs


def is_dotbracket(line: str, length: int) -> bool:
    # all-dots is a valid structure (no pairs), so no bracket is required
    return len(line) == length and all(ch in DB_CHARS for ch in line)


def parse_ct(text: str) -> dict[int, int]:
    """RNAstructure CT: ``i base i j k`` with 1-based indices, j = 0 if unpaired."""
    lines = [ln.split() for ln in text.strip().splitlines() if ln.strip()]
    if not lines:
        return {}
    n = int(lines[0][0])
    partner: dict[int, int] = {}
    for row in lines[1 : n + 1]:
        i, j = int(row[0]), int(row[4])
        if j:
            partner[i - 1] = j - 1
    return partner


def ct_to_pairs(text: str) -> set[tuple[int, int]]:
    partner = parse_ct(text)
    return {(i, j) for i, j in partner.items() if i < j}


SEQ_RE = re.compile(r"^[ACGUacguNnTt]+$")


def extract_predictions(text: str, seq_by_id: dict[str, str]) -> dict[str, str]:
    """Map id -> dot-bracket from tool output, using only unambiguous matches."""
    counts: dict[str, int] = {}
    for seq in seq_by_id.values():
        counts[seq.upper()] = counts.get(seq.upper(), 0) + 1
    by_seq = {
        seq.upper(): rid
        for rid, seq in seq_by_id.items()
        if counts[seq.upper()] == 1
    }
    out: dict[str, str] = {}
    pending_id: str | None = None
    pending_seq: str | None = None
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            pending_id = line[1:].split()[0]
            pending_seq = None
            continue

        tokens = line.split()
        if len(tokens) >= 2:
            head, tail = tokens[0], tokens[-1]
            if head in seq_by_id:
                hit = next(
                    (t for t in tokens[1:] if is_dotbracket(t, len(seq_by_id[head]))),
                    None,
                )
                if hit is not None:
                    out[head] = hit
                    continue
            if tail in seq_by_id and is_dotbracket(head, len(seq_by_id[tail])):
                out[tail] = head
                continue

        if pending_seq is not None and is_dotbracket(line, len(pending_seq)):
            rid = by_seq.get(pending_seq.upper())
            if rid is not None:
                out[rid] = line
            pending_seq = None
            continue

        if SEQ_RE.match(line):
            pending_seq = line
            if pending_id is not None and pending_id in seq_by_id:
                by_seq.setdefault(line.upper(), pending_id)
            continue

        pending_seq = None
    return out
