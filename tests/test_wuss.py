"""WUSS parsing, including the letter pairs the v1 builder drops."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rff.data.wuss import (  # noqa: E402
    canonical_pairs,
    load_labels_v2_delta,
    pairs_from_wuss_columns,
    pairs_from_wuss,
    project_pairs,
    project_wuss,
    v1_pairs,
)


def test_letter_pairs_are_parsed():
    # letters stack like any other helix level: innermost pair first
    assert pairs_from_wuss("AAA...aaa") == [(0, 8), (1, 7), (2, 6)]
    assert pairs_from_wuss("AB....ab") == [(0, 6), (1, 7)]


def test_bracket_levels_and_letters_coexist():
    # a '(' level crossing a letter level, which is exactly what WUSS levels are
    # for: the two helices are independent, not nested
    text = "(AA..[..)..]..aa"
    assert pairs_from_wuss(text) == [(0, 8), (1, 15), (2, 14), (5, 11)]


def test_unmatched_closing_letter_is_refused():
    with pytest.raises(ValueError):
        pairs_from_wuss("..a..")


def test_unmatched_opening_is_dropped_not_guessed():
    assert pairs_from_wuss("AAA..") == []


def test_projection_keeps_letters_and_drops_gap_columns():
    ss = "AAA...aaa"
    aligned = "AC-AC-ACA"  # three gap columns
    assert project_wuss(ss, aligned) == "AA..aaa"


def test_projection_maps_unknown_symbols_to_dot():
    assert project_wuss("A,:_x-", "AAAAAA") == "A...x."


def test_column_pairs_keep_alignment_coordinates():
    assert pairs_from_wuss_columns("(AA..aa)") == [(0, 7), (1, 6), (2, 5)]


def test_project_pairs_drops_a_pair_if_either_end_is_gapped():
    # columns 0 and 5 are the outer pair; column 1/4 are the inner pair
    assert project_pairs("((..))", "A-ACG-") == []


def test_project_pairs_never_fabricates_a_pair_after_a_gap():
    # v2 projects the characters first, obtaining "(..)" -> (0, 3), even though
    # column 0 was paired with the gapped column 5 and column 4's partner was
    # gapped. Pair-level projection must return no pair.
    assert project_wuss("((..))", "A-ACG-") == "(..)"
    assert pairs_from_wuss(project_wuss("((..))", "A-ACG-")) == [(0, 3)]
    assert project_pairs("((..))", "A-ACG-") == []


def test_project_pairs_keeps_only_pairs_with_two_ungapped_columns():
    # outer pair survives (columns 0,7), inner pair loses column 1
    assert project_pairs("((....))", "A-ACGAAC") == [(0, 6)]


def test_canonical_filter_keeps_gc_and_gu_only():
    seq = "GAUC"  # 0 G, 1 A, 2 U, 3 C
    pairs = [(0, 3), (1, 2)]  # GC ok, AU ok
    assert canonical_pairs(seq, pairs) == {(0, 3), (1, 2)}
    assert canonical_pairs("AAAA", [(0, 3)]) == set()


def test_v1_drops_letter_pairs_that_v2_keeps():
    seq = "GGGAAACCC"
    # v1 projection maps every letter to '.', so the same consensus row yields
    # no pairs at all; v2 recovers the three GC pairs
    v1_projection = "".join(ch if ch in "()[]{}<>" else "." for ch in "AAA...aaa")
    assert v1_pairs(seq, v1_projection) == set()
    recovered = canonical_pairs(seq, pairs_from_wuss("AAA...aaa"))
    assert recovered == {(0, 8), (1, 7), (2, 6)}


def test_delta_gzip_round_trips_to_pair_sets(tmp_path):
    import gzip
    import json

    path = tmp_path / "delta.json.gz"
    payload = json.dumps({"RF00001|seq1": [[0, 6], [1, 5]]}).encode()
    path.write_bytes(gzip.compress(payload, mtime=0))
    loaded = load_labels_v2_delta(path)
    assert loaded == {"RF00001|seq1": {(0, 6), (1, 5)}}
