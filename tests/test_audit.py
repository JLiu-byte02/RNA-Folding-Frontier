from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from rff.data.audit import exact_duplicates, near_duplicates, summarize  # noqa: E402
from rff.data.io import Record, pairs_from_dbn  # noqa: E402


def test_pairs_from_dbn_nested():
    assert pairs_from_dbn("((..))") == [(0, 5), (1, 4)]


def test_exact_duplicates_across_splits():
    a = Record(rid="a", seq="ACGU", pairs=[(0, 3)])
    b = Record(rid="b", seq="ACGU", pairs=[(0, 3)])
    c = Record(rid="c", seq="AAAA", pairs=[])
    report = exact_duplicates({"train": [a, c], "test": [b]})
    assert report["n_cross_split_duplicate_sequences"] == 1


def test_near_duplicates_flags_overlap():
    ref = [Record(rid="r1", seq="ACGU" * 10, pairs=[])]
    query = [Record(rid="q1", seq="ACGU" * 10 + "A", pairs=[])]
    report = near_duplicates(ref, query, k=5, min_containment=0.9)
    assert report["n_hits"] == 1


def test_summarize_counts():
    records = [Record(rid="a", seq="ACGU", pairs=[(0, 3)]), Record(rid="b", seq="AA", pairs=[])]
    summary = summarize(records)
    assert summary["n"] == 2
    assert summary["pairs_total"] == 1
