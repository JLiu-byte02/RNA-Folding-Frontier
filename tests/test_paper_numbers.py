"""Every number in the Paper A tables must match a committed artifact.

The 2026-11-14 audits found numbers that had drifted from the JSON they claimed
to summarize. Hand-checking fixed that instance; this test makes the class of
error fail CI instead. It parses the markdown tables directly and compares them
with ``paired_*.json`` / ``band_*.json`` at the precision the paper prints.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT / "docs" / "paper_A_draft_2026-11-13.md"
COHORT = ROOT / "experiments" / "002_e1prime_runs" / "cohort"

LABEL = {
    "RNAfold": "rnafold",
    "EternaFold": "eternafold",
    "RNAstructure": "rnastructure",
    "MXfold2 (retrained)": "mxfold2",
    "head, seed 17": "head17",
    "head, seed 23": "head23",
    "head, seed 41": "head41",
    "head, 3-seed mean": "head_mean",
}

METHOD_ROW = re.compile(
    r"\| (RNAfold|EternaFold|RNAstructure|MXfold2 \(retrained\)|head, seed \d+|head, 3-seed mean) "
    r"\| ([.\d]+) \[([-.\d]+), ([-.\d]+)\] \| ([.\d]+) "
    r"\| ([.\d]+) \[([-.\d]+), ([-.\d]+)\] \| ([.\d]+) \|"
)
DELTA_ROW = re.compile(
    r"\| (head_mean|head, seed 17|seed17) - (\w+) "
    r"\| ([-+.\d]+) \[([-+.\d]+), ([-+.\d]+)\] \| [<.\d]+ "
    r"\| ([-+.\d]+) \[([-+.\d]+), ([-+.\d]+)\] \| [<.\d]+ \|"
)


def _paired(split: str) -> dict:
    return json.loads((COHORT / f"paired_{split}.json").read_text(encoding="utf-8"))


def _method(table: dict, key: str) -> dict:
    for name, value in table["methods"].items():
        if name.lower() == key.lower():
            return value
    raise KeyError(key)


def test_method_table_matches_paired_artifacts():
    text = PAPER.read_text(encoding="utf-8")
    rows = METHOD_ROW.findall(text)
    assert len(rows) == 8, f"expected 8 method rows, parsed {len(rows)}"
    for row in rows:
        name = row[0]
        key = LABEL[name]
        for split, want in (("test_withheld", row[1:5]), ("id_test", row[5:9])):
            ref = _method(_paired(split), key)
            got = (ref["family_macro_f1"], *ref["family_macro_f1_ci95"], ref["micro_f1"])
            assert tuple(round(float(x), 4) for x in want) == tuple(round(x, 4) for x in got), (
                name,
                split,
                want,
                got,
            )


def test_paired_delta_table_matches_paired_artifacts():
    text = PAPER.read_text(encoding="utf-8")
    rows = DELTA_ROW.findall(text)
    assert len(rows) == 6, f"expected 6 delta rows, parsed {len(rows)}"
    for a_label, b_label, d_w, lo_w, hi_w, d_i, lo_i, hi_i in rows:
        a = {"head_mean": "head_mean", "head, seed 17": "head17", "seed17": "head17"}[a_label]
        b = {"RNAfold": "rnafold", "EternaFold": "eternafold", "RNAstructure": "rnastructure", "MXfold2": "mxfold2"}[b_label]
        for split, want in (("test_withheld", (d_w, lo_w, hi_w)), ("id_test", (d_i, lo_i, hi_i))):
            entry = next(
                p for p in _paired(split)["paired"]
                if p["a"].lower() == a.lower() and p["b"].lower() == b.lower()
            )
            got = (entry["delta"], *entry["ci95"])
            assert tuple(round(float(x), 4) for x in want) == tuple(round(x, 4) for x in got), (
                a,
                b,
                split,
                want,
                got,
            )


def test_band_table_matches_band_artifact():
    text = PAPER.read_text(encoding="utf-8")
    band = json.loads((COHORT / "band_test_withheld.json").read_text(encoding="utf-8"))
    by_band = {tuple(b["band"]): b for b in band["bands"] if b.get("n_sequences", 0)}
    pattern = re.compile(
        r"\| (\d+)-(\d+) \| ([\d,]+) / (\d+) "
        r"\| ([+\-.0-9]+) \[([+\-.0-9]+), ([+\-.0-9]+)\] "
        r"\| ([+\-.0-9]+) \[([+\-.0-9]+), ([+\-.0-9]+)\] "
        r"\| ([+\-.0-9]+) \[([+\-.0-9]+), ([+\-.0-9]+)\] \|"
    )
    rows = pattern.findall(text)
    assert len(rows) == 4, f"expected 4 band rows, parsed {len(rows)}"
    for row in rows:
        lo, hi = int(row[0]), int(row[1])
        entry = by_band[(lo, hi)]
        assert int(row[2].replace(",", "")) == entry["n_sequences"]
        assert int(row[3]) == entry["n_families"]
        for idx, tool in enumerate(("RNAfold", "EternaFold", "RNAstructure")):
            paired = next(p for p in entry["paired"] if p["b"] == tool)
            want = tuple(round(float(x), 3) for x in row[4 + 3 * idx : 7 + 3 * idx])
            got = (round(paired["delta"], 3), *(round(x, 3) for x in paired["ci95"]))
            assert want == got, (tool, (lo, hi), want, got)
