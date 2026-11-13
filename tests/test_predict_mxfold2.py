"""Contract tests for the MXfold2 inference wrapper.

The wrapper is the only place that turns MXfold2's stdout into the
``{"id", "dbn"}`` rows the external-baseline scorer consumes, so a parsing
regression would silently corrupt the four-method table.  These tests pin the
stdout contract and the embargo refusal without needing the vendor package.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]

_spec = importlib.util.spec_from_file_location(
    "predict_mxfold2", ROOT / "scripts" / "predict_mxfold2.py"
)
mx = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(mx)


def test_parse_stdout_reads_triples_with_scores():
    text = ">a\nGGGAAACCC\n(((...)))\t(-4.2)\n>b\nACGU\n....\t(0.0)\n"
    assert mx._parse_stdout(text) == [("a", "(((...)))"), ("b", "....")]


def test_parse_stdout_rejects_a_truncated_record():
    text = ">a\nGGGAAACCC\n(((...)))\t(-4.2)\n>b\nACGU\n"
    with pytest.raises(ValueError):
        mx._parse_stdout(text)


def test_parse_stdout_rejects_missing_header():
    with pytest.raises(ValueError):
        mx._parse_stdout("GGG\n...\n")


def test_embargo_refusal_names_the_path():
    with pytest.raises(SystemExit) as exc:
        mx.refuse_embargoed(Path("splits/s1prime_filtered/test_novel.jsonl"))
    assert "test_novel" in str(exc.value)


def test_embargo_allows_ordinary_splits():
    mx.refuse_embargoed(Path("splits/s1prime_filtered/id_test.jsonl"))
