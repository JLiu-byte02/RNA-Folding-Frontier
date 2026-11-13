"""Regression tests for common-cohort external scoring."""

from __future__ import annotations

import gzip
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import cohort_paired_eval as cpe  # noqa: E402


SEQ = "ACGUGCAUGCAU"
TRUTH = [[0, 11], [1, 10], [2, 9], [3, 8]]


def _jsonl(path: Path, rows):
    with path.open("w", encoding="utf-8", newline="\n") as fh:
        for row in rows:
            fh.write(json.dumps(row) + "\n")


def _labels(path: Path, ids):
    with gzip.open(path, "wt", encoding="utf-8") as fh:
        json.dump({rid[5:]: TRUTH for rid in ids}, fh)


def _split(path: Path, ids):
    _jsonl(
        path,
        [
            {"id": rid, "seq": SEQ, "family": "RF1", "pairs": TRUTH}
            for rid in ids
        ],
    )


def _cohort(path: Path, ids):
    _jsonl(path, [{"id": rid, "family": "RF1"} for rid in ids])


def test_pairs_prediction_is_scored_on_the_common_cohort(tmp_path, monkeypatch):
    rid = "rfam:RF00001|seedA"
    split, cohort, labels, preds = (tmp_path / n for n in ("split.jsonl", "cohort.jsonl", "labels.json.gz", "pred.jsonl"))
    out = tmp_path / "paired.json"
    _split(split, [rid])
    _cohort(cohort, [rid])
    _labels(labels, [rid])
    _jsonl(preds, [{"id": rid, "pairs": TRUTH}])
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "cohort_paired_eval.py",
            "--split", str(split),
            "--cohort", str(cohort),
            "--labels", str(labels),
            "--method", f"m={preds}",
            "--min-true-pairs", "4",
            "--n-boot", "50",
            "--out", str(out),
        ],
    )
    assert cpe.main() == 0
    result = json.loads(out.read_text(encoding="utf-8"))
    assert result["n_sequences"] == 1
    assert result["methods"]["m"]["family_macro_f1"] == 1.0
    assert result["methods"]["m"]["tp"] == 4


def test_missing_prediction_in_common_cohort_is_an_error(tmp_path, monkeypatch):
    ids = ["rfam:RF00001|seedA", "rfam:RF00001|seedB"]
    split, cohort, labels, preds = (tmp_path / n for n in ("split.jsonl", "cohort.jsonl", "labels.json.gz", "pred.jsonl"))
    _split(split, ids)
    _cohort(cohort, ids)
    _labels(labels, ids)
    _jsonl(preds, [{"id": ids[0], "pairs": TRUTH}])
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "cohort_paired_eval.py",
            "--split", str(split),
            "--cohort", str(cohort),
            "--labels", str(labels),
            "--method", f"m={preds}",
            "--min-true-pairs", "4",
            "--n-boot", "50",
            "--out", str(tmp_path / "paired.json"),
        ],
    )
    with pytest.raises(SystemExit, match="lack at least one prediction"):
        cpe.main()
