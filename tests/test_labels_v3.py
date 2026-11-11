"""The frozen labels_v3 artifact is a committed, hash-checked record."""

from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_labels_v3_artifact_and_summary_agree():
    summary = json.loads(
        (ROOT / "experiments/011_labels_v3/labels_v3_summary.json").read_text(
            encoding="utf-8"
        )
    )
    path = ROOT / summary["pairs_file"]
    blob = path.read_bytes()
    assert hashlib.sha256(blob).hexdigest() == summary["pairs_sha256"]
    assert len(blob) == summary["pairs_bytes"]

    labels = json.loads(gzip.decompress(blob))
    assert len(labels) == summary["n_sequences"] == 92_781
    assert all(value is not None for value in labels.values())


def test_eval_subset_is_repairable_under_v3():
    summary = json.loads(
        (ROOT / "experiments/011_labels_v3/labels_v3_summary.json").read_text(
            encoding="utf-8"
        )
    )
    audit = summary["eval_subset_audit"]
    assert audit["n_records"] == 2_372
    assert audit["n_v1_empty_truth"] == 463
    assert audit["n_v3_parse_errors"] == 0
    assert audit["n_v3_lt4_pairs"] == 0
    assert audit["n_excluded_by_lt4_rule"] == 0
    assert audit["lt4_gate_exceeded_10pct"] is False
