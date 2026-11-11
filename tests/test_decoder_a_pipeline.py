"""End-to-end check of score -> soft constraints -> subgradient -> better structures.

Uses RNAfold MFE labels, so it validates the wiring rather than learning from
real data. It is slow-ish for a unit test (~1 s), which is why the sequence
count is small; the full run lives in `scripts/train_decoder_a_toy.py` and its
record in `experiments/003_decoder_a/`.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]

pytest.importorskip("RNA")


def test_toy_pipeline_improves_f1(tmp_path):
    out = tmp_path / "toy.json"
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts" / "train_decoder_a_toy.py"),
            "--n-sequences", "8",
            "--epochs", "4",
            "--out", str(out),
        ],
        check=True,
        capture_output=True,
    )
    report = json.loads(out.read_text())
    assert report["f1_final"] > report["f1_initial"] + 0.1
