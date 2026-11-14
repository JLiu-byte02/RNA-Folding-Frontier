"""The family diagnostics must regenerate the committed artifacts byte-for-byte."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
COHORT = ROOT / "experiments" / "002_e1prime_runs" / "cohort"
NAMES = (
    "family_delta_structure.json",
    "family_delta_concentration.json",
    "seed_variation.json",
)


def test_family_diagnostics_reproduce_committed_artifacts(tmp_path):
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts" / "family_diagnostics.py"),
            "--cohort",
            str(COHORT),
            "--out-dir",
            str(tmp_path),
        ],
        check=True,
        cwd=ROOT,
    )
    for name in NAMES:
        assert (tmp_path / name).read_bytes() == (COHORT / name).read_bytes(), name
