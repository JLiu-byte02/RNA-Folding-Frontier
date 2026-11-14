"""The paper figures must be regenerable from committed artifacts.

The rendering needs matplotlib, which the fresh-clone CI installs but the
Windows workstation does not have; the test skips there and runs in CI.
"""

from __future__ import annotations

import subprocess
import sys

import pytest

pytest.importorskip("matplotlib.pyplot")

from pathlib import Path  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]


def test_figures_render_from_artifacts(tmp_path):
    out = tmp_path / "figures"
    subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "make_paper_figures.py"), "--out-dir", str(out)],
        check=True,
        cwd=ROOT,
    )
    for name in (
        "fig1_overlap_and_labels.png",
        "fig1_overlap_and_labels.svg",
        "fig2_label_projection.png",
        "fig2_label_projection.svg",
        "fig6_coverage_and_control.png",
        "fig6_coverage_and_control.svg",
        "fig7_four_method_cohort.png",
        "fig7_four_method_cohort.svg",
    ):
        path = out / name
        assert path.exists() and path.stat().st_size > 1000, name
