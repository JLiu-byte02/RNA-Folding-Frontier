"""The embargo guard must fail closed (Amendment 6 §A6.5)."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rff.analysis.embargo import (  # noqa: E402
    EmbargoViolation,
    assert_allowed,
    check_amendments,
    sha256_file,
)


def _fixture(tmp_path: Path, tamper: bool = False):
    f = tmp_path / "AMENDMENT_X.md"
    f.write_text("frozen text\n", encoding="utf-8")
    locked = sha256_file(f)
    if tamper:
        f.write_text("frozen text, but edited\n", encoding="utf-8")
    lockfile = tmp_path / "lock.json"
    lockfile.write_text(json.dumps({"files": {"AMENDMENT_X.md": locked}}), encoding="utf-8")
    return f, lockfile


def test_ordinary_paths_are_never_blocked(tmp_path):
    f, lock = _fixture(tmp_path)
    assert_allowed([str(tmp_path / "id_test.jsonl")], confirmatory=False,
                   root=tmp_path, lockfile=lock.relative_to(tmp_path))


def test_test_novel_without_flag_is_refused(tmp_path):
    f, lock = _fixture(tmp_path)
    with pytest.raises(EmbargoViolation):
        assert_allowed([str(tmp_path / "test_novel.jsonl")], confirmatory=False,
                       root=tmp_path, lockfile=lock.relative_to(tmp_path))


def test_flag_with_matching_hashes_is_allowed(tmp_path):
    f, lock = _fixture(tmp_path)
    assert_allowed([str(tmp_path / "test_novel.jsonl")], confirmatory=True,
                   root=tmp_path, lockfile=lock.relative_to(tmp_path))


def test_flag_with_edited_amendment_is_refused(tmp_path):
    f, lock = _fixture(tmp_path, tamper=True)
    with pytest.raises(EmbargoViolation) as err:
        assert_allowed([str(tmp_path / "test_novel.jsonl")], confirmatory=True,
                       root=tmp_path, lockfile=lock.relative_to(tmp_path))
    assert "changed since their hashes" in str(err.value)


def test_missing_locked_file_is_refused_not_ignored(tmp_path):
    lock = tmp_path / "lock.json"
    lock.write_text(json.dumps({"files": {"GONE.md": "0" * 64}}), encoding="utf-8")
    with pytest.raises(EmbargoViolation):
        check_amendments(tmp_path, lock.relative_to(tmp_path))


def test_the_real_lockfile_matches_the_real_amendments():
    """If this fails, an amendment was edited after its hash was quoted."""
    verified = check_amendments(ROOT, Path("state/amendments.lock.json"))
    assert verified["experiments/001_e1prime/AMENDMENT_4.md"].startswith("08d6771a")


CONFIRMATORY_SCRIPTS = [
    (
        "h24_decoder_compare.py",
        ["--dump", "dumps/test_novel.jsonl", "--split", "splits/test_novel.jsonl",
         "--out", "h24.json"],
    ),
    (
        "run_baselines.py",
        ["--split", "splits/test_novel.jsonl", "--tool", "rnafold",
         "--predictions", "preds.txt", "--out", "p5.json"],
    ),
]


@pytest.mark.parametrize("script,args", CONFIRMATORY_SCRIPTS)
def test_confirmatory_scripts_refuse_test_novel_without_the_flag(script, args):
    """End to end: the guard must fire before any file is read."""
    proc = subprocess.run(
        [sys.executable, str(ROOT / "scripts" / script), *args],
        cwd=ROOT, capture_output=True, text=True,
    )
    assert proc.returncode != 0
    assert "embargoed" in (proc.stderr + proc.stdout).lower()


@pytest.mark.parametrize("script,args", CONFIRMATORY_SCRIPTS)
def test_confirmatory_scripts_pass_the_guard_when_allowed(script, args):
    """Same invocation plus --confirmatory: the refusal must be gone.

    These paths do not exist, so the script still fails afterwards; what this
    checks is that the failure is no longer the embargo, i.e. the guard is not a
    blanket refusal.
    """
    proc = subprocess.run(
        [sys.executable, str(ROOT / "scripts" / script), *args, "--confirmatory"],
        cwd=ROOT, capture_output=True, text=True,
    )
    assert "embargoed" not in (proc.stderr + proc.stdout).lower()
