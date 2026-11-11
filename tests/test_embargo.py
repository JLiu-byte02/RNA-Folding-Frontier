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
    blob_at_head,
    check_amendments,
    sha256_bytes,
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


def _git_repo(tmp_path: Path, text: str = "frozen text\n") -> Path:
    """A throwaway repo with one committed amendment; returns its path."""
    repo = tmp_path / "repo"
    repo.mkdir()

    def git(*args):
        subprocess.run(
            ["git", "-c", "user.name=t", "-c", "user.email=t@local", *args],
            cwd=repo, check=True, capture_output=True,
        )

    git("init", "-q")
    (repo / "AMENDMENT_X.md").write_bytes(text.encode())
    git("add", "AMENDMENT_X.md")
    git("commit", "-q", "-m", "freeze")
    return repo


def _blob_lock(repo: Path, path: str = "AMENDMENT_X.md") -> Path:
    digest = sha256_bytes(blob_at_head(repo, path))
    lock = repo / "lock.json"
    lock.write_text(
        json.dumps({"hash_basis": "git-blob", "files": {path: {"blob_sha256": digest}}}),
        encoding="utf-8",
    )
    return lock


def test_blob_basis_ignores_an_uncommitted_edit_but_refuses_it(tmp_path):
    """The lock is over HEAD bytes; a dirty worktree must not slip through."""
    repo = _git_repo(tmp_path)
    lock = _blob_lock(repo)
    assert_allowed([str(repo / "test_novel.jsonl")], confirmatory=True,
                   root=repo, lockfile=lock.relative_to(repo))
    (repo / "AMENDMENT_X.md").write_text("edited after the freeze\n", encoding="utf-8")
    with pytest.raises(EmbargoViolation) as err:
        check_amendments(repo, lock.relative_to(repo))
    assert "working tree" in str(err.value)


def test_blob_basis_refuses_a_line_ending_rewrite(tmp_path):
    """The exact autocrlf failure mode: same text, different bytes."""
    repo = _git_repo(tmp_path)
    lock = _blob_lock(repo)
    (repo / "AMENDMENT_X.md").write_bytes(b"frozen text\r\n")
    with pytest.raises(EmbargoViolation) as err:
        check_amendments(repo, lock.relative_to(repo))
    assert "line-ending" in str(err.value)


def test_blob_basis_refuses_when_head_moved(tmp_path):
    repo = _git_repo(tmp_path)
    lock = _blob_lock(repo)
    (repo / "AMENDMENT_X.md").write_text("a later amendment\n", encoding="utf-8")
    subprocess.run(
        ["git", "-c", "user.name=t", "-c", "user.email=t@local", "commit", "-qam", "edit"],
        cwd=repo, check=True, capture_output=True,
    )
    with pytest.raises(EmbargoViolation) as err:
        check_amendments(repo, lock.relative_to(repo))
    assert "blob" in str(err.value)


def test_untracked_locked_file_is_refused(tmp_path):
    repo = _git_repo(tmp_path)
    lock = _blob_lock(repo)
    (repo / "EXTRA.md").write_text("never committed\n", encoding="utf-8")
    lock.write_text(
        json.dumps({"hash_basis": "git-blob",
                    "files": {"EXTRA.md": {"blob_sha256": sha256_file(repo / "EXTRA.md")}}}),
        encoding="utf-8",
    )
    with pytest.raises(EmbargoViolation) as err:
        check_amendments(repo, lock.relative_to(repo))
    assert "not present in HEAD" in str(err.value)


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
