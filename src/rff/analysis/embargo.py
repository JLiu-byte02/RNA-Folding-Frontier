"""Access guard for the confirmatory (test_novel) evaluation set.

Amendment 6 §A6.5. The confirmatory scripts must refuse a path under
``test_novel`` unless two independent conditions hold:

1. the caller passed ``--confirmatory``;
2. every pre-registration file listed in ``state/amendments.lock.json`` still
   hashes to the recorded value.

The second condition is what makes the flag meaningful: a flag can be typed at
any time, but it cannot make an edited amendment match its old hash. The lockfile
is committed, so editing it is itself visible in the history.

Hash basis (PI Round 4g §1)
---------------------------
A sha256 is only meaningful together with the bytes it was taken over. On
Windows with ``core.autocrlf=true`` a file can be LF in the git blob and CRLF in
the working tree, so a working-tree hash would pass here and fail on a Linux
checkout of the same commit. Locks therefore record ``blob_sha256`` and the guard
verifies ``git cat-file blob HEAD:<path>``. A working tree that no longer equals
that blob is *also* refused: it means either an uncommitted edit or an EOL
rewrite, and in both cases the file being read is not the file that was frozen.
"""

from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path
from typing import Iterable


class EmbargoViolation(RuntimeError):
    """Raised when a confirmatory read is attempted without the required state."""


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def blob_at_head(root: Path, rel: str) -> bytes | None:
    """Bytes of ``HEAD:<rel>``, or None when the path is not in the commit."""
    proc = subprocess.run(
        ["git", "cat-file", "blob", f"HEAD:{rel}"],
        cwd=root,
        capture_output=True,
        check=False,
    )
    return proc.stdout if proc.returncode == 0 else None


def _touches_test_novel(paths: Iterable[str | Path]) -> bool:
    return any("test_novel" in str(p).replace("\\", "/").lower() for p in paths)


def check_amendments(root: Path, lockfile: Path) -> dict[str, str]:
    """Return the verified hashes, or raise if any file has moved."""
    lock = json.loads((root / lockfile).read_text(encoding="utf-8"))
    default_basis = lock.get("hash_basis", "worktree")
    mismatches = []
    verified: dict[str, str] = {}
    for rel, spec in lock["files"].items():
        if isinstance(spec, dict):
            expected = spec.get("blob_sha256") or spec.get("sha256") or ""
            basis = spec.get("basis", default_basis)
        else:
            expected = spec
            basis = default_basis

        if basis == "git-blob":
            blob = blob_at_head(root, rel)
            if blob is None:
                mismatches.append(f"{rel}: not present in HEAD")
                continue
            actual = sha256_bytes(blob)
            if actual != expected:
                mismatches.append(f"{rel}: blob {actual} != {expected}")
                continue
            target = root / rel
            if not target.exists():
                mismatches.append(f"{rel}: missing from the working tree")
                continue
            on_disk = sha256_bytes(target.read_bytes())
            if on_disk != actual:
                mismatches.append(
                    f"{rel}: working tree {on_disk} != committed blob {actual} "
                    "(uncommitted edit or line-ending rewrite)"
                )
                continue
            verified[rel] = actual
        else:
            target = root / rel
            if not target.exists():
                mismatches.append(f"{rel}: missing")
                continue
            actual = sha256_file(target)
            if actual != expected:
                mismatches.append(f"{rel}: {actual} != {expected}")
                continue
            verified[rel] = actual
    if mismatches:
        raise EmbargoViolation(
            "pre-registration files changed since their hashes were recorded: "
            + "; ".join(mismatches)
        )
    return verified


def assert_allowed(
    paths: Iterable[str | Path],
    *,
    confirmatory: bool,
    root: Path,
    lockfile: Path = Path("state/amendments.lock.json"),
) -> None:
    """No-op for ordinary paths; hard failure for an unguarded confirmatory read."""
    if not _touches_test_novel(paths):
        return
    if not confirmatory:
        raise EmbargoViolation(
            "this path is under test_novel, which is embargoed; "
            "pass --confirmatory only after H2/P5 are pre-registered and the "
            "amendment hashes still match (Amendment 6 A6.5)"
        )
    check_amendments(root, lockfile)


__all__ = [
    "EmbargoViolation",
    "assert_allowed",
    "blob_at_head",
    "check_amendments",
    "sha256_bytes",
    "sha256_file",
]
