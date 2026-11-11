"""Access guard for the confirmatory (test_novel) evaluation set.

Amendment 6 §A6.5. The confirmatory scripts must refuse a path under
``test_novel`` unless two independent conditions hold:

1. the caller passed ``--confirmatory``;
2. every pre-registration file listed in ``state/amendments.lock.json`` still
   hashes to the recorded value.

The second condition is what makes the flag meaningful: a flag can be typed at
any time, but it cannot make an edited amendment match its old hash. The lockfile
is committed, so editing it is itself visible in the history.
"""

from __future__ import annotations

import hashlib
import json
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


def _touches_test_novel(paths: Iterable[str | Path]) -> bool:
    return any("test_novel" in str(p).replace("\\", "/").lower() for p in paths)


def check_amendments(root: Path, lockfile: Path) -> dict[str, str]:
    """Return the verified hashes, or raise if any file has moved."""
    lock = json.loads((root / lockfile).read_text(encoding="utf-8"))
    mismatches = []
    for rel, expected in lock["files"].items():
        target = root / rel
        if not target.exists():
            mismatches.append(f"{rel}: missing")
            continue
        actual = sha256_file(target)
        if actual != expected:
            mismatches.append(f"{rel}: {actual} != {expected}")
    if mismatches:
        raise EmbargoViolation(
            "pre-registration files changed since their hashes were recorded: "
            + "; ".join(mismatches)
        )
    return dict(lock["files"])


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


__all__ = ["EmbargoViolation", "assert_allowed", "check_amendments", "sha256_file"]
