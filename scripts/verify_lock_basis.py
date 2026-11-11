#!/usr/bin/env python3
"""Which bytes do the locked hashes actually cover?

Amendment 6 §A6.5 pins a sha256 per pre-registration file in
``state/amendments.lock.json``. On Windows with ``core.autocrlf=true`` a file can
be CRLF in the working tree and LF in the git blob, so a hash computed on the
wrong one would look fine here and fail on a Linux checkout. This script answers
the question three ways for every locked file:

1. the working-tree bytes,
2. the blob at HEAD,
3. every blob that ever existed for that path, with the commits where it was
   the version on disk.

The report says whether all three agree, which commit first carried the locked
bytes, and whether the working tree is dirty. A hash that only ever matched a
working-tree file is flagged ``worktree-only`` — that is the failure mode the PI
asked to rule out.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def git(*args: str) -> bytes:
    proc = subprocess.run(
        ["git", *args], cwd=ROOT, capture_output=True, check=False
    )
    if proc.returncode != 0:
        raise RuntimeError(
            f"git {' '.join(args)} failed: {proc.stderr.decode(errors='replace').strip()}"
        )
    return proc.stdout


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def blob_at(commit: str, path: str) -> bytes | None:
    proc = subprocess.run(
        ["git", "show", f"{commit}:{path}"], cwd=ROOT, capture_output=True, check=False
    )
    return proc.stdout if proc.returncode == 0 else None


def commits_touching(path: str) -> list[str]:
    out = git("log", "--format=%H", "--", path).decode().split()
    return out


def check(path: str, locked: str) -> dict:
    entry: dict = {"path": path, "locked": locked}
    worktree = ROOT / path
    wt_hash = sha256_bytes(worktree.read_bytes()) if worktree.exists() else None
    entry["worktree_sha256"] = wt_hash
    entry["worktree_eol_crlf"] = (
        b"\r\n" in worktree.read_bytes() if worktree.exists() else None
    )

    head_blob = blob_at("HEAD", path)
    entry["head_blob_sha256"] = sha256_bytes(head_blob) if head_blob is not None else None

    matches = []
    history = []
    for commit in commits_touching(path):
        blob = blob_at(commit, path)
        if blob is None:
            continue
        digest = sha256_bytes(blob)
        history.append({"commit": commit[:12], "sha256": digest, "crlf": b"\r\n" in blob})
        if digest == locked:
            matches.append(commit[:12])

    entry["matching_commits"] = matches
    entry["history"] = history
    entry["blob_matches"] = bool(matches)
    entry["worktree_matches"] = wt_hash == locked
    entry["worktree_equals_head"] = wt_hash is not None and wt_hash == entry["head_blob_sha256"]
    if entry["blob_matches"] and entry["worktree_matches"]:
        entry["verdict"] = "both"
    elif entry["blob_matches"]:
        entry["verdict"] = "blob-only"      # working tree is dirty or mis-hashed
    elif entry["worktree_matches"]:
        entry["verdict"] = "worktree-only"  # the failure mode the PI asked about
    else:
        entry["verdict"] = "neither"
    return entry


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--lockfile", default="state/amendments.lock.json")
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    lock = json.loads((ROOT / args.lockfile).read_text(encoding="utf-8"))
    results = []
    for rel, value in lock["files"].items():
        if isinstance(value, dict):
            value = value.get("blob_sha256") or value.get("sha256") or ""
        results.append(check(rel, value))

    ok = all(r["verdict"] == "both" for r in results)
    report = {
        "hash_basis_question": "do working-tree bytes and git blob bytes give the same sha256?",
        "lockfile": args.lockfile,
        "all_agree": ok,
        "files": results,
    }
    text = json.dumps(report, indent=2) + "\n"
    if args.out:
        out = Path(args.out)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(text, encoding="utf-8")
        print(json.dumps({"out": str(out), "all_agree": ok}))
    else:
        print(text)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
