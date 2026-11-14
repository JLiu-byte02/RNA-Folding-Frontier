"""The public repository must never contain credentials.

The repo is public and now carries many generated artifacts, so a single pasted
token would be permanent. This test scans every tracked text file for private
key blocks and literal credential assignments, and rejects key-like filenames.
It runs in the fresh-clone CI job, i.e. against exactly what a reader can see.

Non-secret infrastructure identifiers (cluster hostname, campus VPN endpoints)
are documented separately in ``docs/public_repo_hygiene.md``; they are not
credentials and are intentionally kept for reproducibility.
"""

from __future__ import annotations

import re
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]

SKIP_SUFFIXES = {".gz", ".pt", ".pth", ".png", ".jpg", ".jpeg", ".pdf", ".zip", ".npz", ".npy"}

PATTERNS = [
    ("private key block", re.compile(r"BEGIN [A-Z ]*PRIVATE KEY")),
    ("literal SSHPASS assignment", re.compile(r"\bSSHPASS\s*=\s*[^\s\"']+")),
    ("literal password assignment", re.compile(r"\bpassword\s*=\s*[\"'][^\"']+[\"']", re.I)),
    ("literal api key", re.compile(r"\bapi[_-]?key\s*=\s*[\"'][A-Za-z0-9_\-]{16,}[\"']", re.I)),
    ("literal token", re.compile(r"\btoken\s*=\s*[\"'][A-Za-z0-9_\-\.]{20,}[\"']", re.I)),
]

KEY_FILE = re.compile(r"(^|/)(id_(rsa|dsa|ecdsa|ed25519)|\.env|secrets?\.(json|ya?ml|txt))$")


def tracked_files() -> list[Path]:
    out = subprocess.run(
        ["git", "ls-files", "-z"], cwd=ROOT, capture_output=True, check=True
    ).stdout.decode()
    return [ROOT / p for p in out.split("\0") if p]


def test_no_key_like_files_are_tracked():
    bad = [str(p.relative_to(ROOT)) for p in tracked_files() if KEY_FILE.search(str(p.relative_to(ROOT)).replace("\\", "/"))]
    assert not bad, f"key-like files are tracked: {bad}"


def test_no_credential_literals_in_tracked_files():
    findings = []
    for path in tracked_files():
        if path.suffix.lower() in SKIP_SUFFIXES or not path.is_file():
            continue
        try:
            text = path.read_text(encoding="utf-8", errors="ignore")
        except OSError:
            continue
        for label, pattern in PATTERNS:
            match = pattern.search(text)
            if match:
                line = text[: match.start()].count("\n") + 1
                findings.append(f"{path.relative_to(ROOT)}:{line} ({label})")
    assert not findings, "credential-like literals found: " + "; ".join(findings)
