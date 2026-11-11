"""Dry run of the D1 → K1 and D2 diagnostic chain on synthetic dumps.

The confirmatory path cannot be allowed to fail on a schema bug the hour the
real dumps land. This builds dumps with exactly the documented fields
(`id`, `family`, `length`, `kept`, `truth`, `scored`), runs the three analysis
scripts as subprocesses, and checks that each output has the keys the next step
and the paper table consume.

It is deliberately not a scientific test: the probabilities are noise, so no
sign of any effect is asserted. Only wiring, schema and non-crashing behaviour.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def _make_sequence(n: int, pairs):
    seq = ["A"] * n
    for i, j in pairs:
        seq[i], seq[j] = "G", "C"
    return "".join(seq)


def _build(tmp_path: Path, arm: str, seed: int, n_rows: int = 12):
    import random

    rng = random.Random(seed)
    dump = tmp_path / f"{arm}_test_withheld.jsonl"
    split_rows = []
    with open(dump, "w") as fh:
        for k in range(n_rows):
            n = 90 + 10 * k
            i, j = 10 + k, n - 15
            truth = [[i, j]]
            seq = _make_sequence(n, [tuple(p) for p in truth])
            scored = {
                f"{i},{j}": round(rng.uniform(0.45, 0.95), 3),
                f"20,{n - 30}": round(rng.uniform(0.2, 0.8), 3),
                f"5,{n - 8}": round(rng.uniform(0.2, 0.8), 3),
            }
            fh.write(
                json.dumps(
                    {
                        "id": f"s{k}",
                        "family": f"f{k % 3}",
                        "length": n,
                        "kept": [0, n],
                        "truth": truth,
                        "scored": scored,
                    }
                )
                + "\n"
            )
            split_rows.append({"id": f"s{k}", "seq": seq, "family": f"f{k % 3}", "pairs": truth})
    split = tmp_path / f"split_{arm}.jsonl"
    split.write_text("".join(json.dumps(r) + "\n" for r in split_rows))
    return dump, split


def _run(*args):
    proc = subprocess.run(
        [sys.executable, *args], cwd=ROOT, capture_output=True, text=True
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
    return proc


def test_d1_writes_summary_and_per_sequence_for_k1(tmp_path):
    dump_dir = tmp_path / "dumps"
    dump_dir.mkdir()
    for arm in ("a", "e"):
        dump, _ = _build(dump_dir, arm, seed=3)
        assert dump.exists()

    out = tmp_path / "d1.json"
    _run(
        str(ROOT / "scripts" / "d1_scoring.py"),
        "--dump-dir", str(dump_dir),
        "--out", str(out),
    )
    summary = json.loads(out.read_text())
    assert set(summary) >= {"a_test_withheld", "e_test_withheld"}
    entry = summary["a_test_withheld"]
    assert "arms" in entry
    for arm in ("bare", "partner", "nested"):
        assert arm in entry["arms"]
        assert "80-299" in entry["arms"][arm]["bands"]
        band = entry["arms"][arm]["bands"]["80-299"]
        if band.get("n"):
            assert "at_fixed" in band and "at_oracle" in band

    per_seq = out.with_suffix(out.suffix + ".per_sequence")
    rows = [json.loads(l) for l in (per_seq / "a_test_withheld.jsonl").read_text().splitlines()]
    assert rows and {"id", "family", "length", "f1_nested_oracle_80_299"} <= set(rows[0])


def test_k1_consumes_d1_per_sequence_files(tmp_path):
    dump_dir = tmp_path / "dumps"
    dump_dir.mkdir()
    for arm in ("a", "e"):
        _build(dump_dir, arm, seed=11)
    out = tmp_path / "d1.json"
    _run(
        str(ROOT / "scripts" / "d1_scoring.py"),
        "--dump-dir", str(dump_dir), "--out", str(out),
    )
    per_seq = out.with_suffix(out.suffix + ".per_sequence")
    a_file = per_seq / "a_test_withheld.jsonl"
    e_file = per_seq / "e_test_withheld.jsonl"

    k1_out = tmp_path / "k1.json"
    _run(
        str(ROOT / "scripts" / "k1_check.py"),
        "--a", str(a_file), str(a_file),
        "--e", str(e_file), str(e_file),
        "--n-boot", "200",
        "--out", str(k1_out),
    )
    report = json.loads(k1_out.read_text())
    text = json.dumps(report)
    assert "e_minus_a" in text or "delta" in text


def test_d2_needs_only_dump_and_split(tmp_path):
    dump_dir = tmp_path / "dumps"
    dump_dir.mkdir()
    _, split = _build(dump_dir, "a", seed=5)
    out = tmp_path / "d2.json"
    _run(
        str(ROOT / "scripts" / "d2_structure_errors.py"),
        "--dump-dir", str(dump_dir),
        "--names", "a_test_withheld",
        "--split", str(split),
        "--n-boot", "200",
        "--null-reps", "2",
        "--out", str(out),
    )
    report = json.loads(out.read_text())
    assert report
    text = json.dumps(report)
    assert "both_kept" in text
