"""Contract tests for baseline running/scoring.

The parsers are the risky part: a dot-bracket assigned to the wrong sequence
would silently corrupt every baseline number. These tests cover the three tool
output shapes and, end to end, check that a perfect predictor scores 1.0 and a
predictor missing whole families is not credited for them.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from rff.analysis.dotbracket import (  # noqa: E402
    ct_to_pairs,
    extract_predictions,
    parse_dotbracket,
    parse_ct,
)


def test_parse_dotbracket_nested():
    assert parse_dotbracket("((..))") == {(0, 5), (1, 4)}
    assert parse_dotbracket("....") == set()


def test_parse_dotbracket_pseudoknot_brackets():
    pairs = parse_dotbracket("([)]")
    assert pairs == {(0, 2), (1, 3)}


def test_parse_dotbracket_rejects_mismatches():
    with pytest.raises(ValueError):
        parse_dotbracket("((..)")
    with pytest.raises(ValueError):
        parse_dotbracket("..))")
    with pytest.raises(ValueError):
        parse_dotbracket("((..))", length=5)


def test_ct_roundtrip():
    ct = "\n".join(
        [
            "6",
            "1 G 0 2 5 0",
            "2 C 1 3 4 0",
            "3 A 2 0 0 0",
            "4 U 3 0 0 0",
            "5 G 4 1 2 0",
            "6 C 5 0 0 0",
        ]
    )
    assert parse_ct(ct) == {0: 4, 1: 3, 4: 1}
    assert ct_to_pairs(ct) == {(0, 4), (1, 3)}


def test_extract_rnafold_layout():
    seqs = {"a": "GGGAAACCC"}
    text = "GGGAAACCC\n(((...)))\n"
    assert extract_predictions(text, seqs) == {"a": "(((...)))"}


def test_extract_id_tab_db_and_id_space_db():
    seqs = {"x": "GGGAAACCC"}
    assert extract_predictions("x\t(((...)))\n", seqs) == {"x": "(((...)))"}
    assert extract_predictions("x (((...)))  -3.2\n", seqs) == {"x": "(((...)))"}


def test_extract_does_not_guess_unknown_structures():
    seqs = {"x": "GGGAAACCC", "y": "GGGAAACCC"}
    # identical sequences: matching by sequence would be ambiguous, so the
    # unlabelled RNAfold layout must not assign either id
    assert extract_predictions("GGGAAACCC\n(((...)))\n", seqs) == {}
    # a labelled line still works
    assert extract_predictions("y\t(((...)))\n", seqs) == {"y": "(((...)))"}


def _write_split(path: Path, rows) -> None:
    path.write_text("".join(json.dumps(r) + "\n" for r in rows))


def _rows():
    seq = "GGGAAACCC"
    return [
        {"id": "a1", "seq": seq, "family": "A", "pairs": [[0, 8], [1, 7]]},
        {"id": "a2", "seq": seq, "family": "A", "pairs": [[0, 8]]},
        {"id": "b1", "seq": seq, "family": "B", "pairs": [[0, 8]]},
        {"id": "b2", "seq": seq, "family": "B", "pairs": []},
    ]


def _run(tmp_path, extra_args):
    import subprocess

    split = tmp_path / "split.jsonl"
    _write_split(split, _rows())
    out = tmp_path / "out.json"
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts" / "run_baselines.py"),
            "--split", str(split),
            "--out", str(out),
            "--tool", "toy",
            *extra_args,
        ],
        check=True,
        capture_output=True,
    )
    return json.loads(out.read_text())


def test_perfect_fake_tool_scores_one(tmp_path):
    fake = tmp_path / "fake.py"
    fake.write_text(
        "import sys, json\n"
        "rows = [json.loads(l) for l in open(sys.argv[1]) if l.strip()]\n"
        "for r in rows:\n"
        "    db = ['.'] * len(r['seq'])\n"
        "    for i, j in r['pairs']:\n"
        "        db[i], db[j] = '(', ')'\n"
        "    print(r['id'] + '\\t' + ''.join(db))\n"
    )
    # the fake tool writes a split file of its own; feeding it the fasta is not
    # the point here, so the command ignores {fasta} on purpose
    split_as_jsonl = tmp_path / "rows.jsonl"
    _write_split(split_as_jsonl, _rows())
    report = _run(
        tmp_path,
        ["--cmd", f'"{sys.executable}" "{fake}" "{split_as_jsonl}"'],
    )
    assert report["n_missing"] == 0
    assert report["family_macro_f1"]["mean"] == pytest.approx(1.0)
    assert report["micro"]["f1"] == pytest.approx(1.0)


def test_missing_family_is_not_credited(tmp_path):
    preds = tmp_path / "preds.txt"
    preds.write_text("a1\t((.....))\na2\t(.......)\n")
    report = _run(tmp_path, ["--predictions", str(preds)])
    assert report["n_missing"] == 2
    # family A perfect, family B absent -> macro over scored families only
    assert report["family_macro_f1"]["mean"] == pytest.approx(1.0)
    assert report["family_macro_f1"]["n_families"] == 1
