"""Contract tests for the D1/K1 measurement.

The important ones are (a) the sparse decoder agrees with the production dense
decoder, (b) the single-pass threshold sweep agrees with decoding once per cut,
and (c) the K1 verdict actually responds to a calibration-only difference.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from rff.analysis.scoring import (  # noqa: E402
    decode_pairs,
    increment_sweep,
    length_band,
)
from rff.models.dp import pair_map_to_pairs  # noqa: E402


def test_length_bands_match_the_p2_table():
    assert length_band(79) == "L<80"
    assert length_band(80) == "80-149"
    assert length_band(149) == "80-149"
    assert length_band(150) == "150-299"
    assert length_band(300) == "L>=300"


def test_sparse_decoder_matches_production_dense_decoder():
    rng = np.random.default_rng(0)
    for _ in range(20):
        n = 60
        prob = rng.random((n, n))
        prob = 0.5 * (prob + prob.T)
        np.fill_diagonal(prob, 0.0)
        dense = pair_map_to_pairs(prob, threshold=0.5)
        sparse = {
            (i, j): float(prob[i, j])
            for i in range(n)
            for j in range(i + 4, n)
            if prob[i, j] >= 0.5
        }
        assert decode_pairs(sparse, 0.5, n) == dense


def test_bare_arm_does_not_enforce_partners_or_nesting():
    scored = {(0, 20): 0.9, (1, 21): 0.8, (2, 19): 0.7, (0, 5): 0.6}
    assert decode_pairs(scored, 0.5, 40, partner=False) == [
        (0, 5),
        (0, 20),
        (1, 21),
        (2, 19),
    ]
    partner_only = decode_pairs(scored, 0.5, 40, partner=True, nesting=False)
    assert (0, 20) in partner_only and (0, 5) not in partner_only
    nested = decode_pairs(scored, 0.5, 40, partner=True, nesting=True)
    # 0.9 pair first; the 0.8 pair crosses it, the 0.7 pair nests inside it.
    assert nested == [(0, 20), (2, 19)]


def test_increment_sweep_matches_one_decode_per_cut():
    rng = np.random.default_rng(3)
    rows = []
    for _ in range(5):
        scored = {}
        for i in range(40):
            for j in range(i + 4, 40):
                p = float(rng.random())
                if p >= 0.2:
                    scored[(i, j)] = p
        rows.append({"length": 40, "_scored": scored})
    cuts = [0.2, 0.35, 0.5, 0.7]
    for flags in (
        dict(partner=False, nesting=False),
        dict(partner=True, nesting=False),
        dict(partner=True, nesting=True),
    ):
        swept = increment_sweep(rows, cuts, **flags)
        for k, row in enumerate(rows):
            for c, cut in enumerate(cuts):
                assert swept[k][c] == decode_pairs(row["_scored"], cut, 40, **flags)


TRUTH = [(0, 20), (1, 19), (2, 18)]
FAKE = [(40, 60), (41, 59), (42, 58)]


def _dump_pair(dirpath: Path, a_scores: dict, e_scores: dict, n_families: int = 3):
    dirpath.mkdir(parents=True, exist_ok=True)
    for name, scores in (("a", a_scores), ("e", e_scores)):
        with open(dirpath / f"{name}.jsonl", "w") as fh:
            for f in range(n_families):
                fh.write(
                    json.dumps(
                        {
                            "id": f"{name}_seq{f}",
                            "family": f"FAM{f}",
                            "source": "synthetic",
                            "length": 100,
                            "kept": [43, 57],
                            "truth": [list(p) for p in TRUTH],
                            "scored": {f"{i},{j}": v for (i, j), v in scores.items()},
                        }
                    )
                    + "\n"
                )


def _run(monkeypatch, module, argv):
    monkeypatch.setattr(sys, "argv", argv)
    return module.main()


def test_k1_fires_when_the_gap_is_purely_a_scale_shift(tmp_path, monkeypatch):
    import d1_scoring
    import k1_check

    a = {p: 0.45 for p in TRUTH} | {p: 0.10 for p in FAKE}
    e = {p: 0.90 for p in TRUTH} | {p: 0.30 for p in FAKE}
    dumps = tmp_path / "dump"
    _dump_pair(dumps, a, e)
    out = tmp_path / "d1.json"
    assert (
        _run(
            monkeypatch,
            d1_scoring,
            ["d1", "--dump-dir", str(dumps), "--out", str(out)],
        )
        == 0
    )
    summary = json.loads(out.read_text())
    band = summary["a"]["arms"]["partner"]["bands"]["80-299"]
    assert band["at_fixed"]["family_macro_f1"] == 0.0
    assert band["at_oracle"]["family_macro_f1"] == 1.0
    assert summary["e"]["arms"]["partner"]["bands"]["80-299"]["at_fixed"]["family_macro_f1"] == 1.0

    per_seq = out.with_suffix(out.suffix + ".per_sequence")
    result = tmp_path / "k1.json"
    assert (
        _run(
            monkeypatch,
            k1_check,
            [
                "k1",
                "--a",
                str(per_seq / "a.jsonl"),
                "--e",
                str(per_seq / "e.jsonl"),
                "--n-boot",
                "500",
                "--out",
                str(result),
            ],
        )
        == 0
    )
    verdict = json.loads(result.read_text())
    assert verdict["k1_fires"] is True
    assert abs(verdict["arms"]["partner"]["delta"]) < 1e-9


def test_k1_does_not_fire_when_e_is_genuinely_better(tmp_path, monkeypatch):
    import d1_scoring
    import k1_check

    a = {(0, 20): 0.60, (1, 19): 0.55, (2, 18): 0.50}
    a |= {(40, 60): 0.58, (41, 59): 0.57, (42, 58): 0.52}
    e = {p: v for p, v in zip(TRUTH, (0.90, 0.85, 0.80))}
    e |= {p: v for p, v in zip(FAKE, (0.20, 0.15, 0.10))}
    dumps = tmp_path / "dump"
    _dump_pair(dumps, a, e)
    out = tmp_path / "d1.json"
    _run(monkeypatch, d1_scoring, ["d1", "--dump-dir", str(dumps), "--out", str(out)])
    summary = json.loads(out.read_text())
    a_best = summary["a"]["arms"]["partner"]["bands"]["80-299"]["at_oracle"]["family_macro_f1"]
    e_best = summary["e"]["arms"]["partner"]["bands"]["80-299"]["at_oracle"]["family_macro_f1"]
    assert a_best < e_best - 0.01, (a_best, e_best)

    per_seq = out.with_suffix(out.suffix + ".per_sequence")
    result = tmp_path / "k1.json"
    _run(
        monkeypatch,
        k1_check,
        [
            "k1",
            "--a",
            str(per_seq / "a.jsonl"),
            "--e",
            str(per_seq / "e.jsonl"),
            "--n-boot",
            "500",
            "--out",
            str(result),
        ],
    )
    verdict = json.loads(result.read_text())
    assert verdict["k1_fires"] is False


def test_v3_label_policy_excludes_rows_and_counts_their_predictions(tmp_path, monkeypatch):
    import d1_scoring

    dumps = tmp_path / "dump"
    dumps.mkdir()
    rows = []
    for idx, truth_v3 in enumerate(
        (
            [[0, 10], [1, 9], [2, 8], [3, 7]],
            [[0, 12], [1, 11], [2, 10], [3, 9]],
            [],
            None,
        )
    ):
        rows.append(
            {
                "id": f"seq{idx}",
                "family": f"FAM{idx}",
                "source": "synthetic",
                "length": 100,
                "kept": [0, 100],
                "truth": [[0, 20], [1, 19], [2, 18]],
                "truth_v3": truth_v3,
                "scored": {"0,10": 0.9, "5,25": 0.8},
            }
        )
    for name in ("a", "e"):
        with open(dumps / f"{name}.jsonl", "w") as fh:
            for row in rows:
                fh.write(json.dumps(row) + "\n")

    out = tmp_path / "d1_v3.json"
    _run(
        monkeypatch,
        d1_scoring,
        [
            "d1",
            "--dump-dir",
            str(dumps),
            "--label-version",
            "v3",
            "--min-true-pairs",
            "4",
            "--out",
            str(out),
        ],
    )
    summary = json.loads(out.read_text())["a"]
    assert summary["n_sequences"] == 2
    policy = summary["label_policy"]
    assert policy["n_excluded_below_min"] == 1
    assert policy["n_excluded_label_parse_failure"] == 1
    assert policy["excluded_false_positive_diagnostic"]["n_rows"] == 2
    assert policy["excluded_false_positive_diagnostic"]["predicted_pairs_total"] > 0

