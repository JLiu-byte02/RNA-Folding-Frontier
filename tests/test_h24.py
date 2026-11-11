"""Contract tests for the H2.4 decoder comparison.

The point of H2.4 is that consistency can remove a long-range pair that the
greedy decoder accepts because it outranks the true local pairs. The test builds
exactly that situation on a toy sequence and checks the script's accounting.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _seq_with_pairs(length: int, pairs, extra=()) -> str:
    seq = ["A"] * length
    for i, j in list(pairs) + list(extra):
        seq[i], seq[j] = "G", "C"
    return "".join(seq)


def _run(tmp_path: Path, seq: str, truth, scored: dict[str, float], name="toy") -> dict:
    split = tmp_path / "split.jsonl"
    split.write_text(json.dumps({"id": "s1", "seq": seq, "pairs": []}) + "\n")
    dump = tmp_path / "dump.jsonl"
    dump.write_text(
        json.dumps(
            {
                "id": "s1",
                "family": "f1",
                "length": len(seq),
                "kept": [0, 1],
                "truth": truth,
                "scored": scored,
            }
        )
        + "\n"
    )
    out = tmp_path / "out.json"
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts" / "h24_decoder_compare.py"),
            "--dump", str(dump),
            "--split", str(split),
            "--name", name,
            "--out", str(out),
            "--n-boot", "100",
        ],
        check=True,
        capture_output=True,
    )
    return json.loads(out.read_text())


def test_dp_drops_long_range_pair_that_outranks_true_pairs(tmp_path):
    # The spurious pair crosses both true pairs, so the three are mutually
    # exclusive. Greedy takes it on its higher probability (2 true pairs lost);
    # the DP compares total log-odds (2 x 1.10 > 1.39) and keeps the truth.
    truth = [[5, 12], [20, 160]]
    spurious = (10, 140)
    seq = _seq_with_pairs(200, [tuple(p) for p in truth], [spurious])
    scored = {f"{i},{j}": 0.75 for i, j in truth}
    scored[f"{spurious[0]},{spurious[1]}"] = 0.8

    report = _run(tmp_path, seq, truth, scored)
    fm = report["family_macro_per_sequence"]
    assert fm["greedy"]["long_range_fp"] == 1
    assert fm["greedy"]["tp"] == 0
    assert fm["dp"]["long_range_fp"] == 0
    assert fm["dp"]["tp"] == 2
    assert report["H2_4"]["mean_d_long_range_fp"] == -1.0
    assert report["H2_4"]["mean_d_tp"] == 2.0
    # one family only: the contrast cannot be bootstrapped, and the script says so
    assert report["H2_4"]["contrast_drop_lrfp_minus_drop_tp"]["contrast"] is None


def test_identical_decoders_give_zero_contrast(tmp_path):
    truth = [[2, 20], [30, 48]]
    seq = _seq_with_pairs(80, [tuple(p) for p in truth])
    scored = {f"{i},{j}": 0.9 for i, j in truth}
    report = _run(tmp_path, seq, truth, scored)
    assert report["H2_4"]["mean_d_tp"] == 0.0
    assert report["H2_4"]["mean_d_long_range_fp"] == 0.0
    # a single family cannot be bootstrapped; the script must say so, not crash
    assert report["H2_4"]["contrast_drop_lrfp_minus_drop_tp"]["contrast"] is None
    assert report["H2_4"]["passes"] is False


def test_dp_output_is_a_valid_structure(tmp_path):
    truth = [[2, 20]]
    seq = _seq_with_pairs(60, [(2, 20)])
    scored = {"2,20": 0.9, "4,30": 0.95, "6,40": 0.85}
    report = _run(tmp_path, seq, truth, scored)
    assert report["n_sequences"] == 1
    assert report["dp_config"]["hairpin_penalty"] == 0.0
    assert report["family_macro_per_sequence"]["dp"]["tp"] >= 1
