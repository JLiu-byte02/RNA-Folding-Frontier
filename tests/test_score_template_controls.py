"""Contract test for the fast template-control scorer."""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_template_controls import score_dump  # noqa: E402


def test_score_dump_uses_v3_primary_rows(tmp_path):
    seq = "GGGAAACCCGGGAAACCC"
    rows = [
        {
            "id": "s1",
            "family": "F1",
            "length": len(seq),
            "truth": [[0, 7]],
            "truth_v3": [[0, 7], [1, 6], [8, 15], [9, 14]],
            "scored": {"0,7": 0.9, "1,6": 0.9},
        },
        {
            "id": "s2",
            "family": "F2",
            "length": len(seq),
            "truth": [],
            "truth_v3": [[2, 5]],
            "scored": {},
        },
    ]
    path = tmp_path / "dump.jsonl"
    path.write_text("".join(json.dumps(r) + "\n" for r in rows))
    result = score_dump(path, 0.5, 4)
    assert result["n_scored"] == 1
    assert result["n_excluded"] == 1
    assert result["family_macro_f1"] is not None
