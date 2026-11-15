"""Exploratory FP anatomy must describe the deployed head, not a map heuristic."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
COHORT = ROOT / "experiments" / "002_e1prime_runs" / "cohort"
SCRIPTS = ROOT / "scripts"


def _json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def test_fp_artifacts_match_frozen_head17() -> None:
    for split in ("test_withheld", "id_test"):
        paired = _json(COHORT / f"paired_{split}.json")
        anatomy = _json(COHORT / f"fp_anatomy_{split}.json")
        oracle = _json(COHORT / f"shift_oracle_{split}.json")
        head = paired["methods"]["head17"]

        assert anatomy["n_sequences"] == paired["n_sequences"]
        assert oracle["n_families"] == paired["n_families"]
        assert anatomy["methods"]["head"]["n_fp"] == head["fp"] == oracle["n_fp"]
        assert round(oracle["baseline"]["macro_f1"], 4) == round(head["family_macro_f1"], 4)
        assert round(oracle["baseline"]["micro_f1"], 4) == round(head["micro_f1"], 4)

        for method in anatomy["methods"].values():
            assert abs(sum(method["frac"].values()) - 1.0) <= 0.0005
            assert method["near_miss_frac"] == method["frac"]["1-2"]
            assert 0 <= method["n_sequences_with_fp"] <= anatomy["n_sequences"]

        assert 0 <= oracle["n_shiftable"] <= oracle["n_fp"]
        assert oracle["oracle_shift"]["macro_f1"] >= oracle["baseline"]["macro_f1"]


def test_fp_scripts_reject_threshold_map_inputs() -> None:
    for name in ("fp_anatomy.py", "shift_oracle.py"):
        source = (SCRIPTS / name).read_text(encoding="utf-8")
        assert '"--head-preds"' in source
        assert '"--threshold"' not in source
        assert '"--dump"' not in source
