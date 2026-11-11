"""Contract tests for the sequence-level positive-control split."""

from __future__ import annotations

import json
import random
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from build_seqlvl_split import (  # noqa: E402
    cluster_clan,
    deal_clusters,
    identity,
    kmer_hashes,
    minhash,
)


def _aligner():
    from Bio import Align

    al = Align.PairwiseAligner()
    al.mode = "global"
    al.match_score = 2.0
    al.mismatch_score = -1.0
    al.open_gap_score = -6.0
    al.extend_gap_score = -1.0
    return al


def _mutate(seq: str, n: int, rng: random.Random) -> str:
    out = list(seq)
    for pos in rng.sample(range(len(seq)), n):
        out[pos] = rng.choice([b for b in "ACGU" if b != out[pos]])
    return "".join(out)


def _rand_seq(n: int, rng: random.Random) -> str:
    return "".join(rng.choice("ACGU") for _ in range(n))


def test_identity_of_identical_and_mutated_sequences():
    rng = random.Random(0)
    seq = _rand_seq(120, rng)
    al = _aligner()
    assert identity(al, seq, seq, 0.8) == 1.0
    near = _mutate(seq, 2, rng)
    assert identity(al, seq, near, 0.8) >= 0.95


def test_far_sequences_are_rejected_by_the_kmer_filter():
    rng = random.Random(1)
    a, b = _rand_seq(150, rng), _rand_seq(150, rng)
    assert identity(_aligner(), a, b, 0.8) == 0.0
    assert len(kmer_hashes(a) & kmer_hashes(b)) < 0.2 * len(kmer_hashes(a))


def test_clustering_keeps_near_duplicates_together():
    rng = random.Random(2)
    base = _rand_seq(140, rng)
    recs = [
        {"id": "a1", "seq": base},
        {"id": "a2", "seq": _mutate(base, 3, rng)},
        {"id": "a3", "seq": _mutate(base, 6, rng)},
        {"id": "b1", "seq": _rand_seq(140, rng)},
        {"id": "b2", "seq": _rand_seq(140, rng)},
    ]
    cache: dict[str, list[int]] = {}
    clusters = cluster_clan(recs, _aligner(), 0.8, cache)
    by_id = {i: r["id"] for i, r in enumerate(recs)}
    groups = sorted(sorted(by_id[i] for i in g) for g in clusters)
    assert ["a1", "a2", "a3"] in groups
    assert len(groups) == 3


def test_deal_clusters_never_splits_a_cluster_and_fills_required_splits():
    counts = {"train": 0, "val": 0, "test": 0}
    clusters = [[0, 1], [2], [3], [4, 5, 6]]
    rng = np.random.default_rng(0)
    assigned = deal_clusters(clusters, counts, {"train": 0.75, "val": 0.1, "test": 0.15}, rng)
    assert set(assigned) == {0, 1, 2, 3}
    assert counts["train"] > 0 and counts["test"] > 0
    assert counts["train"] + counts["val"] + counts["test"] == 7


def test_end_to_end_split_is_cluster_disjoint(tmp_path, monkeypatch):
    import build_seqlvl_split

    rng = random.Random(3)
    rows = []
    for clan in ("CL1", "CL2"):
        base = _rand_seq(150, rng)
        for k in range(3):
            rows.append(
                {
                    "id": f"{clan}_dup{k}",
                    "seq": _mutate(base, k, rng),
                    "pairs": [[0, 10]],
                    "family": f"{clan}_FAM",
                    "clan": clan,
                    "source": "synthetic",
                }
            )
        for k in range(3):
            rows.append(
                {
                    "id": f"{clan}_uniq{k}",
                    "seq": _rand_seq(150, rng),
                    "pairs": [[1, 11]],
                    "family": f"{clan}_FAM",
                    "clan": clan,
                    "source": "synthetic",
                }
            )
    src = tmp_path / "train.jsonl"
    with open(src, "w") as fh:
        for row in rows:
            fh.write(json.dumps(row) + "\n")
    out = tmp_path / "split"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "build_seqlvl_split",
            "--train",
            str(src),
            "--out",
            str(out),
            "--identity",
            "0.8",
            "--val-frac",
            "0.15",
            "--test-frac",
            "0.25",
        ],
    )
    assert build_seqlvl_split.main() == 0
    manifest = json.loads((out / "manifest.json").read_text())
    assert manifest["max_train_test_identity"] < 0.8
    assert sum(manifest["counts"].values()) == len(rows)
    assert manifest["counts"]["train"] > 0 and manifest["counts"]["test"] > 0
    # every clan with >=2 clusters must appear on both sides
    for clan, info in manifest["per_clan"].items():
        if info["n_clusters"] >= 2:
            assert info["splits"].get("train", 0) > 0, clan
            assert info["splits"].get("test", 0) > 0, clan
    # no sequence id appears twice across splits
    ids = []
    for split in ("train", "val", "test"):
        f = out / f"{split}.jsonl"
        if f.exists():
            ids += [json.loads(line)["id"] for line in open(f)]
    assert len(ids) == len(set(ids)) == len(rows)
