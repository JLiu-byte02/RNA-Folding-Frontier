"""Contract tests for the P2-block geometry and full-length donors."""

from __future__ import annotations

import random
import sys
from collections import Counter
from pathlib import Path

import torch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from build_block_cache import (  # noqa: E402
    assemble_sequence,
    block_context,
    block_spans,
)
from build_block_donors import build_donors  # noqa: E402


def test_block_spans_cover_each_position_once():
    spans = block_spans(300, 16)
    covered = [p for s, e in spans for p in range(s, e)]
    assert covered == list(range(300))
    assert len(spans) == 19
    assert block_spans(16, 16) == [(0, 16)]
    assert block_spans(1, 16) == [(0, 1)]


def test_block_context_keeps_native_flank_only():
    native = "A" * 200
    donor = "N" * 200
    comp = block_context(native, donor, 48, 64, 32)
    assert comp[:16] == "N" * 16
    assert comp[16:96] == "A" * 80  # 48-32 .. 64+32
    assert comp[96:] == "N" * 104


def test_block_context_clips_at_sequence_ends():
    native = "A" * 100
    donor = "N" * 100
    first = block_context(native, donor, 0, 16, 32)
    assert first[:48] == "A" * 48 and first[48:] == "N" * 52
    last = block_context(native, donor, 96, 100, 32)
    assert last[:64] == "N" * 64 and last[64:] == "A" * 36


def test_block_context_rejects_length_mismatch():
    try:
        block_context("ACGU", "ACG", 0, 2, 1)
    except ValueError:
        return
    raise AssertionError("expected ValueError")


class _SelfCheckBackbone:
    """Fake backbone: column 0 = position, column 1 = base identity vs native."""

    def __init__(self, native: str) -> None:
        self.native = native

    def embed_strings(self, strings, batch_size=64):
        out = []
        for text in strings:
            idx = torch.arange(len(self.native), dtype=torch.float32)
            same = torch.tensor(
                [float(a == b) for a, b in zip(text, self.native)]
            )
            out.append(torch.stack([idx, same], dim=1))
        return out


def test_assembly_takes_each_position_from_its_own_native_pass():
    native = ("ACGU" * 40)[:150]
    donor = "N" * 150
    emb = assemble_sequence(_SelfCheckBackbone(native), native, donor, 16, 32)
    assert emb.shape == (150, 2)
    assert emb[:, 0].tolist() == list(range(150))
    assert emb[:, 1].min().item() == 1.0, emb[:, 1].tolist()


def test_donor_e_is_all_n_and_full_length():
    recs = [{"id": "a", "seq": "ACGUGCAUGCUGCAUGCAUG", "clan": "X"}]
    donors = build_donors("e", recs, rng=random.Random(0))
    assert donors["a"] == "N" * len(recs[0]["seq"])


def test_donor_b_preserves_dinucleotide_composition():
    rng = random.Random(0)
    seqs = ["".join(rng.choice("ACGU") for _ in range(120)) for _ in range(5)]
    recs = [{"id": str(n), "seq": s, "clan": "X"} for n, s in enumerate(seqs)]
    donors = build_donors("b", recs, rng=random.Random(1))
    for rec in recs:
        d = donors[rec["id"]]
        assert len(d) == len(rec["seq"])
        assert Counter(zip(rec["seq"], rec["seq"][1:])) == Counter(zip(d, d[1:]))


def test_donor_d_uses_a_different_clan():
    recs = [
        {"id": "t", "seq": "A" * 50, "clan": "X"},
        {"id": "o", "seq": "G" * 80, "clan": "Y"},
    ]
    donors = build_donors("d", recs, rng=random.Random(0))
    assert donors["t"] == "G" * 50


def test_donor_c_takes_aligned_bases_and_keeps_native_at_gaps():
    from build_conditions import positional_map

    target = "ACGUGCAUGCGGG"
    donor = "ACGUGCAAGC"
    recs = [
        {"id": "t", "seq": target, "clan": "X"},
        {"id": "h", "seq": donor, "clan": "X"},
    ]
    donors = build_donors(
        "c", recs[:1], rng=random.Random(0), donor_pool=recs,
        homolog_map={"t": {"homolog": "h", "identity": 0.9}},
    )
    out = donors["t"]
    mapping = positional_map(target, donor)
    assert len(out) == len(target)
    assert any(j == -1 for j in mapping), "test needs a gapped column"
    for pos, j in enumerate(mapping):
        expected = target[pos] if j < 0 else donor[j]
        assert out[pos] == expected, (pos, j)
    assert sum(a != b for a, b in zip(out, target)) > 0
    assert any(j < 0 and out[pos] == target[pos] for pos, j in enumerate(mapping))
