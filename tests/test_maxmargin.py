"""Contract tests for the structured max-margin step.

The claim these tests defend is narrow but essential: the subgradient moves the
score matrix in the direction that makes the true structure decode next time.
A sign error here would train a head to predict the complement of the truth and
would look like "the method does not work".
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rff.decoders.nn_dp import NEG, NNDecoder, PairScoreModel, mask_allowed  # noqa: E402
from rff.training.maxmargin import (  # noqa: E402
    indicator,
    loss_augmented,
    margin_violation,
    structured_step,
    train_epoch,
)


def decode_with(scores):
    def decode(seq, S):
        return NNDecoder(PairScoreModel(pair_score=S)).decode(seq, S)

    return decode


def test_indicator_and_violation_zero_when_equal():
    scores = np.zeros((6, 6))
    pairs = [(0, 5)]
    assert indicator(pairs, scores.shape).sum() == 1
    assert margin_violation(scores, pairs, pairs) == 0.0


def test_loss_augmented_skips_truth_and_disallowed_cells():
    scores = np.zeros((6, 6))
    mask = mask_allowed("GGGAAA", 3)
    truth = [(0, 5)]
    aug = loss_augmented(scores, truth, mask, margin=2.0)
    assert aug[0, 5] == 0.0
    for i in range(6):
        for j in range(i + 1, 6):
            if (i, j) == (0, 5) or not mask[i, j]:
                continue
            assert aug[i, j] == pytest.approx(2.0)


def test_one_step_increases_truth_minus_prediction_score():
    seq = "GGGAAACCC"
    mask = mask_allowed(seq)
    truth = [(0, 8), (1, 7)]
    theta = np.full((len(seq), len(seq)), NEG)
    theta[mask] = 0.5  # ties everywhere: the decoder picks something arbitrary
    before = margin_violation(theta, decode_with(None)(seq, theta)[0], truth)
    info = structured_step(theta, seq, truth, mask, decode_with(None), margin=3.0, lr=0.5)
    after = margin_violation(theta, info["predicted"], truth)
    assert after < before or info["loss"] == 0
    assert theta[0, 8] > 0.5 or (0, 8) in info["predicted"]


def test_training_reduces_loss_on_a_learnable_toy():
    # one sequence, one consistent label: with tied scores the DP starts at the
    # empty structure, so the first step must reward exactly the true pairs
    seqs = [("GGGAAACCC", [(1, 7), (2, 6)])]
    theta = np.full((9, 9), NEG)
    theta[mask_allowed(seqs[0][0])] = 0.0
    losses = []
    for _ in range(12):
        info = train_epoch(theta, seqs, decode_with(None), lambda s: mask_allowed(s), margin=2.0, lr=0.3)
        losses.append(info["mean_loss"])
    assert losses[-1] < losses[0]
    assert losses[-1] == 0.0


def test_theta_stays_disallowed_outside_the_mask():
    seq = "GGGAAACCC"
    mask = mask_allowed(seq)
    theta = np.full((9, 9), NEG)
    theta[mask] = 0.0
    structured_step(theta, seq, [(0, 8)], mask, decode_with(None), margin=1.0, lr=0.5)
    assert np.all(theta[~mask] <= NEG / 2)
