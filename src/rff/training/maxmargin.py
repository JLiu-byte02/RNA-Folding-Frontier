"""Structured max-margin training for a pair-score head (MXfold2-style).

The head produces one score per candidate pair. Decoding is the structured
prediction, so the loss is not differentiable through the argmax; the standard
remedy is a structured SVM step:

    y_hat = argmax_y  score_theta(y) + loss(y, y_true)      (loss-augmented)
    grad  = pair_indicator(y_hat) - pair_indicator(y_true)
    theta -= lr * grad

The subgradient is the difference of two pair-indicator matrices, which is why
this works with any exact decoder: our DP here, ViennaRNA soft constraints later.
``loss`` is the symmetric-difference size, so the augmented term rewards a
decoder that only wins by over-predicting.

Everything here is plain numpy and takes the decoder as a callable, so the
gradient direction can be tested without a GPU, a backbone, or ViennaRNA.
"""

from __future__ import annotations

from typing import Callable

import numpy as np

from rff.decoders.nn_dp import NEG

Decoder = Callable[[str, np.ndarray], tuple[list[tuple[int, int]], float]]


def indicator(pairs, shape: tuple[int, int]) -> np.ndarray:
    m = np.zeros(shape)
    for i, j in pairs:
        m[i, j] = 1.0
    return m


def loss_augmented(
    scores: np.ndarray, truth, mask: np.ndarray, margin: float
) -> np.ndarray:
    """Add ``margin`` to every allowed cell that is not a true pair."""
    aug = np.array(scores, dtype=float, copy=True)
    penalise = mask & ~indicator(truth, scores.shape).astype(bool)
    aug[penalise] += margin
    return aug


def margin_violation(scores: np.ndarray, predicted, truth) -> float:
    """score(predicted) - score(truth), restricted to real pair cells."""
    s = np.where(scores <= NEG, 0.0, scores)
    return float(s[indicator(predicted, scores.shape).astype(bool)].sum()
                 - s[indicator(truth, scores.shape).astype(bool)].sum())


def structured_step(
    theta: np.ndarray,
    seq: str,
    truth,
    mask: np.ndarray,
    decode: Decoder,
    margin: float = 1.0,
    lr: float = 0.1,
) -> dict:
    """One loss-augmented subgradient step; returns diagnostics, mutates ``theta``."""
    aug = loss_augmented(theta, truth, mask, margin)
    aug = np.where(mask, aug, NEG)
    pred, _ = decode(seq, aug)
    pred_set, truth_set = set(pred), set(truth)
    loss = len(pred_set ^ truth_set)
    grad = indicator(pred_set, theta.shape) - indicator(truth_set, theta.shape)
    theta -= lr * grad
    theta[...] = np.where(mask, theta, NEG)
    return {
        "predicted": sorted(pred_set),
        "loss": loss,
        "violation": margin_violation(theta, pred_set, truth_set),
        "n_pred": len(pred_set),
    }


def train_epoch(
    theta: np.ndarray,
    sequences,
    decode: Decoder,
    mask_fn: Callable[[str], np.ndarray],
    margin: float = 1.0,
    lr: float = 0.1,
) -> dict:
    """One pass over ``[(seq, truth_pairs), ...]``; returns aggregate loss."""
    total = 0
    n = 0
    for seq, truth in sequences:
        mask = mask_fn(seq)
        info = structured_step(theta, seq, truth, mask, decode, margin, lr)
        total += info["loss"]
        n += 1
    return {"mean_loss": total / n if n else float("nan"), "n": n}
