"""Nearest-neighbour dynamic programming decoder ("decoder A").

Purpose
-------
P2's decoder reads pairs off the probability map under a greedy
one-partner/non-crossing rule, which is not a structure model: a pair is kept on
its own evidence with no notion of what the pair closes. Decoder A replaces it
with exact DP over the standard nearest-neighbour decomposition (hairpin,
interior/bulge, multiloop, exterior), so that a pair is only accepted when the
substructure it closes is consistent — the H2.4 comparison in Amendment 2.

What this module is *not*
-------------------------
It is not yet a ViennaRNA-parameter-faithful energy model. The loop terms are
explicit, swappable arrays, and the default is the zero-penalty / learned-score
configuration used when decoding a probability map. Turner-2004 parity (special
hairpins, terminal AU, Ninio/logistic multiloop, dangle models) is a separate,
later task with its own cross-check against ``RNAfold``; claiming physics
before that check would be exactly the kind of unverified claim this project
forbids.

Score convention
----------------
All quantities are **scores to maximise**. ``pair_score[i, j]`` is the score of
the pair itself (for a probability map: the log-odds of that cell), loop terms
are non-positive penalties in the same units, and ``-inf`` marks a disallowed
configuration. With ``-pair_score`` and Turner tables this is the usual
MFE-as-minimisation formulation.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

NEG = -1e9
VALID = NEG / 2.0


def pairs_to_dotbracket(length: int, pairs) -> str:
    """Convert a pair list to dot-bracket notation."""
    chars = ["."] * length
    for i, j in pairs:
        if i >= j:
            raise ValueError(f"pair is not ordered: {(i, j)}")
        chars[i] = "("
        chars[j] = ")"
    return "".join(chars)


@dataclass
class PairScoreModel:
    """Loop penalties plus an optional per-pair score matrix.

    ``hairpin[k]`` is the penalty for a hairpin loop of ``k`` unpaired
    nucleotides; ``interior[k, m]`` is the penalty for a loop of ``k`` unpaired
    nucleotides split ``m`` / ``k - m`` between the two sides. Multiloop terms
    are added once per closing pair, once per branch, and once per unpaired
    position.
    """

    pair_score: np.ndarray | None = None
    hairpin: np.ndarray | None = None
    interior: np.ndarray | None = None
    ml_closing: float = 0.0
    ml_branch: float = 0.0
    ml_unpaired: float = 0.0
    min_loop: int = 3
    max_loop: int = 30

    def hairpin_penalty(self, loop_len: int) -> float:
        if self.hairpin is None:
            return 0.0
        if loop_len >= len(self.hairpin):
            return float(self.hairpin[-1])
        return float(self.hairpin[loop_len])

    def interior_penalty(self, loop_len: int, asymmetry: int) -> float:
        if self.interior is None:
            return 0.0
        loop_len = min(loop_len, self.interior.shape[0] - 1)
        asymmetry = min(asymmetry, self.interior.shape[1] - 1)
        return float(self.interior[loop_len, asymmetry])


class NNDecoder:
    """O(n^3) nearest-neighbour DP with traceback.

    Disallowed pairs must be ``-inf`` in ``pair_score`` (or the model decides);
    positions whose pair score is ``-inf`` are never paired.
    """

    def __init__(self, model: PairScoreModel):
        self.model = model

    def _offsets(self):
        m = self.model
        dks, dls = [], []
        # dk, dl >= 1: the inner pair must be strictly inside (i, j).
        for dk in range(1, m.max_loop + 2):
            for dl in range(1, m.max_loop + 2):
                if dk + dl - 2 > m.max_loop:
                    continue
                dks.append(dk)
                dls.append(dl)
        dks = np.array(dks, dtype=np.int64)
        dls = np.array(dls, dtype=np.int64)
        pens = np.array(
            [
                m.interior_penalty(dk + dl - 2, abs(dk - dl))
                for dk, dl in zip(dks, dls)
            ],
            dtype=np.float64,
        )
        return dks, dls, pens

    def decode(self, seq: str, pair_score: np.ndarray | None = None) -> tuple[list[tuple[int, int]], float]:
        n = len(seq)
        model = self.model
        S = model.pair_score if pair_score is None else pair_score
        if S is None:
            raise ValueError("no pair score matrix: pass pair_score= or set model.pair_score")
        if S.shape != (n, n):
            raise ValueError(f"pair_score must be {n}x{n}, got {S.shape}")
        min_loop = model.min_loop
        if n <= min_loop + 1:
            return [], 0.0

        dks, dls, pens = self._offsets()
        V = np.full((n, n), NEG)
        A0 = np.full((n, n), NEG)  # multiloop content, >= 0 branches
        A1 = np.full((n, n), NEG)  # multiloop content, >= 1 branch
        W = np.zeros((n, n))
        Wtb = np.ones((n, n), dtype=np.int64)  # default: i unpaired
        Vtb = np.zeros((n, n), dtype=np.int64)
        A0tb = np.zeros((n, n), dtype=np.int64)
        A1tb = np.zeros((n, n), dtype=np.int64)
        Vk = np.zeros((n, n), dtype=np.int64)
        Vl = np.zeros((n, n), dtype=np.int64)
        A0k = np.zeros((n, n), dtype=np.int64)
        A1k = np.zeros((n, n), dtype=np.int64)
        Wk = np.zeros((n, n), dtype=np.int64)

        # Base case for the multiloop content: every position unpaired.
        for i in range(n):
            for j in range(i, n):
                A0[i, j] = (j - i + 1) * model.ml_unpaired
                A0tb[i, j] = 4

        for d in range(min_loop + 1, n):
            for i in range(n - d):
                j = i + d
                s = float(S[i, j])
                if s > NEG:
                    # --- V: structure closed by (i, j) ---------------------
                    best = model.hairpin_penalty(d - 1)
                    code, tk, tl = 0, 0, 0

                    mask = (dks + dls) < d
                    if mask.any():
                        kk = i + dks[mask]
                        ll = j - dls[mask]
                        vals = np.where(
                            V[kk, ll] > VALID, V[kk, ll] + pens[mask], NEG
                        )
                        pos = int(np.argmax(vals))
                        if vals[pos] > best:
                            best = float(vals[pos])
                            idx = int(np.flatnonzero(mask)[pos])
                            code, tk, tl = 1, int(dks[idx]), int(dls[idx])

                    if j - 1 > i + 1 and A1[i + 1, j - 1] > VALID:
                        cand = model.ml_closing + A1[i + 1, j - 1]
                        if cand > best:
                            best = cand
                            code, tk, tl = 2, 0, 0

                    V[i, j] = s + best
                    Vtb[i, j] = code
                    Vk[i, j] = tk
                    Vl[i, j] = tl

                # --- A0/A1: multiloop content on [i, j] -------------------
                # branch term: (i, k) is the first branch, then any content
                ks = np.arange(i + min_loop + 1, j + 1)
                branch_best, branch_k = NEG, 0
                if ks.size:
                    rest = np.where(
                        ks < j, A0[np.minimum(ks + 1, j), j], 0.0
                    )
                    vals = np.where(
                        V[i, ks] > VALID,
                        V[i, ks] + model.ml_branch + rest,
                        NEG,
                    )
                    pos = int(np.argmax(vals))
                    branch_best, branch_k = float(vals[pos]), int(ks[pos])

                # i unpaired, then content on [i+1, j]
                cand0 = A0[i + 1, j] + model.ml_unpaired if i + 1 <= j else NEG
                if branch_best >= cand0:
                    A0[i, j], A0tb[i, j], A0k[i, j] = branch_best, 2, branch_k
                else:
                    A0[i, j], A0tb[i, j], A0k[i, j] = cand0, 1, 0

                cand1 = (
                    A1[i + 1, j] + model.ml_unpaired
                    if i + 1 <= j and A1[i + 1, j] > VALID
                    else NEG
                )
                if branch_best >= cand1:
                    A1[i, j], A1tb[i, j], A1k[i, j] = branch_best, 2, branch_k
                else:
                    A1[i, j], A1tb[i, j], A1k[i, j] = cand1, 1, 0

                # --- W: exterior on [i, j] --------------------------------
                cand = W[i + 1, j] if j > i else 0.0
                wcode, wk = 1, 0
                ks = np.arange(i + min_loop + 1, j + 1)
                if ks.size:
                    rest = np.where(ks < j, W[np.minimum(ks + 1, j), j], 0.0)
                    vals = np.where(V[i, ks] > VALID, V[i, ks] + rest, NEG)
                    pos = int(np.argmax(vals))
                    if vals[pos] > cand:
                        cand = float(vals[pos])
                        wcode, wk = 2, int(ks[pos])
                W[i, j] = cand
                Wtb[i, j] = wcode
                Wk[i, j] = wk

        self.tables = {"V": V, "A0": A0, "A1": A1, "W": W}

        pairs: list[tuple[int, int]] = []
        stack = [(0, n - 1, "W")]
        steps = 0
        while stack:
            steps += 1
            if steps > 20 * n:
                raise RuntimeError(
                    f"traceback did not terminate after {steps} steps; "
                    f"next state {stack[-1]}"
                )
            i, j, kind = stack.pop()
            if i > j:
                continue
            if kind == "W":
                if Wtb[i, j] == 1:
                    stack.append((i + 1, j, "W"))
                else:
                    k = int(Wk[i, j])
                    stack.append((i, k, "V"))
                    stack.append((k + 1, j, "W"))
                continue
            if kind in ("A0", "A1"):
                code = (A0tb if kind == "A0" else A1tb)[i, j]
                kk = int((A0k if kind == "A0" else A1k)[i, j])
                if code == 4:
                    pass
                elif code == 1:
                    stack.append((i + 1, j, kind))
                else:
                    stack.append((i, kk, "V"))
                    stack.append((kk + 1, j, "A0"))
                continue
            pairs.append((i, j))
            code = Vtb[i, j]
            if code == 1:
                stack.append((i + int(Vk[i, j]), j - int(Vl[i, j]), "V"))
            elif code == 2:
                stack.append((i + 1, j - 1, "A1"))

        pairs.sort()
        return pairs, float(W[0, n - 1])


def logit_score(prob: np.ndarray, eps: float = 1e-6) -> np.ndarray:
    """Log-odds transform of a probability map, symmetric."""
    p = np.clip(0.5 * (prob + prob.T), eps, 1 - eps)
    return np.log(p / (1 - p))


def mask_allowed(seq: str, min_loop: int = 3) -> np.ndarray:
    """Boolean matrix of canonical or GU pairs with ``j - i > min_loop``."""
    n = len(seq)
    mask = np.zeros((n, n), dtype=bool)
    allowed = {"AU", "UA", "GC", "CG", "GU", "UG"}
    for i in range(n):
        for j in range(i + min_loop + 1, n):
            mask[i, j] = seq[i] + seq[j] in allowed
    return mask
