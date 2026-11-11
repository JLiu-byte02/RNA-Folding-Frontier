#!/usr/bin/env python3
"""Verify the three assumptions in `rff.decoders.vienna_soft` against real ViennaRNA.

Amendment 3 §A3.1 makes A-energy (Turner + `sc_add_bp`) the confirmatory decoder,
so its wiring has to be checked before any number depends on it. Three things
would silently corrupt results, and each has a direct test:

1. **Zero corrections must reproduce RNAfold.** A logic bug in the wrapper, a
   wrong index base or an accidental energy offset shows up immediately.
2. **A bonus must attract the pair it is given.** A sign error would push the
   decoder away from exactly the pairs the head likes — the failure mode that
   would look like "the method does not work".
3. **Indices must be 0-based.** If the bindings were 1-based, a bonus on (i, j)
   would move a neighbouring pair instead; the test asks for the exact pair.

Writes a JSON record; exits non-zero if any check fails.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import RNA  # noqa: E402

from rff.decoders.nn_dp import NEG, mask_allowed  # noqa: E402
from rff.decoders.vienna_soft import ViennaSoftDecoder  # noqa: E402


def random_sequences(rng: np.random.Generator, n: int, lengths) -> list[str]:
    return ["".join(rng.choice(list("ACGU"), int(L))) for L in rng.choice(lengths, n)]


def gc_rich(rng: np.random.Generator, n: int, length: int) -> list[str]:
    out = []
    for _ in range(n):
        out.append("".join(rng.choice(list("GGCCAU"), length)))
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="experiments/003_decoder_a/parity_vienna.json")
    ap.add_argument("--n-random", type=int, default=40)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    seqs = random_sequences(rng, args.n_random, [20, 30, 40, 60, 80, 120])
    seqs += gc_rich(rng, 20, 60)

    dec = ViennaSoftDecoder(scale=1.0)
    parity = {"n": len(seqs), "mismatched": 0, "examples": []}
    attraction = {"n": 0, "failures": 0, "examples": [], "moved_structure": 0}
    index_base = {"n": 0, "wrong_pair": 0, "examples": []}

    for seq in seqs:
        n = len(seq)
        ref_struct, ref_energy = RNA.fold(seq)
        ref_pairs = set(dec.decode(seq, np.zeros((n, n)))[0])
        ref_pairs_ref = RNA.ptable(ref_struct)
        ref_set = {(i - 1, ref_pairs_ref[i] - 1) for i in range(1, n + 1) if ref_pairs_ref[i] > i}
        if ref_pairs != ref_set or abs(dec.last["energy"] - ref_energy) > 1e-6:
            parity["mismatched"] += 1
            if len(parity["examples"]) < 5:
                parity["examples"].append(
                    {"seq": seq, "rnafold": ref_struct, "ours": sorted(ref_pairs), "energy": dec.last["energy"]}
                )

        mask = mask_allowed(seq)
        if not mask.any():
            continue
        ii, jj = np.where(mask)
        pick = int(rng.integers(0, len(ii)))
        i, j = int(ii[pick]), int(jj[pick])
        scores = np.full((n, n), NEG)
        scores[i, j] = 20.0  # 20 kcal/mol of attraction: must dominate Turner terms
        pairs, energy = dec.decode(seq, scores)
        attraction["n"] += 1
        if (i, j) in set(pairs):
            if set(pairs) != ref_pairs:
                attraction["moved_structure"] += 1
        else:
            attraction["failures"] += 1
            if len(attraction["examples"]) < 5:
                attraction["examples"].append({"seq": seq, "asked": [i, j], "got": sorted(pairs), "energy": energy})

        # a bonus one position off must not attract the asked pair as strongly:
        # if the bindings were 1-based, the previous test would have moved (i-1,j-1)
        index_base["n"] += 1
        wrong = np.full((n, n), NEG)
        if i + 1 < j and mask[i + 1, j - 1]:
            wrong[i + 1, j - 1] = 20.0
            wrong_pairs, _ = dec.decode(seq, wrong)
            if (i, j) in set(wrong_pairs) and (i + 1, j - 1) not in set(wrong_pairs):
                index_base["wrong_pair"] += 1

    ok = (
        parity["mismatched"] == 0
        and attraction["failures"] == 0
        and index_base["wrong_pair"] == 0
    )
    report = {
        "vienna_version": getattr(RNA, "__version__", "unknown"),
        "scale": dec.scale,
        "zero_score_parity": parity,
        "bonus_attracts_asked_pair": attraction,
        "index_base_probe": index_base,
        "passes": ok,
    }
    out = ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: report[k] for k in ("vienna_version", "passes")}))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
