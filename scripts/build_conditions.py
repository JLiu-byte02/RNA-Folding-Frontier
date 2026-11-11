#!/usr/bin/env python3
"""Generate the E1' distal-context swap conditions (b, d, e).

Condition (a) is the native sequence and needs no file: it is the cache built
from the split itself.

Region definition (documented because the PI's text leaves it implicit): the
local window is the ``--window`` nt on either side of the **sequence centre**,
``[L//2 - w, L//2 + w]``, and everything outside it is distal. The centre is
used because it is independent of the true structure; a truth-dependent
definition (distal = far from any true pair) is available with
``--window-mode paired`` as a sensitivity analysis. Only the distal region is
replaced, so every condition keeps the sequence length, the local window and
therefore the positional distribution fixed.

Conditions:

``b``  distal replaced by a dinucleotide-composition-preserving shuffle
       (Altschul-Erickson Eulerian trail over the dinucleotide multigraph), so
       nucleotide *and* dinucleotide statistics are kept but no family signal is.
``c``  distal replaced by the *aligned* distal region of a same-family homolog
       below 60% identity (``--homolog-map`` from ``build_homolog_pool.py``).
       Positions whose column aligns to a gap in the homolog keep the native
       base; the fraction of distal positions that actually came from the
       homolog is written per record as ``distal_mapped_frac``. Only records with
       a verified homolog are written, and ``csubset_native.jsonl`` holds their
       unmodified sequences so the (a) row can be scored on exactly that subset.
``d``  distal replaced by a same-length slice from a sequence of a different
       clan (neither the instance nor its family).
``e``  distal replaced by ``N`` (unknown token); this is the distortion anchor,
       not a natural-input condition.
"""

from __future__ import annotations

import argparse
import json
import random
from collections import Counter, defaultdict
from pathlib import Path

ALPHABET = "ACGU"


def positional_map(target: str, donor: str) -> list[int]:
    """``map[i]`` = donor index aligned to target index ``i``, or -1 for a gap.

    Match/mismatch scoring is used rather than a substitution matrix because
    Rfam seed rows contain IUPAC ambiguity codes that a matrix would reject.
    """
    from Bio import Align

    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -6.0
    aligner.extend_gap_score = -1.0
    aln = aligner.align(target, donor)[0]
    mapping = [-1] * len(target)
    for (ts, te), (ds, de) in zip(aln.aligned[0], aln.aligned[1]):
        if te - ts != de - ds:
            raise AssertionError("aligned blocks must have equal length")
        for off in range(te - ts):
            mapping[ts + off] = ds + off
    return mapping


def dinuc_shuffle(seq: str, rng: random.Random) -> str:
    """Altschul-Erickson: random Eulerian trail over the dinucleotide graph.

    Hierholzer's algorithm with a randomised edge order. The graph built from a
    sequence always admits a trail from its first base, so the trail exists; the
    randomised order is what makes different draws differ.
    """
    if len(seq) <= 2:
        return seq
    # Keyed by whatever symbols occur: Rfam seed alignments use IUPAC ambiguity
    # codes (S, W, R, ...), so a fixed ACGU alphabet would fail on them.
    adjacency: dict[str, list[str]] = defaultdict(list)
    for a, b in zip(seq, seq[1:]):
        adjacency[b]  # ensure the sink vertex exists
        adjacency[a].append(b)
    for targets in adjacency.values():
        rng.shuffle(targets)
    ptr = {b: 0 for b in adjacency}
    stack = [seq[0]]
    trail: list[str] = []
    while stack:
        v = stack[-1]
        if ptr[v] < len(adjacency[v]):
            w = adjacency[v][ptr[v]]
            ptr[v] += 1
            stack.append(w)
        else:
            trail.append(stack.pop())
    trail.reverse()
    assert len(trail) == len(seq), (len(trail), len(seq))
    return "".join(trail)


def distal_mask(length: int, window: int, mode: str, pairs) -> list[bool]:
    mask = [True] * length
    if mode == "centre":
        lo, hi = max(0, length // 2 - window), min(length, length // 2 + window + 1)
        for k in range(lo, hi):
            mask[k] = False
        return mask
    if mode == "paired":
        for i, j in pairs:
            for k in range(max(0, i - window), min(length, i + window + 1)):
                mask[k] = False
            for k in range(max(0, j - window), min(length, j + window + 1)):
                mask[k] = False
        return mask
    raise ValueError(mode)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--split", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--conditions", default="b,d,e")
    ap.add_argument("--window", type=int, default=32)
    ap.add_argument("--window-mode", default="centre", choices=["centre", "paired"])
    ap.add_argument("--seed", type=int, default=17)
    ap.add_argument("--donor-pool", default=None, help="jsonl of sequences for (d)")
    ap.add_argument(
        "--homolog-map",
        default=None,
        help="json from build_homolog_pool.py; required for condition c",
    )
    args = ap.parse_args()

    rng = random.Random(args.seed)
    recs = [json.loads(line) for line in open(Path(args.split)) if line.strip()]
    donors_by_clan: dict[str, list[str]] = defaultdict(list)
    pool_path = Path(args.donor_pool) if args.donor_pool else Path(args.split)
    if pool_path == Path(args.split):
        donor_recs = recs
    else:
        donor_recs = [json.loads(l) for l in open(pool_path) if l.strip()]
    for rec in donor_recs:
        donors_by_clan[rec.get("clan", "?")].append(rec["seq"])
    clan_keys = sorted(donors_by_clan)
    homolog: dict[str, dict] = (
        json.loads(Path(args.homolog_map).read_text()) if args.homolog_map else {}
    )
    seq_by_id = {r["id"]: r["seq"] for r in donor_recs}

    for cond in [c.strip() for c in args.conditions.split(",") if c.strip()]:
        rows = []
        native_rows: list[dict] = []
        for rec in recs:
            seq = rec["seq"]
            pairs = [tuple(p) for p in rec.get("pairs", [])]
            mask = distal_mask(len(seq), args.window, args.window_mode, pairs)
            out = list(seq)
            extra: dict = {}
            if cond == "b":
                idx = [k for k, d in enumerate(mask) if d]
                if idx:
                    sub = "".join(seq[k] for k in idx)
                    shuffled = dinuc_shuffle(sub, rng)
                    for pos, ch in zip(idx, shuffled):
                        out[pos] = ch
            elif cond == "c":
                if not args.homolog_map:
                    raise SystemExit("condition c needs --homolog-map")
                hit = homolog.get(rec["id"])
                if not hit:
                    continue
                donor = seq_by_id.get(hit["homolog"])
                if donor is None:
                    continue
                mapping = positional_map(seq, donor)
                moved = 0
                distal_n = 0
                for pos, is_distal in enumerate(mask):
                    if not is_distal:
                        continue
                    distal_n += 1
                    j = mapping[pos]
                    if 0 <= j < len(donor):
                        out[pos] = donor[j]
                        moved += 1
                extra = {
                    "homolog": hit["homolog"],
                    "homolog_identity": hit["identity"],
                    "distal_mapped_frac": round(moved / max(1, distal_n), 4),
                }
                native_rows.append(rec)
            elif cond == "d":
                own = rec.get("clan", "?")
                other_clans = [c for c in clan_keys if c != own]
                if other_clans:
                    donor = rng.choice(donors_by_clan[rng.choice(other_clans)])
                    n_distal = sum(mask)
                    if len(donor) > n_distal:
                        start = rng.randrange(0, len(donor) - n_distal)
                        piece = donor[start : start + n_distal]
                    else:
                        piece = (donor * (n_distal // len(donor) + 1))[:n_distal]
                    k = 0
                    for pos, is_distal in enumerate(mask):
                        if is_distal:
                            out[pos] = piece[k]
                            k += 1
            elif cond == "e":
                for pos, is_distal in enumerate(mask):
                    if is_distal:
                        out[pos] = "N"
            else:
                raise ValueError(f"unknown condition {cond!r}")
            rows.append({**rec, "seq": "".join(out), "condition": cond})
            if extra:
                rows[-1].update(extra)

        path = Path(args.out) / f"{cond}.jsonl"
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w") as fh:
            for row in rows:
                fh.write(json.dumps(row) + "\n")
        print(f"{cond}: wrote {len(rows)} -> {path}")
        if cond == "c":
            sub = Path(args.out) / "csubset_native.jsonl"
            with open(sub, "w") as fh:
                for row in native_rows:
                    fh.write(json.dumps(row) + "\n")
            print(f"c: wrote {len(native_rows)} native rows -> {sub}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
