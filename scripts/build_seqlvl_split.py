#!/usr/bin/env python3
"""Sequence-level split inside training clans: the positive control for Q2.

The clan-held-out splits (S1') ask whether the model still works when the whole
fold is new. This split asks the opposite question, and it exists to prove the
instrument can see a shortcut at all: families are shared between train and
test, only near-duplicate *sequences* are kept apart.

Construction, in the order the check has to be read:

1. Sequences are grouped by clan. Everything here stays inside clans that are
   already in the training pool, so no family is ever held out.
2. Inside a clan, sequences are clustered at ``--identity`` (default 0.8) by
   exact global alignment. Clusters are atomic, so a near-duplicate never
   straddles two splits.
3. Clusters are dealt to train/val/test so each split gets whole clusters, with
   every clan that has at least two clusters represented on both train and test.
   Clans with a single cluster can only go to train and are counted.
4. The manifest reports the largest train-test identity actually observed inside
   a clan, so "below threshold" is verified rather than assumed.

Candidate generation uses minhash buckets (k=8, 32 hashes) so a clan with
thousands of members does not become a quadratic alignment job. The bucket rule
is a filter only: the accept/reject decision always comes from the alignment.

Prediction this split is built to test (PI Round 4, Q2): condition (c), a
low-identity same-family donor, should help more than condition (d), a
different-clan donor — the reverse of S1'. If it does not, the instrument
cannot detect family shortcuts and the S1' null is uninformative.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402


def read_jsonl(path: Path) -> list[dict]:
    with open(path) as fh:
        return [json.loads(line) for line in fh if line.strip()]


def kmer_hashes(seq: str, k: int = 8) -> set[int]:
    if len(seq) < k:
        return {hash(seq) & 0xFFFFFFFF}
    return {
        int.from_bytes(
            hashlib.blake2b(seq[i : i + k].encode(), digest_size=4).digest(), "big"
        )
        for i in range(len(seq) - k + 1)
    }


def minhash(hashes: set[int], n_hashes: int = 32) -> list[int]:
    ordered = sorted(hashes)
    if len(ordered) <= n_hashes:
        return ordered
    idx = np.linspace(0, len(ordered) - 1, n_hashes).astype(int)
    return [ordered[i] for i in idx]


class UnionFind:
    def __init__(self, n: int) -> None:
        self.parent = list(range(n))

    def find(self, a: int) -> int:
        while self.parent[a] != a:
            self.parent[a] = self.parent[self.parent[a]]
            a = self.parent[a]
        return a

    def union(self, a: int, b: int) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[rb] = ra


def identity(aligner, a: str, b: str, band: float) -> float:
    """Alignment identity over the aligned columns, or 0 if clearly below band.

    The k-mer Jaccard is used as a fast rejection test first: for sequences that
    share almost no 8-mers, an alignment cannot reach the threshold.
    """
    sa, sb = kmer_hashes(a), kmer_hashes(b)
    union = len(sa | sb)
    if union and len(sa & sb) / union < max(0.05, band - 0.5):
        return 0.0
    aln = aligner.align(a, b)[0]
    same = 0
    cols = 0
    for (ts, te), (ds, de) in zip(aln.aligned[0], aln.aligned[1]):
        width = te - ts
        cols += width
        same += sum(
            1 for off in range(width) if a[ts + off] == b[ds + off]
        )
    return same / cols if cols else 0.0


def candidate_pairs(sigs: list[list[int]], full_below: int = 60, min_shared: int = 3):
    n = len(sigs)
    if n <= full_below:
        for i in range(n):
            for j in range(i + 1, n):
                yield i, j
        return
    buckets: dict[int, list[int]] = defaultdict(list)
    for idx, sig in enumerate(sigs):
        for h in sig:
            buckets[h].append(idx)
    seen: set[tuple[int, int]] = set()
    for members in buckets.values():
        if len(members) > 200:
            # A very common hash is uninformative; keep it but cap the blow-up.
            members = members[:200]
        for a_i in range(len(members)):
            for b_i in range(a_i + 1, len(members)):
                i, j = members[a_i], members[b_i]
                key = (i, j) if i < j else (j, i)
                if key in seen:
                    continue
                seen.add(key)
                shared = sum(
                    1 for x, y in zip(sigs[key[0]], sigs[key[1]]) if x == y
                )
                if shared >= min_shared:
                    yield key


def cluster_clan(recs: list[dict], aligner, threshold: float, cache: dict) -> list[list[int]]:
    sigs = []
    for rec in recs:
        sig = cache.get(rec["id"])
        if sig is None:
            sig = minhash(kmer_hashes(rec["seq"]))
            cache[rec["id"]] = sig
        sigs.append(sig)
    uf = UnionFind(len(recs))
    for i, j in candidate_pairs(sigs):
        if uf.find(i) == uf.find(j):
            continue
        if identity(aligner, recs[i]["seq"], recs[j]["seq"], threshold) >= threshold:
            uf.union(i, j)
    groups: dict[int, list[int]] = defaultdict(list)
    for idx in range(len(recs)):
        groups[uf.find(idx)].append(idx)
    return sorted(groups.values(), key=lambda g: (-len(g), g[0]))


def deal_clusters(
    clusters: list[list[int]],
    counts: dict[str, int],
    targets: dict[str, float],
    rng: np.random.Generator,
    require: tuple[str, ...] = ("train", "test"),
) -> dict[int, str]:
    """Assign whole clusters to splits, balancing sizes against the targets."""
    order = list(range(len(clusters)))
    rng.shuffle(order)
    order.sort(key=lambda c: -len(clusters[c]))
    assigned: dict[int, str] = {}
    for pos, c in enumerate(order):
        remaining = len(order) - pos
        missing = [s for s in require if counts[s] == 0]
        if remaining == len(missing) and missing:
            name = missing[0]
        else:
            name = max(
                targets,
                key=lambda s: targets[s] * (sum(counts.values()) + len(clusters[c]))
                - counts[s],
            )
        assigned[c] = name
        counts[name] += len(clusters[c])
    return assigned


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--train", default="splits/s1prime_filtered/train.jsonl")
    ap.add_argument("--out", default="splits/seqlvl_positive_control")
    ap.add_argument("--identity", type=float, default=0.80)
    ap.add_argument("--val-frac", type=float, default=0.10)
    ap.add_argument("--test-frac", type=float, default=0.15)
    ap.add_argument("--seed", type=int, default=17)
    ap.add_argument("--max-clans", type=int, default=None)
    args = ap.parse_args()

    from Bio import Align

    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -6.0
    aligner.extend_gap_score = -1.0

    recs = read_jsonl(ROOT / args.train)
    by_clan: dict[str, list[dict]] = defaultdict(list)
    for rec in recs:
        by_clan[rec["clan"]].append(rec)
    clans = sorted(by_clan)
    if args.max_clans:
        clans = clans[: args.max_clans]

    targets = {
        "train": 1.0 - args.val_frac - args.test_frac,
        "val": args.val_frac,
        "test": args.test_frac,
    }
    rng = np.random.default_rng(args.seed)
    cache: dict[str, list[int]] = {}
    total = {"train": 0, "val": 0, "test": 0}
    split_rows: dict[str, list[dict]] = {"train": [], "val": [], "test": []}
    info: dict[str, dict] = {}
    single_cluster_clans = 0
    max_cross = 0.0

    for clan in clans:
        members = by_clan[clan]
        clusters = cluster_clan(members, aligner, args.identity, cache)
        info[clan] = {
            "n_sequences": len(members),
            "n_clusters": len(clusters),
            "cluster_sizes": [len(c) for c in clusters],
        }
        if len(clusters) == 1:
            single_cluster_clans += 1
            for idx in clusters[0]:
                split_rows["train"].append(members[idx])
            total["train"] += len(members)
            info[clan]["splits"] = {"train": len(members)}
            continue

        counts = {"train": 0, "val": 0, "test": 0}
        assigned = deal_clusters(clusters, counts, targets, rng)
        per_split = defaultdict(int)
        for c_idx, split in assigned.items():
            for idx in clusters[c_idx]:
                split_rows[split].append(members[idx])
                per_split[split] += 1
        for split, n in per_split.items():
            total[split] += n
        info[clan]["splits"] = dict(per_split)

        # Verify the threshold holds across train and test. The same candidate
        # rule as clustering is used, so this verifies the split given the
        # filter; the accept/reject decision is still the aligner's.
        side = {}
        for c_idx, s in assigned.items():
            for idx in clusters[c_idx]:
                side[idx] = s
        sigs = [cache[m["id"]] for m in members]
        for i, j in candidate_pairs(sigs):
            if side[i] == side[j]:
                continue
            v = identity(aligner, members[i]["seq"], members[j]["seq"], 0.0)
            if v > max_cross:
                max_cross = v
                info[clan]["max_cross_identity"] = {
                    "value": round(float(v), 4),
                    "a": members[i]["id"],
                    "b": members[j]["id"],
                }
            if v >= args.identity:
                raise AssertionError(
                    f"cluster split leaked a >= {args.identity} pair in clan {clan}: "
                    f"{members[i]['id']} vs {members[j]['id']} at {v:.3f}"
                )

    out_dir = ROOT / args.out
    out_dir.mkdir(parents=True, exist_ok=True)
    for split, rows in split_rows.items():
        if not rows:
            continue
        with open(out_dir / f"{split}.jsonl", "w") as fh:
            for row in rows:
                fh.write(json.dumps(row) + "\n")
        print(f"{split}: {len(rows)} sequences, {len({r['family'] for r in rows})} families")
    manifest = {
        "source": args.train,
        "identity_threshold": args.identity,
        "seed": args.seed,
        "n_clans": len(clans),
        "n_single_cluster_clans": single_cluster_clans,
        "counts": total,
        "families": {s: len({r["family"] for r in rows}) for s, rows in split_rows.items()},
        "max_train_test_identity": round(float(max_cross), 4),
        "note": (
            "Positive control: families are shared across splits by construction; "
            "only sequences at or above the identity threshold are kept apart."
        ),
        "per_clan": info,
    }
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({k: manifest[k] for k in ("counts", "families", "max_train_test_identity")}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
