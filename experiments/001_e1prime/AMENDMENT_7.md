# Amendment 7 — lock basis, nested-subset solver, labels_v2, and the (s, b) fitting protocol

Status: written 2026-11-12, before any `(s, b)` fit and before any scoring of
`id_test` or `withheld`. It carries PI Round 4g rulings §1–§3. Amendment 6 is
unchanged; this file only adds.

## A7.1 — A locked hash is a hash of git blob bytes

Amendment 6 §A6.5 pinned sha256 values, but did not say over which bytes. On this
machine `core.autocrlf=true`, so a working-tree file can be CRLF while its git
blob is LF; a lock taken over the wrong one would pass here and fail on a Linux
checkout of the same commit.

Ruling: every lock is the sha256 of `git cat-file blob HEAD:<path>`. The guard
(`src/rff/analysis/embargo.py`) additionally refuses when the working tree no
longer equals that blob, because in that case the file being read is not the file
that was frozen — either it was edited after the fact or a line-ending rewrite
changed it.

Verified 2026-11-12 with `scripts/verify_lock_basis.py`
(`experiments/006_lock_basis/lock_basis_check.json`). For all five locked files,
the working-tree hash, the HEAD blob hash and every historical blob hash agree,
and no file contains CRLF:

| file | blob sha256 (first 8) | first quoted in |
|---|---|---|
| `AMENDMENT_2.md` | `4a9d2f48` | `a8e4d48242d6` |
| `AMENDMENT_3.md` | `7e035123` | `207a8eb560cd` |
| `AMENDMENT_4.md` | `08d6771a` | `3a16258146fc`, restored in `d6a2abb` |
| `AMENDMENT_5.md` | `cba17608` | `0117538d80a7` |
| `test_novel_eval.jsonl` | `9d1ab36c` | `b40d22e551e3` |

The commit in which each hash was first quoted is recorded in the lockfile, so
the claim is checkable rather than asserted. The same check must be re-run from a
fresh Linux clone when the cluster returns; that run is part of A6.1.

## A7.2 — Primary-F1 reduction is a polynomial DP over a conflict graph

The previous `max_nested_subset` enumerated every subset of a crossing component
of 16 pairs or fewer and, above that, fell back to greedy by shortest span. Both
halves are replaced.

**Conflict, not just crossing.** Two pairs are incompatible if they interleave
*or* share an endpoint. Only the first is a pseudoknot; the second is a base with
two partners, which no secondary structure permits. `crosses` alone let the
exhaustive reference keep `(17, 27)` and `(18, 27)` together, and it was the DP's
disagreement with that reference which exposed the hole.

**Solver.** Pairs that conflict with nothing are kept. The rest form a conflict
graph whose components are independent, and each component is solved exactly by
Nussinov DP restricted to the given arcs:

```
f[i][j] = max( f[i][j-1],
               max_{(k,j) in A, k >= i} f[i][k-1] + 1 + f[k+1][j-1] )
```

with `O(n^2 + n|A|)` cells. Tie-break unchanged: among maximum-size solutions,
the lexicographically smallest sorted pair list, obtained by walking arcs in
ascending order and keeping an arc when the enclosing interval's optimum is
preserved.

**Cross-check on the real labels** (`experiments/007_nested_dp/nested_subset_check.json`):

| source | sequences | conflicting sequences | components | compared with exhaustive | mismatches | components > 16 |
|---|---:|---:|---:|---:|---:|---:|
| bpRNA-1m mirror | 66,715 | 15,794 | 36,423 | 33,714 | **0** | 2,709 |
| Rfam 15.0 seeds (cap 20/family) | 30,370 | 5 | 5 | 5 | **0** | 0 |

The largest component has **326** arcs, so "pseudoknot components are tiny" was
wrong as a general assumption and the old path was being exercised more often
than believed. On the 2,709 components above the old limit, the greedy fallback
was suboptimal on **87**, losing **692** pairs in total (12 vs 10 on one
component of 18 arcs); it never produced an invalid matching on this data.

## A7.3 — `labels_v2`: corrected seed labels as a separate artifact

Membership stays frozen. `scripts/build_labels_v2.py` reads the local
`Rfam.seed.gz` and writes a delta keyed `<RF accession>|<seed name>`:

```
pairs_v2 = v1_pairs if id not in delta else delta[id]
```

| quantity | value |
|---|---:|
| seed rows | 92,781 |
| rows whose v2 pair set differs from v1 | 20,931 |
| pairs v1 | 1,659,775 |
| pairs v2 | 2,020,347 |
| pairs added | 360,572 |
| rows whose v1 projection does not parse (empty truth in the splits) | 16,540 |
| of those, parsed by v2 | 7,670 |

Artifact: `experiments/008_labels_v2/labels_v2_delta.json.gz`, sha256
`27bbaa32d11837c77b95207c8eac00dbc81890d55eb72e4f1f509f0626deb0af`, gzipped with
`mtime=0` so the bytes are reproducible (re-running gives the same hash).

Usage rule: primary F1 keeps using the frozen v1 labels; secondary F1 on seeds
may use v2; every number must state which version it used. Per-split added-pair
counts need the membership lists and are produced by
`scripts/build_labels_v2.py --split-dir splits/s1prime` under A6.1.

## A7.4 — `(s, b)` fitting protocol, pre-registered before any fit

The A-energy-fixed arm turns head scores into ViennaRNA pseudo-energies as
`e = -(s * score + b)` kcal/mol. Fitting `(s, b)` means running many
soft-constrained folds, so the search is pre-registered here, before the first
run, rather than chosen after seeing the cost or the curve.

1. **Objective.** Family-macro primary F1 on the training split, 5-fold
   cross-validation grouped by clan so no clan is split across folds.
2. **Coarse grid.** `s ∈ {0.25, 0.5, 1, 2, 4}` (kcal/mol per log-odds unit) ×
   `b ∈ {-6, -3, 0, 3, 6}` (kcal/mol). The offset range covers ±2·ln 20 ≈ ±6,
   because the head is trained with `pos_weight = 20`.
3. **Refinement.** A second 5×5 grid spanning ± one coarse step around the coarse
   optimum on each axis.
4. **Cost control.** First measure seconds per fold by length bin. If the
   projected total exceeds 24 CPU-hours, fit on a family-stratified subsample of
   3,000 training sequences drawn with seed 17, and record the subsample hash and
   the decision. The same subsample is then used for both folds and refinement.
5. **Freeze.** The fitted `(s, b)` is written to JSON, hashed, and added to
   `state/amendments.lock.json` **before** any scoring of `id_test` or
   `withheld`. No refit after seeing those numbers.

## A7.5 — Defects recorded while writing this amendment

- The split manifest's "pseudoknots dropped" sentence is false; corrected for
  future rebuilds in `build_s1prime.py` and recorded in `docs/corrections.md`.
- 16,540 seed rows carry an empty truth in the frozen splits because their v1
  projection does not parse; `labels_v2` repairs 7,670 of them. Not fixed in
  place, because that would change frozen split hashes.
- WUSS letter pairs were dropped by the v1 projection: 59,632 pairs on 13,280
  rows that otherwise parse.
