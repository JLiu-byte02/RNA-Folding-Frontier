# Amendment 8 — pair-level seed labels and primary-scoring eligibility

Status: written 2026-11-12, before any primary F1 has been computed on
`test_novel`, `id_test`, or `withheld`. It carries PI Round 4h §2–§5 and §8.
Amendment 6 and Amendment 7 remain unchanged; this file only adds.

## A8.1 — Empty truth is a label defect, not an unstructured RNA

The scorer's edge case is explicit: `prf(0, 0, 0)["f1"] == 1.0`, and the
sequence-level `pair_prf` has the same vacuous convention. That is a defensible
mathematical convention for an empty prediction against an empty label, but it
is not evidence about RNA structure and must never enter a benchmark average.

The frozen v1 split builder sets `pairs = []` when its bracket projection does
not parse (`build_s1prime` catches `ValueError`). This produced empty truth on
16,540 of 92,781 seed rows. The Round 4 withheld comparison contains 11,527 seed
rows, so masking a predictor toward no output could be rewarded by a vacuous
F1 rather than penalised. The Round 4 numbers are therefore secondary-style
diagnostics, not primary evidence.

## A8.2 — `labels_v3`: parse the alignment columns, then map to sequence

v2 projects `SS_cons` character by character and parses the projected string.
That is not equivalent to pair-level projection. If a gap removes one endpoint
of a helix, the projected string can pair the surviving endpoint with a
different closing column, fabricating a pair. A concrete example is the v2
smoke case `SS_cons="((..))"` with `aligned="A-ACG-"`: v2 returns `(0, 3)`,
while both consensus pairs have a gapped endpoint and v3 returns no pair.

`rff.data.wuss.project_pairs` now:

1. parses `SS_cons` in alignment-column coordinates, including all bracket
   types and WUSS letter pairs;
2. maps both endpoints through the ungapped positions of the row;
3. **drops** the pair if either column is a gap;
4. applies the frozen canonical AU/UA/GC/CG/GU/UG filter.

A corrupt consensus row raises `ValueError`; it is never silently turned into a
plausible label.

The artifact is a full replacement map, not a delta against a delta:

```
labels_v3_pairs = {
  "<RF accession>|<seed name>": [[i, j], ...] or null
}
```

`null` means the row failed pair-level projection or the alignment lengths did
not match. A missing key is therefore never ambiguous.

Artifact: `experiments/011_labels_v3/labels_v3_pairs.json.gz`, sha256
`b8fe71abfab30c97298bc90140bd82fc7b787f2cc439b838d86a1b0b6b895ec1`,
4,236,726 compressed bytes, 25,233,881 uncompressed bytes, gzip `mtime=0`.
Re-running `scripts/build_labels_v3.py` reproduces the hash.

### Parse and label counts

| quantity | value |
|---|---:|
| seed rows | 92,781 |
| v3 pair-level projection parses | 92,781 (100.00%) |
| v3 parse failures | 0 |
| pairs, v1 effective | 1,669,314 |
| pairs, v2 effective | 2,029,886 |
| pairs, v3 | 2,517,623 |
| pairs added v3 vs v2 effective | +487,737 |
| rows changed v3 vs v2 effective | 16,039 |
| rows with 0 v3 pairs | 4,245 |
| rows with 1–3 v3 pairs | 1,127 |
| rows with <4 v3 pairs | 5,372 (5.79%) |

The v2 “effective” count includes v1 fallback on the 9,322 rows where v2 failed
to parse. The v2 artifact itself reports 2,020,347 pairs on rows where v2
parsed; these are different denominators and must not be quoted interchangeably.

### Agreement with v2

On the 83,459 rows where both v2 and v3 parse, the pair sets differ on 6,757
rows: 27,367 pairs are v2-only and 58,637 pairs are v3-only. The PI's prior
expectation that both-parsing rows would be identical is not true on the real
data. This is a v2 projection bug, not a v3 parsing change: v2 can invent a
pair after a gap removes one endpoint, and it can also lose a genuine pair when
the projected stack is disturbed. v3 is the specification; v2 is provenance.

## A8.3 — Frozen scoring policy

For Rfam seed rows, **both primary and secondary labels are v3**. v1 and v2 are
provenance only and may not be substituted into a primary result. bpRNA labels
are unchanged by this amendment.

1. **Pair policy.** Canonical AU/UA/GC/CG/GU/UG only. Pseudoknots remain in the
   full canonical v3 set and are available to the secondary metric.
2. **Primary truth.** The maximum nested subset of the v3 canonical set, using
   the Amendment 7 polynomial DP and its deterministic lexicographic tie-break.
3. **Primary row eligibility.** A sequence with fewer than four true pairs under
   v3 is excluded from primary F1. It is not scored as empty truth, and it is
   not silently averaged as F1=1. This replaces the old family-level
   `consensus_pairs >= 4` condition with a per-sequence scoring-time rule.
   Predictions emitted on excluded rows are reported separately as a
   predicted-pair / false-positive diagnostic.
4. **Parse failures.** A row with `null` v3 labels is excluded from both
   metrics and counted in a separate parse-failure line. There are currently
   zero such rows.
5. **Secondary metric.** Secondary F1 uses the full canonical v3 set over all
   rows with v3 labels. A secondary score restricted to the primary-eligible
   rows is reported alongside it for like-for-like comparisons, never in place
   of the full-set secondary score.
6. **No outcome-dependent correction.** This label correction was made before
   any primary F1 existed. No primary score has been inspected to choose the
   parser, the eligibility threshold, or the metric.

## A8.4 — Round 4 numbers are relabelled

Every Round 4 number that used the frozen splits is labelled
**“secondary-style, v1 labels including empty-truth rows”**. These numbers may
be used to diagnose training dynamics, but:

- they are not primary F1;
- they must never appear in the same table as primary results;
- they must not be quoted as evidence for a method comparison until rerun on
  v3 with the eligibility rule above.

## A8.5 — Frozen evaluation-subset audit

`experiments/001_e1prime/test_novel_eval.jsonl` is unchanged; its locked hash
remains `9d1ab36c…` and its 2,372 records reproduce v1 exactly (0 mismatches
between the JSONL pairs and a fresh v1 reprojection).

| quantity | v1 | v3 |
|---|---:|---:|
| records | 2,372 | 2,372 |
| parse success | 2,372 | 2,372 |
| empty truth | 463 | 0 |
| fewer than 4 true pairs | 463 | 0 |
| excluded by the <4 rule | — | 0 |

All 463 empty v1 rows receive non-empty v3 labels. The `<4` rule removes **0 of
2,372 (0.00%)** evaluation records, below the 10% escalation threshold in PI
Round 4h §5, so there is no reason to delay pinning this amendment.
