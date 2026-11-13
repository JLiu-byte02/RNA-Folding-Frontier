# External baselines on the v3 protocol and the common-cohort policy - 2026-11-13

Status: measurement only. RNAfold, EternaFold and RNAstructure are complete on a
single frozen common cohort. MXfold2 (job 15631, retrained on our train split)
is still running, so the final five-method cohort is not yet frozen. Nothing in
this file is a new method, and no row here is evidence that our head beats a
thermodynamic baseline.

## 1. The defect this policy fixes

`scripts/score_external_baseline.py` skips a sequence when a prediction is
missing. Each tool has its own practical limit, so each baseline was scored on a
different population:

- EternaFold: 1,200 nt cap (CONTRAfold is O(L^3)); 578 withheld rows dropped;
- RNAstructure: 111 rows with IUPAC ambiguity codes are rejected outright by
  `Fold` ("Invalid nucleobase"), of which 104 are inside our label-eligible dump;
- MXfold2: trained and folded with a 600 nt cap;
- our own heads: a 1,024 nt context limit.

A family-macro column computed over different populations is not comparable.
The fix is one explicit id list that every method - our heads included - is
rescored against with `--restrict-ids`.

## 2. Cohort definition

A sequence enters the external-baseline cohort iff all of the following hold:

1. it appears in our model dump `dump_p2fix_seed17/a_<split>.jsonl`, so our
   heads can be scored on it;
2. its `truth_v3` is not null (Amendment 8 parse-failure rule) and has at least
   `--min-true-pairs 4` pairs **before** the nested projection, exactly as
   `rff.analysis.block_eval.score_dump` applies it;
3. a prediction exists in **every** prediction file handed to
   `scripts/build_common_cohort.py`.

Conditions 1 and 2 are label/coverage properties fixed before any score was
computed. Condition 3 is a mechanical intersection. No condition uses an
observed score, and no per-method subset is ever reported as a cohort.

## 3. Coverage accounting

The binding limit is **our 1,024 nt context**, not any tool's cap. Sequences are
dropped for length only when they exceed it, and the length distributions of the
kept and dropped sets do not overlap.

| split | split rows | dump rows | missing | missing % | missing len min / median / max | longest present |
|---|---:|---:|---:|---:|---|---:|
| test_withheld | 17,371 | 16,793 | 578 | 3.327% | 1,275 / 2,990 / 4,381 | 1,014 |
| id_test | 1,402 | 1,392 | 10 | 0.713% | 1,037 / 2,531.5 / 2,729 | 990 |

Families:

| split | split families | lost to coverage | lost to v3 `<4` pairs | surviving |
|---|---:|---:|---:|---:|
| test_withheld | 444 | 2 | 33 | 409 |
| id_test | 168 | 0 | 13 | 155 |

- Fully lost to coverage (withheld): **RF02540** (169/169) and **RF02543**
  (174/174). Partially lost: RF01999 (1/234), RF02541 (234/358).
- id_test loses **no family** to coverage; RF00029 loses 10/65 rows. The 13
  families that disappear from the id_test table disappear because every member
  has fewer than four v3 pairs, not because of length.
- After the RNAstructure IUPAC rule, the withheld cohort drops a further 104
  sequences (16,163 -> 16,059). id_test loses none.

v3 label eligibility inside the dump:

| split | empty truth | 1-3 pairs | `<4` total | `>=4` pairs (label-eligible) | pair p25 / median / p75 |
|---|---:|---:|---:|---:|---|
| test_withheld | 118 | 512 | 630 | 16,163 | 14 / 24 / 34 |
| id_test | 5 | 100 | 105 | 1,287 | 16 / 24 / 32 |

The `<4` class is **not** "rows with 1-3 pairs only": it is the union of
`empty_truth` and `lt_min_true_pairs`, and the two are reported separately in
`cohort/<split>.report.json` (`excluded_by_reason`). An empty truth row is a
label defect (no parseable `SS_cons` projection), not unstructured RNA, and it is
excluded rather than scored as all-negative.

Tool-cap overlap with the 1,024 nt gap:

- withheld: 578 rows are longer than 1,024 nt and 578 are longer than 1,200 nt.
  The sets coincide exactly, so EternaFold's cap removes nothing our context
  did not already remove.
- id_test: 10 rows are longer than 1,024 nt but only 7 are longer than 1,200 nt.
  EternaFold folds the three rows in (1,024, 1,200] and they are then dropped by
  the restrict set. Our context is still the binding limit.

## 4. Per-tool exclusions

| tool | rule | withheld rows | id_test rows | where recorded |
|---|---|---:|---:|---|
| RNAfold 2.8.0 | none | 0 | 0 | `rnafold_*_score.json` |
| EternaFold (CONTRAfold-SE) | `--max-length 1200` | 578 | 7 (+3 dropped by restrict set) | job 15653 log, `n_outside_restrict_set` |
| RNAstructure 6.6 | `--skip-ambiguous` (IUPAC) + `--allow-partial` timeouts | 111 skipped, 104 in cohort | 0 | `<out>.skipped.jsonl`, `<out>.errors.jsonl` |
| MXfold2 (MixC, retrained) | 600 nt train cap; test cap set at export | pending | pending | `experiments/013_mxfold2/run/export_manifest.json` |
| P2-fix heads (ours) | 1,024 nt context | 578 | 10 | `coverage/{withheld,id_test}.json` |

RNAstructure's 111 ambiguous rows are a *tool* limitation, not a length limit
and not a timeout. The 41 timeouts are all longer than 3,337 nt and were never in
the dump, so they are outside this cohort by construction. Do not describe the
104-row gap as slowness.

## 5. Common-cohort results

One cohort per split, identical id lists for every method: withheld 16,059
sequences / 409 families, id_test 1,287 / 155. v3 labels, `--min-true-pairs 4`,
nested truth, 10,000 family-block bootstrap, seed 20261113, 95% intervals.

### Withheld

| method | n | families | family-macro F1 | 95% CI | micro F1 |
|---|---:|---:|---:|---|---:|
| RNAfold | 16,059 | 409 | 0.6004 | [0.5744, 0.6263] | 0.5560 |
| EternaFold | 16,059 | 409 | 0.6050 | [0.5792, 0.6301] | 0.5908 |
| RNAstructure | 16,059 | 409 | 0.5901 | [0.5644, 0.6159] | 0.5482 |
| head seed17 | 16,059 | 409 | 0.6054 | [0.5808, 0.6290] | 0.6040 |
| head seed23 | 16,059 | 409 | 0.5785 | [0.5539, 0.6020] | 0.5897 |
| head seed41 | 16,059 | 409 | 0.5854 | [0.5616, 0.6088] | 0.5933 |
| head 3-seed mean | 16,059 | 409 | 0.5898 | [0.5660, 0.6127] | 0.5955 |

Paired family-block deltas (same resamples, same families):

| comparison | delta macro F1 | 95% CI | p (two-sided) | families better |
|---|---:|---|---:|---:|
| head mean - RNAfold | -0.0106 | [-0.0312, +0.0095] | 0.3092 | 0.406 |
| head mean - EternaFold | -0.0152 | [-0.0350, +0.0046] | 0.1256 | 0.386 |
| head mean - RNAstructure | -0.0003 | [-0.0204, +0.0196] | 0.9840 | 0.430 |
| head seed17 - RNAfold | +0.0050 | [-0.0171, +0.0272] | 0.6580 | 0.443 |
| head seed17 - EternaFold | +0.0004 | [-0.0214, +0.0224] | 0.9726 | 0.433 |

### id_test

| method | n | families | family-macro F1 | 95% CI | micro F1 |
|---|---:|---:|---:|---|---:|
| RNAfold | 1,287 | 155 | 0.5955 | [0.5520, 0.6387] | 0.5295 |
| EternaFold | 1,287 | 155 | 0.5916 | [0.5479, 0.6338] | 0.5465 |
| RNAstructure | 1,287 | 155 | 0.5868 | [0.5430, 0.6302] | 0.5169 |
| head seed17 | 1,287 | 155 | 0.6009 | [0.5594, 0.6407] | 0.6057 |
| head seed23 | 1,287 | 155 | 0.5784 | [0.5359, 0.6194] | 0.5839 |
| head seed41 | 1,287 | 155 | 0.5772 | [0.5366, 0.6172] | 0.5872 |
| head 3-seed mean | 1,287 | 155 | 0.5855 | [0.5453, 0.6243] | 0.5920 |

| comparison | delta macro F1 | 95% CI | p (two-sided) | families better |
|---|---:|---|---:|---:|
| head mean - RNAfold | -0.0101 | [-0.0391, +0.0180] | 0.4884 | 0.439 |
| head mean - EternaFold | -0.0061 | [-0.0311, +0.0185] | 0.6218 | 0.432 |
| head mean - RNAstructure | -0.0013 | [-0.0296, +0.0264] | 0.9332 | 0.458 |
| head seed17 - RNAfold | +0.0053 | [-0.0273, +0.0375] | 0.7518 | 0.503 |
| head seed17 - EternaFold | +0.0092 | [-0.0229, +0.0409] | 0.5708 | 0.510 |

### Reading

- The seed spread (withheld 0.5785-0.6054, id_test 0.5772-0.6009) is as large
  as any head-baseline gap. Only **paired** family-block intervals support a
  comparative claim, and every one of them includes zero.
- The 3-seed mean is below all three thermodynamic tools on both splits. Seed17
  alone is nominally above EternaFold/RNAfold by +0.004 to +0.009, with CIs
  spanning zero and roughly 44-51% of families better - indistinguishable from
  noise.
- Micro F1 is higher for the heads than for the tools (0.5955 vs 0.5560/0.5908
  withheld). That is a length/family-weighting effect and cannot be quoted as a
  headline without the macro column beside it.
- Verdict: **no statistically supported advantage over any external baseline.**
  Any statement that the head "beats thermodynamics" is unsupported and must not
  appear in a draft.
- E-C2 gate target stays the higher baseline per split: EternaFold withheld
  0.6050, RNAfold id_test 0.5955.

`head 3-seed mean` is the arithmetic mean of the per-family F1 of three separately
trained heads, not a single deployable model; its micro is computed from the
summed tp/fp/fn of the three runs. Treat it as a variance-reduction summary only.

## 6. Freezing procedure (to run when 15631 finishes)

```
python scripts/build_common_cohort.py \
  --split splits/s1prime_filtered/<split>.jsonl \
  --dump experiments/002_e1prime_runs/dump_p2fix_seed17/a_<split>.jsonl \
  --preds experiments/002_e1prime_runs/{rnafold,eternafold,rnastructure,mxfold2}_<split>.jsonl \
  --out experiments/002_e1prime_runs/cohort/<split>.jsonl \
  --report experiments/002_e1prime_runs/cohort/<split>.report.json \
  --excluded-out experiments/002_e1prime_runs/cohort/<split>.excluded.jsonl
```

Then rescore every baseline and every head with `--restrict-ids`, and rerun
`scripts/cohort_paired_eval.py` on the new cohort. Adding MXfold2 shrinks the
cohort (600 nt cap) and every number in section 5 will change; the current table
must then be marked superseded, not silently reused. The run is accepted only if
`n_missing_prediction == 0` for every tool and the exclusion accounting is
identical across rows.

## 7. Limits of this table

- RNA is at most 1,014 nt (withheld) / 990 nt (id_test). Two withheld families
  are absent entirely. Conclusions do not extend to long RNA, where RNAstructure
  and MXfold2 are comparatively strongest.
- Family-macro F1 weights 409 (withheld) / 155 (id_test) families equally, but
  only families with at least one eligible sequence are represented. Small
  families with 1-2 members carry the same weight as large ones, so a single
  hard family moves the macro number.
- The nested projection is applied to both truth and prediction, so crossing
  pairs are not scored here.
- It cannot be compared to any older table in this repository that used the
  pre-fix `pair_prf` convention or a per-tool population. Those numbers are
  marked as such in `docs/corrections.md`.

## 8. Artifacts

- `experiments/002_e1prime_runs/cohort/{id_test,test_withheld}.jsonl`
- `experiments/002_e1prime_runs/cohort/{id_test,test_withheld}.report.json`
- `experiments/002_e1prime_runs/cohort/paired_{id_test,test_withheld}.json`
- `experiments/002_e1prime_runs/{rnafold,eternafold,rnastructure}_*_score.json`
- `experiments/002_e1prime_runs/{rnafold,eternafold,rnastructure}_*_by_family.json`
- `experiments/002_e1prime_runs/coverage/{withheld,id_test}.json`
- `scripts/build_common_cohort.py`, `scripts/score_external_baseline.py`,
  `scripts/cohort_paired_eval.py`, `scripts/export_head_preds.py`,
  `scripts/predict_rnastructure.py`, `scripts/predict_eternafold.py`,
  `scripts/export_mxfold2_bpseq.py`
- `docs/rnafold_baseline_2026-11-13.md`
