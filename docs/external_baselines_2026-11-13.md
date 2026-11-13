# External baselines on the v3 protocol and the common-cohort policy — 2026-11-13

Status: measurement only. RNAfold and EternaFold are complete; RNAstructure
(job 15649) and MXfold2 (job 15631, retrained on our train split) are still
running, so the final common cohort is not yet frozen. Nothing in this file is
a new method.

## 1. The defect this policy fixes

`scripts/score_external_baseline.py` skips a sequence when a prediction is
missing. Each tool has its own practical limit, so each baseline was scored on a
different population:

- EternaFold: 1,200 nt cap (CONTRAfold is O(L^3)); 578 withheld rows dropped;
- RNAstructure: 111 rows with IUPAC ambiguity codes are rejected outright by
  `Fold` ("Invalid nucleobase") and skipped by declared rule;
- MXfold2: trained and folded with a 600 nt cap;
- our own heads: a 1,024 nt context limit.

A family-macro column computed over different populations is not comparable.
The fix is one explicit id list that every method — our heads included — is
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
- This is the corrected statement. Earlier drafts of the id_test coverage note
  mixed the two loss modes; the split has 168 families and the coverage loss is
  zero.

v3 label eligibility inside the dump:

| split | empty truth | `<4` pairs | `>=4` pairs (scored) | pair p25 / median / p75 |
|---|---:|---:|---:|---|
| test_withheld | 118 | 630 | 16,163 | 14 / 24 / 34 |
| id_test | 5 | 105 | 1,287 | 16 / 24 / 32 |

Tool-cap overlap with the 1,024 nt gap:

- withheld: 578 rows are longer than 1,024 nt and 578 are longer than 1,200 nt.
  The sets coincide exactly, so EternaFold's cap removes nothing our context
  did not already remove.
- id_test: 10 rows are longer than 1,024 nt but only 7 are longer than 1,200 nt.
  EternaFold folds the three rows in (1,024, 1,200] and they are then dropped by
  the restrict set. Our context is still the binding limit; EternaFold merely
  fails to cover 7 of the 10.

## 4. Per-tool exclusions

| tool | rule | withheld rows | id_test rows | where recorded |
|---|---|---:|---:|---|
| RNAfold 2.8.0 | none | 0 | 0 | `rnafold_*_score.json` |
| EternaFold (CONTRAfold-SE) | `--max-length 1200` | 578 | 7 (+3 dropped by restrict set) | job 15653 log, `n_outside_restrict_set` |
| RNAstructure 6.6 | `--skip-ambiguous` (IUPAC) + `--allow-partial` timeouts | 111 (pending job) | pending | `<out>.skipped.jsonl`, `<out>.errors.jsonl` |
| MXfold2 (MixC, retrained) | 600 nt train cap; test cap set at export | pending | pending | `experiments/013_mxfold2/run/export_manifest.json` |
| P2-fix heads (ours) | 1,024 nt context | 578 | 10 | `coverage/{withheld,id_test}.json` |

RNAstructure's 111 ambiguous rows are a *tool* limitation, not a length limit.
They are inside our dump, so they will shrink the final cohort by up to 111 rows
and some families, which is exactly why the cohort must be built after every
prediction file exists.

## 5. Completed results (before the common cohort is frozen)

Scored on the P2-fix `a` dump ids (id_test 1,287; withheld 16,163), v3 labels,
`--min-true-pairs 4`, nested truth, 10,000 family-block bootstrap, seed
20261113.

| split | tool | n seqs | families | family-macro F1 | 95% CI | micro F1 | mean pred pairs |
|---|---|---:|---:|---:|---|---:|---:|
| withheld | RNAfold | 16,163 | 409 | **0.6003** | [0.5733, 0.6265] | 0.5548 | 37.24 |
| withheld | EternaFold | 16,163 | 409 | **0.6051** | [0.5792, 0.6302] | 0.5901 | 36.24 |
| id_test | RNAfold | 1,287 | 155 | **0.5955** | [0.5530, 0.6375] | 0.5295 | 46.97 |
| id_test | EternaFold | 1,287 | 155 | **0.5916** | [0.5479, 0.6329] | 0.5465 | 46.14 |

For reference, the P2-fix `a` head on the same population (3-seed mean of the
fixed-0.5 all-pairs readout):

| split | head macro | RNAfold macro | d(macro) | head micro | RNAfold micro | d(micro) |
|---|---:|---:|---:|---:|---:|---:|
| withheld | 0.5911 | 0.6003 | -0.0092 | 0.5942 | 0.5548 | +0.0394 |
| id_test | 0.5837 | 0.5955 | -0.0118 | 0.5937 | 0.5295 | +0.0642 |

The head is below plain RNAfold in family-macro F1 and above it in micro F1 on
both splits. The reading is that the head wins on large families (which carry
micro weight) and loses on the many small ones (which the macro weights
equally). This is consistent with E-C1: the per-sequence FM advantage over
RNAfold is small, and it does not survive family re-weighting.

Any statement that the head "beats thermodynamics" is unsupported. The
E-C2 gate target is therefore the higher of the two baselines, not the head's own
level: EternaFold withheld 0.6051, id_test RNAfold 0.5955.

## 6. Freezing procedure (to run when 15631 and 15649 finish)

```
python scripts/build_common_cohort.py   --split splits/s1prime_filtered/<split>.jsonl   --dump experiments/002_e1prime_runs/dump_p2fix_seed17/a_<split>.jsonl   --preds experiments/002_e1prime_runs/{rnafold,eternafold,rnastructure,mxfold2}_<split>.jsonl   --out experiments/002_e1prime_runs/cohort/<split>.jsonl   --report experiments/002_e1prime_runs/cohort/<split>.report.json   --excluded-out experiments/002_e1prime_runs/cohort/<split>.excluded.jsonl
```

Then rescore every baseline and every head with `--restrict-ids`. The run is
accepted only if `n_missing_prediction == 0` for every tool and the exclusion
accounting is identical across rows; otherwise the cohort is wrong and the table
is not printed.

## 7. Limits of this table

- It is restricted to RNA of at most 1,014 nt (withheld) / 990 nt (id_test).
  Two withheld families are absent entirely. Conclusions do not extend to long
  RNA, which is where RNAstructure and MXfold2 are comparatively strongest.
- Family-macro F1 weights 409 (withheld) / 155 (id_test) families equally, but
  only families with at least one eligible sequence are represented.
- It cannot be compared to any older table in this repository that used the
  pre-fix `pair_prf` convention or a per-tool population. Those numbers are
  marked as such in `docs/corrections.md`.

## 8. Artifacts

- `experiments/002_e1prime_runs/{rnafold,eternafold}_*_score.json`
- `experiments/002_e1prime_runs/{rnafold,eternafold}_*_by_family.json`
- `experiments/002_e1prime_runs/coverage/{withheld,id_test}.json`
- `scripts/build_common_cohort.py`, `scripts/score_external_baseline.py`,
  `scripts/predict_rnastructure.py`, `scripts/predict_eternafold.py`,
  `scripts/export_mxfold2_bpseq.py`
- `docs/rnafold_baseline_2026-11-13.md`
