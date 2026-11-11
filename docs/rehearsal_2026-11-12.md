# Confirmatory-path rehearsal, 2026-11-12 (non-protocol)

PI Round 4g §3 asked for the confirmatory chain to be exercised end to end on
real sequences before the cluster returns, with stand-in scores. This is a
**pipeline test**. Nothing here is a result, and nothing here may be quoted as
one.

## Stand-in scores

`scripts/make_standin_dump.py` sampled 300 local ArchiveII sequences
(`max-length 400`, `seed 17`) and wrote the dump schema D1/H2.4 consume. Scores
are ViennaRNA 2.8.0 partition-function base-pair probabilities in logit form.

This is circular: ViennaRNA probabilities fed back into a ViennaRNA soft-constrained
fold. EternaFold and CONTRAfold were the preferred stand-ins and are not
installed locally. The circularity makes the energy arm *harder*, not easier,
because the pseudo-energies and the Turner model then agree by construction.

## Chain, in order

1. `scripts/fit_energy_calibration.py`, 5-fold, A7.4 coarse grid
   (`s ∈ {0.25,0.5,1,2,4}`, `b ∈ {-6,-3,0,3,6}`):
   selected `s = 0.5`, `b = 6.0`, mean fold F1 0.5817.
   **`b` landed on the edge of the coarse grid**, so the A7.4 refinement step is
   not optional: on real data the refinement must run, and if the refined optimum
   is also on an edge that has to be reported rather than silently accepted.
2. `scripts/h24_decoder_compare.py --energy struct` (A-struct, zero penalties):
   basis `recall_matched`, status `confirmed`.
   8,092 true positives for **both** arms, 12,666 predicted pairs each, one
   fewer long-range false positive for the DP (678 vs 679).
   The arms are effectively identical on these scores.
3. `scripts/h24_decoder_compare.py --energy vienna --energy-scale 0.5
   --energy-offset 6.0` (A-energy-fixed): basis `null`, status
   `inconclusive_no_recall_match`, `pair_budget_feasible: false`
   (pair-count gap 1,523 against a tolerance of 160).
   The verdict machinery refused to produce a verdict rather than inventing one,
   which is the behaviour it was written for.
4. Embargo check: the same invocation with a `test_novel` path and no
   `--confirmatory` exits 1 with `EmbargoViolation` before reading anything.

## What this tells us

- The chain runs end to end and writes machine-readable verdicts.
- The verdict logic fails closed: no recall match and no pair-budget match gives
  `null`, not a default.
- **A-struct adds nothing on physically consistent probability maps.** Greedy
  nested matching and the max-weight nested DP return the same structure here.
  That is the opposite of the H2 premise, and it is a rehearsal warning: the
  H2.4 contrast may only bite on over-confident or adversarially inconsistent
  maps. This is not evidence either way — the stand-in is circular — but it says
  the real dumps must be examined for whether greedy is already optimal.
- The strong attraction implied by `b = 6` produces a far larger pair budget than
  greedy, so `(s, b)` at the coarse edge is not a usable operating point.

## Cost measurement

`scripts/grid_timing.py` (`experiments/009_grid_timing/grid_timing.json`), 8
sequences per length bin, production `ViennaSoftDecoder` which calls `sc_add_bp`
once per candidate cell:

| length | median s/fold | median constraints |
|---|---:|---:|
| 0–79 | 0.0032 | 969 |
| 80–149 | 0.0098 | 2,610 |
| 150–299 | 0.0497 | 11,422 |
| 300–599 | 0.1153 | 24,185 |
| 600–1999 | 0.4649 | 81,352 |

Projected over 17,700 training sequences with the A7.4 50-point grid: **≈12.4
CPU-hours** on this length mix — below the 24-hour subsample threshold, but
ArchiveII is short-biased. The decision must be recomputed on the real split's
length distribution before the fit, and A7.4 allows the 3,000-sequence subsample
(seed 17) if that projection crosses 24 CPU-hours.

## Files

- `experiments/010_rehearsal/dump.jsonl`, `split.jsonl`, `standin_meta.json`
- `experiments/010_rehearsal/fit_sb.json`
- `experiments/010_rehearsal/h24_struct.json`, `h24_energy.json`
- `experiments/009_grid_timing/grid_timing.json`
