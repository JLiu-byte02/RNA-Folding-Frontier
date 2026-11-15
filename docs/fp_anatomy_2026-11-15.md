# False-positive anatomy and shift-oracle bound (exploratory, 2026-11-15)

**Question.** When the head is wrong, is it wrong in the same way as RNAfold?
For every false-positive pair emitted by the deployed decoder, measure the
Manhattan distance in pair coordinates to the nearest true pair,
`d(i,j) = min |i-i'| + |j-j'|`, and bin it. Near-miss FPs (d = 1-2) mean the
helix is placed but offset; far FPs (d > 10) mean structure is invented
elsewhere.

**Data.** The frozen four-method cohort, truth = nested v3 pairs, condition a
(native context), head seed 17, and the deployed greedy head predictions
(`headpreds/seed17_{test_withheld,id_test}.jsonl`), not a thresholded
probability map. Baseline = RNAfold predictions on the same ids.

| split | method | FP pairs | near-miss (1-2 nt) | far (>10 nt) | mean d |
|---|---|---:|---:|---:|---:|
| withheld | head | 188,158 | **.2193** | **.5164** | 32.679 |
| withheld | RNAfold | 311,768 | .1337 | .6431 | 36.095 |
| id_test | head | 19,935 | **.3107** | **.3494** | 18.727 |
| id_test | RNAfold | 30,469 | .2025 | .5133 | 26.145 |

**Reading.** The deployed head emits *fewer* false positives than RNAfold on
both splits, but a larger share of its FPs are near-misses (+8.6 pp on
withheld, +10.8 pp on id_test) and a smaller share are far inventions
(-12.7 pp, -16.4 pp). The pattern is consistent across splits: the head more
often places a helix near the right location with endpoint/register errors,
whereas RNAfold's errors are more often elsewhere. This is a distributional
comparison, not a matched-count comparison; distance in pair coordinates is a
proxy for "same helix, wrong register" and is not causal.

**Shift-oracle upper bound.** For each sequence, count an unpredicted truth pair
as recoverable when at least one head FP lies within Manhattan radius 2, then
pretend an oracle snaps the FP onto that truth pair. This is an upper bound: it
ignores crossing constraints, competing FPs, shift-selection errors, and the
decoder's model score.

| split | base macro | oracle macro | delta | shiftable FPs | share of head FPs |
|---|---:|---:|---:|---:|---:|
| withheld | .6061 | .6312 | **+.0251** | 13,155 | 7.0% |
| id_test | .6057 | .6319 | **+.0263** | 2,621 | 13.2% |

The base macro-F1 values match the frozen head-17 entries in
`paired_{test_withheld,id_test}.json` (.6061/.6057), so the diagnostic is using
the same cohort, labels, and decoder as the main table. The oracle does not
establish that a shift-corrected decoder can recover this gain; it only bounds
the gain that is available before structural constraints and selection costs.

**Correction.** The first pass used thresholded probability-map cells as if
they were head predictions. That pass is withdrawn. The corrected artifacts use
the deployed prediction files and are recorded in `docs/corrections.md`.

**Possible next test (not run).** A shift-corrected decoder: for each predicted
helix, try register shifts of +/-1 and +/-2 and keep a shift only when the model
score improves; evaluate on the frozen cohort. If near-miss FPs convert to TPs,
the head's macro deficit may be partly a decoding artefact rather than a
representation limit. The PI must pre-register the success bar and kill
criterion before this test; the oracle is not a result.

**Artifacts.** `experiments/002_e1prime_runs/cohort/fp_anatomy_{test_withheld,id_test}.json`
and `experiments/002_e1prime_runs/cohort/shift_oracle_{test_withheld,id_test}.json`,
produced by `scripts/fp_anatomy.py` and `scripts/shift_oracle.py`.

**Reproduction (server checkout).**

```bash
for split in test_withheld id_test; do
  python scripts/fp_anatomy.py \
    --split splits/s1prime_filtered/${split}.jsonl \
    --cohort experiments/002_e1prime_runs/cohort/${split}.jsonl \
    --head-preds experiments/002_e1prime_runs/headpreds/seed17_${split}.jsonl \
    --preds experiments/002_e1prime_runs/rnafold_${split}.jsonl \
    --out experiments/002_e1prime_runs/cohort/fp_anatomy_${split}.json
  python scripts/shift_oracle.py \
    --split splits/s1prime_filtered/${split}.jsonl \
    --cohort experiments/002_e1prime_runs/cohort/${split}.jsonl \
    --head-preds experiments/002_e1prime_runs/headpreds/seed17_${split}.jsonl \
    --out experiments/002_e1prime_runs/cohort/shift_oracle_${split}.json
done
```
