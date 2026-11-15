# False-positive anatomy (exploratory, 2026-11-15)

**Question.** When the head is wrong, is it wrong in the same way as RNAfold?
For every false-positive pair emitted at threshold 0.5, measure the Manhattan
distance in pair coordinates to the nearest true pair, `d(i,j) = min |i-i'| +
|j-j'|`, and bin it. Near-miss FPs (d = 1-2) mean the helix is placed but
offset; far FPs (d > 10) mean structure invented elsewhere.

**Data.** The frozen four-method cohort, truth = nested v3 pairs, condition a
(native context), head seed 17; baseline = RNAfold predictions on the same ids.

| split | method | FP pairs | near-miss (1-2 nt) | far (>10 nt) | mean d |
|---|---|---:|---:|---:|---:|
| withheld | head | 359,939 | **.199** | **.570** | 35.7 |
| withheld | RNAfold | 311,768 | .134 | .643 | 36.1 |
| id_test | head | 38,583 | **.295** | **.424** | 23.0 |
| id_test | RNAfold | 30,469 | .203 | .513 | 26.1 |

**Reading.** The head emits more false positives overall, but a larger share of
them are near-misses (+6.5 pp on withheld, +9.3 pp on id_test) and a smaller
share are far inventions (-7.3 pp, -8.9 pp). The pattern is consistent across
splits: the head tends to place helices near the right location with
endpoint/register errors, whereas RNAfold's errors are more often elsewhere.

**Caveat.** The two methods do not emit the same number of FPs, so this compares
distributions, not matched counts. Distance in pair coordinates is a proxy for
"same helix, wrong register" and is not causal.

**Possible next test (not run).** A shift-corrected decoder: for each predicted
helix, try register shifts of +/-1 and +/-2 and keep a shift only when the model
score improves; evaluate on the frozen cohort. If near-miss FPs convert to TPs,
the head's macro deficit may be partly a decoding artefact rather than a
representation limit.

**Artifacts.** `experiments/002_e1prime_runs/cohort/fp_anatomy_{test_withheld,id_test}.json`,
produced by `scripts/fp_anatomy.py`.
