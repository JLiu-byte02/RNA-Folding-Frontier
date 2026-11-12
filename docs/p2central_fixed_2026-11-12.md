# Corrected P2-central after the scorer fix — 2026-11-12

This is the table required by Round 4j item 3. It uses the retrained condition-a
heads `runs/p2fix_seed{17,23,41}/best.pt`, v3 labels, v3 nested primary truth,
and fresh probability dumps for conditions a–e. No H2.4 or P2-block run was
launched.

## Protocol

- Three seeds: 17, 23, 41.
- Splits: id_test and test_withheld.
- Conditions: a native, b dinucleotide shuffle, c same-family homolog below
  60% identity, d other-clan donor, e masked/N distal.
- Decoding: threshold + one partner + non-crossing.
- Fixed threshold: 0.5.
- Oracle threshold: chosen per condition on the 80–299 nt band by
  family-macro all-pairs F1, on the same 15-point grid used by D1.
- Readouts: kept–kept pairs only (both endpoints in the record's kept window),
  all pairs, and pooled micro F1.
- Rows with fewer than four v3 true pairs are excluded. For kept–kept, rows with
  zero truth pairs in the kept window are excluded and counted.

## Results (three-seed mean)

| split | cond | fixed kept | fixed all | fixed micro | oracle cut | oracle kept | oracle all | oracle micro |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| id_test | a | .5759 | .5837 | .5937 | .75* | .5753 | .5964 | .6161 |
| id_test | b | .5203 | .3335 | .2132 | .78* | .5186 | .3425 | .2393 |
| id_test | c | .5616 | .4380 | .3467 | .85* | .5755 | .4541 | .3718 |
| id_test | d | .5035 | .3224 | .2110 | .80* | .5118 | .3317 | .2364 |
| id_test | e | .5142 | .3488 | .2677 | .60* | .5216 | .3517 | .2710 |
| withheld | a | .5910 | .5911 | .5942 | .72* | .5906 | .6019 | .6089 |
| withheld | b | .5275 | .3478 | .2743 | .77* | .5286 | .3565 | .2982 |
| withheld | c | .5927 | .4927 | .4044 | .80* | .6117 | .5141 | .4223 |
| withheld | d | .5231 | .3463 | .2670 | .73* | .5271 | .3535 | .2856 |
| withheld | e | .5171 | .3620 | .3259 | .58* | .5236 | .3647 | .3293 |

\* The per-seed oracle cuts differ; the table reports the seed mean of the
metrics, not a metric computed at the mean cut. The raw per-seed cuts and all
counts are in `experiments/002_e1prime_runs/p2central_fixed_seed*.json`.

## Immediate reading

- **a is the best condition overall**, on both splits and both thresholds.
- **b, d and e lose most of the all-pairs F1**, and a substantial part of the
  kept-window F1. The all-pairs loss is larger than the kept-window loss, as
  expected when distal content is destroyed.
- **c retains much more of a's kept-window performance**, especially on
  withheld: kept–kept .5927 fixed / .6117 oracle versus a .5910/.5906. On
  id_test c is .5616/.5755 versus a .5759/.5753. The c–d contrast is +.058
  (id_test fixed kept) and +.070 (withheld fixed kept).
- This c–d contrast is consistent with a family/template shortcut, but it is
  not yet evidence of one: the positive control and P2-block are explicitly not
  run, per Round 4j.
- e is worse than a on every readout. This is the sanity check required after
  the scorer fix, not a new finding.

### Population correction for the retention ratios

The c–d contrast above is computed from each condition's own row set. Condition
c has a verified homolog only on a subset (449 id_test / 7,572 withheld after
the ≥4-pair rule), so a ratio that divides c's F1 by a and e computed on the
full split mixes populations. On the common a/c/e population the corrected
family-bootstrap estimates are:

| split | readout | R_c | R_d |
|---|---|---:|---:|
| id_test | kept–kept | .612 [.413, .853] | −.173 [−.617, .080] |
| withheld | kept–kept | .747 [.564, .960] | .082 [−.040, .203] |
| withheld | all-pairs | .504 [.390, .606] | −.069 [−.103, −.036] |

So P2-central does **not** meet the H3 90% retention bar. H3 is decided on
P2-block, where all conditions are measured on the same records. Details:
`docs/p2central_ratios_ci_2026-11-12.md`.

## Mandatory invariants

For every condition, the family-macro and pooled-micro values differ by less
than .15; the largest fixed kept-window gap is .12 (id_test b). The corrected
dinucleotide-shuffle floor is .0162 (id_test) / .0173 (withheld), below the .05
failure threshold. Both reporting invariants pass.

## Artifacts

- `experiments/002_e1prime_runs/p2central_fixed_seed17.json`
- `experiments/002_e1prime_runs/p2central_fixed_seed23.json`
- `experiments/002_e1prime_runs/p2central_fixed_seed41.json`
- `experiments/002_e1prime_runs/p2central_fixed_summary.json`
- `experiments/002_e1prime_runs/p2central_ratios_ci.json` (family-bootstrap CIs
  for the c/d retention ratios) and `docs/p2central_ratios_ci_2026-11-12.md`
- `scripts/p2central_eval.py`
- `scripts/aggregate_p2central.py`

`test_novel` remains unscored. No result here is a final test-set claim.
