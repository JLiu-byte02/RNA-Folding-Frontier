# Amendment 15A - Pre-unblinding protocol for H4, shift refinement, refit and five seeds

Status: written 2026-11-15, before any H4 label score, shift result, new-seed
result, or default-weight MXfold2 score on the project cohort. At the time of
writing the only H4 quantities inspected are source counts and train-overlap
counts: 533 `bpRNA_PDB_*` records in the frozen bpRNA-1m mirror, of which 347
remain after the k=15 train-containment screen described below. No H4
difference-in-difference has been computed.

This amendment is the protocol for the work PI Round 4n authorized. It does
**not** authorize reading `test_novel`: `AMENDMENT_14.md` remains in force until
`AMENDMENT_15.md` replaces its list. The earlier amendments remain immutable.

## A15A.1 - H4: annotation convention versus structural truth

### Population

The source is the `bpRNA_PDB_*` subset of the frozen
`data/raw/bprna/bpRNA_1m.json` mirror. The source object has 533 records and
only the fields `sequence` and `structure`.

The population is constructed once, before any predictor is run:

1. Normalise sequence to uppercase U.
2. Drop a record when its exact sequence occurs in
   `splits/s1prime_filtered/train.jsonl`.
3. Build unique 15-mer sets for every train sequence. Drop a record when its
   raw 15-mer containment against the train set is `>= .5`. A df-filtered
   containment (the project's 5% df index) is reported as a descriptive
   sensitivity, not used for the primary exclusion.
4. Run `cmscan --cut_ga --rfam --nohmmonly` against Rfam 15.0. Keep the
   highest-bitscore family hit per chain; break an exact score tie by the
   lexicographically smaller family accession. `NOHIT` chains are dropped.
5. Keep only chains whose selected family CM is available in
   `data/raw/rfam/Rfam.cm`.

No newer PDB chains are added in this round; that would require a separate
frozen source cut. Every dropped row and reason is written to the H4 cohort
manifest.

### Two labels

For each retained chain, build both labels under the same policy:

- **Comparative label:** `cmalign` the chain to the selected family CM, read
  the `#=GC SS_cons` row, parse it in alignment-column coordinates, map both
  endpoints to the ungapped chain with
  `rff.data.wuss.project_pairs`, and drop any pair with a gap at either
  endpoint.
- **Structural label:** the `structure` field from the frozen bpRNA PDB
  record. This is the fixed structural annotator; no second annotator is
  introduced.

For both labels: keep canonical AU/UA/GC/CG/GU/UG pairs only, take the maximum
nested subset, and require at least four pairs. A chain is in the H4 evaluation
set only when **both** labels survive. A missing or unparsable label is
reported as a construction failure, not replaced by the other label.

### Predictors and statistic

Predictors are:

- `p2fix_seed{17,23,41}` with the deployed threshold-.5 greedy/nested decoder;
- `onehot_seed{17,23,41}` with the identical decoder;
- RNAfold and EternaFold with their frozen baseline settings.

Every predictor must have a prediction for every H4 chain; otherwise the chain
is dropped from the common H4 set and counted. No subset comparison is allowed.

For label `l` and predictor `m`, compute `F_m,l(f)` as the mean
family-macro F1 over the eligible sequences in family `f`; learned predictors
are averaged over their seeds inside each family. Define

```
D(f) = [F_head,comparative(f) - F_RNAfold,comparative(f)]
     - [F_head,structural(f)  - F_RNAfold,structural(f)]
D    = mean_f D(f)
```

The 95% interval is a 10,000-resample family bootstrap, seed `20261115`.
The one-hot control uses the same definition with the one-hot seed mean.
Report both label-specific tables, per-seed values, exact family counts, and
all exclusion reasons.

### Pre-registered H4 decision

- **Supported:** `D >= .03`, its family-bootstrap 95% CI excludes zero, and
  the one-hot control has `D > 0` in the same direction.
- **Directionally positive but not supported:** `.01 <= D < .03` with a CI
  excluding zero. Report it as insufficient for the H4 claim; no extension.
- **Killed:** `D < .01`, or the family-bootstrap CI includes zero.
- **Underpowered:** fewer than 40 families have both labels. Report the family
  count and stop; do not add sequences or relax the label rule.

If the head passes but the one-hot control is non-positive, H4 is not
supported; report that pattern as an FM-specific association without calling
it annotation convention.

## A15A.2 - Shift-corrected decoding (secondary analysis)

### Inputs and greedy baseline

Use the frozen common cohort (withheld 16,040 / 408; id_test 1,272 / 153),
the existing `headpreds/seed{17,23,41}_<split>.jsonl` greedy predictions, and
the matching p2fix probability dumps. A probability `p` is converted to
`logit(p)` with `p` clipped to `[1e-6, 1-1e-6]`.

The greedy baseline is the deployed threshold-.5 nested greedy decoder, not a
thresholded probability map.

### Frozen refinement rule

For each seed and sequence:

1. Partition accepted pairs into maximal stacked helices: consecutive
   `(i,j), (i+1,j-1), ...` pairs with no missing rung.
2. Order helices by decreasing summed logit of their original pairs; break ties
   by smaller outer `i`, then smaller outer `j`.
3. Remove the current helix from the accepted set. Try
   `delta in {-2,-1,+1,+2}` on its 3' strand:
   `(i+k, j-k) -> (i+k, j-k+delta)` for every rung.
4. Keep a candidate only if every rung is in bounds, canonical/GU, has no
   shared endpoint with any accepted pair, and does not cross any accepted
   pair. Pre-existing pairs in the helix are excluded from that comparison.
5. Accept the valid candidate with the largest strictly positive increase in
   the helix's summed logit. Ties use the smaller `abs(delta)`, then the order
   `-2,-1,+1,+2`. Otherwise restore the original helix.
6. Do one pass; do not re-process a helix after later shifts.

The secondary exploratory rule is the same geometry with the helix's Turner
stacking energy as the criterion. It is reported but cannot define the gate.

### Gate and reporting

The primary statistic is the three-seed mean paired difference
`shift - greedy` in family-macro F1. The family bootstrap has 10,000 resamples
and seed `20261115`.

- **Pass:** withheld mean gain `>= .008` with CI excluding zero, and id_test
  mean point estimate `>= 0`.
- **Kill:** mean gain `< .004`, or the withheld CI includes zero.
- If `.004 <= gain < .008` with a CI excluding zero, report "below the
  pre-registered pass bar"; it does not enter `test_novel`.

Report every seed, the fraction of the radius-2 oracle gain recovered
separately for withheld/id_test, and changes in near-miss (1-2 nt) and far
(>10 nt) false-positive counts. No threshold or tie-break change is allowed
after seeing the output.

## A15A.3 - Clan-grouped `(s,b)` refit

The fit uses the existing frozen train dump
`experiments/002_e1prime_runs/dump_p2fix_seed17/train.jsonl` and
`splits/s1prime_filtered/train.jsonl`, with v3 nested truth and the `>=4`-pair
exclusion already used by the fitter.

The code change is only the fold assignment:

- group records by their frozen `clan` field;
- assign each clan to `crc32(clan) % 5`;
- never split a clan across folds; if a record lacks a clan, use its family as
  the group and report the count;
- use the same folds for the coarse grid, refinement and any boundary
  extension.

The grids are exactly Amendment 7 A7.4: coarse
`s={.25,.5,1,2,4}` x `b={-6,-3,0,3,6}`, then the second 5x5 grid spanning
one coarse step around the coarse optimum. If the refined winner is on the
minimum or maximum of either axis, extend that axis one step outward once,
using the same family of steps, and rerun the grid once. Do not extend
recursively. Ties select the smaller `s`, then smaller `b`.

Freeze the selected `(s,b)` JSON and its hash before scoring id_test or
withheld. Report A-fixed on withheld under the old and refit values; only the
refit value goes to `test_novel`.

## A15A.4 - Two additional head seeds

Seeds **43 and 47** are pre-declared, giving five seeds total. They use the
existing frozen train split, the same `giga_a_train` cache, the same
`scripts/train_head_s1.py` recipe and the same corrected model-selection
metric. No architecture, learning-rate, validation-family split, or early
stopping change is permitted.

The five-seed table reports:

- the five per-seed family-macro values;
- their mean and the across-seed Student-t interval with four degrees of
  freedom;
- the family-block bootstrap interval for the seed mean;
- the same paired head-minus-baseline deltas as the current table.

No seed may be dropped after seeing a withheld or id_test result.

## A15A.5 - Default-weight MXfold2 row

Run the released `TrainSetAB.pth` through the existing
`scripts/predict_mxfold2.py` wrapper with model `MixAB`, default released
weights, CPU folding and the same 600 nt cap. The common cohort is already
below that cap. This is one descriptive row only:

- flag its training-data overlap with our withheld and id_test sequences;
- it cannot define the best clean baseline;
- the short-schedule retrain remains labelled as such;
- no further MXfold2 training is authorized.

## A15A.6 - Order, provenance and `test_novel`

Order: A15A pin -> H4 -> shift refinement -> refit -> seeds 43/47 -> default
MXfold2 -> `AMENDMENT_15.md` -> one `test_novel` pass.

Every output directory records the input hashes (train split, raw bpRNA,
Rfam seed/CM, label files), the code commit, and the output hashes. Any code
defect after a result is seen is written to `docs/corrections.md`; affected
items are rerun and the rerun is disclosed.

A15A does not authorize a `test_novel` read. Amendment 15 will pin the final
analysis list after these results, explicitly supersede Amendment 14 without
editing it, and only then may `test_novel` be read once.

## A15A.7 - Not authorized

- No new architecture, hyperparameter sweep, alternate H4 label source, or
  additional seed after a result is seen.
- No threshold relaxation for the shift rule.
- No use of `test_novel` for any model, label, calibration, seed or tie-break
  choice.
- No claim that D establishes a general annotation-convention mechanism if
  the family count or one-hot control fails the stated gate.
