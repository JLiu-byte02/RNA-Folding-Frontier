# Amendment 10 — P2-block, P3 train-and-test, D1-core, and H3

Status: written 2026-11-12 **before any P2-block cache exists or has been
read**, and before any P3 head has been trained. It carries PI Round 4k
(`docs/claude/transcripts/2026-11-12_round04k_reply.md`). Amendments 8 and 9
remain in force for labels and scoring; this file only adds the P2-block/P3
decision rules, the H3 estimand, and the calibration order.

At the time of writing, `data/cache/` contains no `giga_block_*` directory and
no P3 checkpoint exists. The only numbers quoted below are already-frozen
central-window results; they are not used to choose a threshold, a block
geometry, or a seed.

## A10.1 — What P2-central did and did not establish

The corrected P2-central table is a **test-time central-window** intervention:
one condition-a-trained head, with the distal context of each sequence replaced
at evaluation. It is exploratory for D1 because a single forward pass cannot
give every position its own native flank. Its pattern is:

- condition c retains about 77% (id_test) and 102% (withheld) of the a−e
  kept–kept difference;
- condition d retains about −17% and 8%;
- a still beats e on withheld kept–kept pairs.

That pattern is consistent with the alternative hypothesis that evolutionarily
consistent distal context transfers, but it does not test whether a head
**trained** without distal context generalises better OOD. P3 does.

The P2-central ratios are exploratory. Bootstrap CIs for them are added to the
report, but no D1 verdict is taken from P2-central.

## A10.2 — P2-block construction (frozen geometry)

The implementation is `scripts/build_block_donors.py` and
`scripts/build_block_cache.py`; the geometry is fixed here and may not be
changed after a cache is read:

- non-overlapping blocks of **16 nt**;
- **32 nt** of native flank on each side of a block, clipped at the sequence
  ends;
- maximum length **1024 nt**;
- donor everywhere except the block plus its native flank, then keep only the
  block's positions; concatenating the blocks covers every position exactly
  once;
- condition **a**: the existing full-context cache (`giga_a_{split}`), which is
  identical to the block composite because the donor equals the native string;
- condition **b**: dinucleotide-preserving shuffle of the whole sequence, seed
  17;
- condition **c**: same-family homolog below 60% identity, projected onto the
  target; positions aligned to a gap keep the native base; only records with a
  verified homolog have a c donor;
- condition **d**: a same-length slice from a sequence of a different clan,
  drawn from the train pool with seed 17;
- condition **e**: `N` at every position outside the block plus its flank.

Condition e-block is therefore **approximately a k=32 crop**, not an idealised
distal-only ablation: each block's positions are embedded with 32 nt of their
own native flank, but the backbone also sees an N-padded sequence outside that
window. Condition b-block preserves local composition and length but removes
the native distal content. The D1-core verdict is about these interventions as
implemented, not about an unobservable counterfactual in which the FM sees no
distal context at all.

Every cache writes a `manifest.json` with `records_cached`, `forward_passes`,
`seconds_per_sequence`, and the geometry. A cache that does not match this
geometry is not used.

## A10.3 — Test-time arm and H3

The test-time arm uses the three corrected condition-a heads
`runs/p2fix_seed{17,23,41}/best.pt` (Amendment 9 §A9.4). The weights are fixed
across conditions: only the input cache changes. Conditions a–e are evaluated
on `id_test` and `test_withheld` at the fixed threshold **0.5**; no threshold
is tuned on either split.

Truth policy is Amendment 8: v3 canonical plus GU pairs, all-pairs for this
readout, and sequences with fewer than four v3 true pairs excluded. bpRNA rows
are unchanged. `pair_prf` follows Amendment 9 §A9.1 (zero true positives with a
non-empty prediction is F1 = 0). Every row reports family-macro F1 and pooled
micro F1 side by side; nested-subset primary F1 and kept–kept F1 are secondary
readouts.

Condition c has a restricted population. **H3 is evaluated on the common
subset of withheld records that have a c donor and at least four v3 true
pairs**, with a, c, and e all restricted to the same record ids. The full-split
a/b/d/e table is reported separately and is not substituted into H3.

For each seed and family, compute the family-mean all-pairs F1 at threshold 0.5
for each condition. Average each family's F1 across the three seeds, then take
the family macro-average. Define

```
R_c = (F1_c - F1_e) / (F1_a - F1_e)
R_d = (F1_d - F1_e) / (F1_a - F1_e)
```

with a 95% percentile bootstrap over families (10,000 resamples,
`numpy.random.default_rng(20261112)`), resampling the same families for all
conditions and seeds. If the point-estimate denominator is zero or negative,
H3 is reported as undefined and not supported; no sign is flipped. If more than
1% of bootstrap resamples have a non-positive denominator, the ratio interval is
reported as unstable and cannot support H3.

**H3 is supported only if the 95% CI lower bound for `R_c` is at least 0.90
and the 95% CI upper bound for `R_d` is at most 0.30.** Otherwise H3 is not
supported. There is no post-hoc relaxation of these bounds and no substitution
of a different metric. `id_test` is reported descriptively with the same
estimand; no H3 verdict is taken there.

## A10.4 — Exploratory dose-response for condition c

On the same common-c withheld subset, split records by the homolog identity
stored in `data/processed/homologs/test_withheld.json`:

- bin 1: `0.40 <= identity < 0.50`;
- bin 2: `0.50 <= identity < 0.60`.

Compute `R_c` within each bin using the same estimator, seeds, and bootstrap as
H3, and report the number of records and families per bin. H3 predicts that
`R_c` rises from bin 1 to bin 2. This is exploratory: no gate, no alternative
bin edges after seeing the result, and no selection of a bin that makes the
trend look stronger.

## A10.5 — P3 train-and-test and the D1-core rule

P3 trains a separate head on each condition's **train-split** cache and
evaluates it on the same condition's withheld cache:

- **a arm:** the existing `runs/p2fix_seed{17,23,41}` heads. They were trained
  on `giga_a_train` with exactly the P3 hyperparameters, so retraining an
  identical arm would only add seed noise. They are reused, not reselected.
- **e-block arm:** train on `giga_block_e_train`, evaluate on
  `giga_block_e_test_withheld`.
- **b-block arm:** train on `giga_block_b_train`, evaluate on
  `giga_block_b_test_withheld`.

Training uses the frozen `splits/s1prime_filtered/train.jsonl` membership and
the pairs already in that file. The train split contains no Rfam seed rows, so
Amendment 8's v3 correction changes no training label. Evaluation swaps in v3
for seed rows at scoring time. No `id_test` or `test_withheld` record is used
for fitting, early stopping, or checkpoint selection; the validation split is
the same family-disjoint 5% of train used by `train_head_s1.py`, with the
condition's own cache.

The three arms are fitted on the **intersection of records covered by all three
train caches**, and evaluated on the intersection of records covered by all
three withheld caches. This prevents a cache that happens to drop long
sequences from giving one arm an easier population. The intersection sizes and
the ids dropped from the frozen membership are recorded with the P3 result; no
arm receives extra records.

The P3 hyperparameters are frozen to the corrected P2 run:

| setting | value |
|---|---:|
| epochs | 12 |
| learning rate | 2e-3 |
| weight decay | 0.01 |
| `pos_weight` | 20 |
| pair budget | 640,000 |
| validation fraction | 0.05 |
| selection metric | validation family-macro F1 |
| evaluation threshold | 0.5 |
| seeds | 17, 23, 41 |

The D1-core readout is the **withheld all-pairs family-macro F1**, using v3
labels, the ≥4-true-pair exclusion, and threshold 0.5. Nested-subset primary F1
and micro F1 are reported alongside but do not decide the rule.

For each block arm and seed, compute the family-level difference

```
Delta_f = F1_f(block-trained) - F1_f(a-trained)
```

on the same families, average `Delta_f` across the three seeds, and bootstrap
families 10,000 times with `numpy.random.default_rng(20261112)` for a 95%
percentile CI. Report the per-seed deltas and both arms' unpaired F1s.

**D1-core is falsified if neither the e-block nor the b-block arm has a
seed-averaged `Delta` of at least 0.02 with a 95% CI lower bound above zero.**
It is supported only if at least one arm meets both conditions. A point
estimate above 0.02 whose CI includes zero is not support under this rule.

Two limits are recorded now so they cannot be discovered later as excuses:

1. e-block changes the input distribution as well as removing distal content;
   b-block preserves local composition but not native distal content. Neither
   is a perfect counterfactual.
2. The FM was pretrained on RNAcentral, so family identity may live in the
   frozen backbone even when the head is trained on a block cache. P3 tests
   head-level reliance, not pretraining contamination; no split of ours can
   remove the latter.

## A10.6 — Calibration and H2.4 order

The old `(s, b)` fit is cancelled under Amendment 9 and is not reused.

1. Dump condition-a train logits with `runs/p2fix_seed17/best.pt`. Seed 17 is
   fixed here because it is the seed used for the existing D1/K1 and decoder
   dumps; it is not chosen by an H2.4 outcome.
2. Fit `(s, b)` under Amendment 7 §A7.4 on the frozen train split: v3 labels,
   maximum-nested-subset truth, ≥4 true pairs, 5-fold clan-grouped CV, the
   pre-registered 5×5 coarse grid, and the same 5×5 refinement if the coarse
   optimum lies on an edge.
3. Freeze the selected values with `scripts/freeze_energy_fit.py
   --write-lock`; the fit JSON is added to `state/amendments.lock.json` and
   committed before any H2.4 number is computed. A coarse-edge optimum is
   refined, never reported as final.
4. Before H2.4, run the D2C false-positive-helix energetics check on `id_test`
   and `test_withheld` with the same head and record the baseline counts,
   energies, and MFE-overlap fractions. These are expectations, not a gate:
   they tell us what H2.4 would have to change.
5. Run H2.4 on `id_test` and `test_withheld` under Amendments 3 and 4, and
   report its basis (`recall_matched`, `pair_budget`, or `null`). The
   `(s, b)` freeze is a precondition; a fit produced after seeing H2.4 is
   invalid.

This CPU work may run in parallel with the GPU block-cache builds, but its
output must not be inspected to choose a P3 seed, threshold, or geometry.

## A10.7 — Positive control

The near-duplicate positive control is **deferred until after P3**, per PI
Round 4k. It is not cancelled. No result from it may be used to reinterpret the
D1-core verdict.

## A10.8 — Reporting invariants and embargo

- Family-macro and pooled micro F1 appear side by side; a gap above 0.15
  without a stated reason fails the report.
- Every split includes a dinucleotide-shuffle floor; a floor above 0.05 fails
  the report as a scoring or cache defect.
- Round 4 numbers remain labelled **secondary-style, v1 labels including
  empty-truth rows** and never share a table with primary results.
- `test_novel` remains embargoed. No scoring output for it exists and none is
  produced until Amendment 10's outcomes are recorded and the PI authorises
  confirmatory scoring.
- No cache, checkpoint, threshold, seed, bin edge, or block geometry is changed
  after seeing a D1-core or H3 outcome. A change requires a new dated
  amendment and a fresh run.

## A10.9 — Freeze

This amendment is append-only after its first commit. Its sha256 is recorded in
`state/amendments.lock.json` before any P2-block cache is read. The
`first_quoted_in` field is `pending-commit` until the commit that carries it;
that field is provenance and does not change the locked bytes.
