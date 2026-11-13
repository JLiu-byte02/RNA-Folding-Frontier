# Amendment 12 — pretraining-coverage hypothesis E-C1

Status: written 2026-11-13 **before computing either coverage proxy or the
coverage–F1 relationship**. It carries PI Round 4l. No one has inspected the
relationship between Rfam full-region membership, RiNALMo pseudo-perplexity,
and family-level FM advantage before this file was pinned.

This amendment is separate from the D1/H2 pre-registrations. It does not
reopen D1, H2, or D2. `test_novel` remains embargoed until the analysis list
in §A12.6 is fixed and this amendment is committed.

## A12.1 — Hypothesis

An FM-based predictor's advantage over a thermodynamic baseline on a family
scales with how well that family is represented in the FM's pretraining
corpus. The useful family signal lives in pretraining, not in head training.

The test is on the 409 families of the P3/D1 common withheld population, where
the head has no family shortcut but the FM may have seen the clan during
pretraining.

## A12.2 — Unit and outcome

**Unit:** family. The family set is the 409 families in the frozen P3 common
withheld population: records present for a, e and b in all three seeds, with
at least four v3 true pairs. The family list is written to
`experiments/002_e1prime_runs/ec1_families.json` after this amendment is
pinned and is not selected by any outcome.

**Per-sequence scores:**

- FM arm: `dump_p2fix_seed17/a_test_withheld.jsonl`, greedy nested decoder at
  fixed threshold 0.5, truth = the maximum nested subset of the v3 canonical
  set. Per-sequence F1 uses `pair_prf` with the Amendment 9 convention.
- Thermodynamic arm: ViennaRNA `RNA.fold(seq)` MFE structure, parsed to pairs;
  per-sequence F1 against the same nested v3 truth.

**Outcome:** for each family,

```
DeltaF1_f = mean_seq F1_FM(f) - mean_seq F1_RNAfold(f)
```

The family-macro FM advantage is the mean of `DeltaF1_f` across families.
This controls for how intrinsically hard each family is.

## A12.3 — Coverage proxies, fixed in advance

### P1 (primary): Rfam full-region membership

`P1_f = log(1 + full_region_membership_f)`, read from column 16 of
`data/raw/rfam/family.txt`, keyed by Rfam accession. Column 16 is the
full-region sequence count (for RF00001 it is 594,154; column 15 is the seed
count 712). No FM output enters P1.

### P2 (secondary): RiNALMo-giga masked pseudo-perplexity

The released checkpoint has `RiNALMoForPreTraining` and `lm_head.*` weights
(504 tensor keys total), so P2 is available.

For each sequence:

1. tokenize with the released RiNALMo tokenizer;
2. mask 15% of non-special positions, using four fixed masks with seeds
   17, 23, 41 and 53;
3. compute the mean cross-entropy (natural log) on masked positions only;
4. `pppl = exp(mean NLL)`;
5. average `pppl` over the four masks.

The family value `P2_f` is the median over its sequences. The model is frozen;
no fine-tuning and no head is trained.

If the MLM forward cannot be reproduced from the released weights for any
reason, P2 is dropped and only P1 is reported; no substitute proxy is invented
after seeing outcomes.

## A12.4 — Covariates and partial correlation

Family-level covariates, fixed before looking at the outcome:

- median sequence length;
- mean GC content;
- mean v3 pair density (`pairs / length`);
- multiloop presence: whether the maximum nested v3 structure contains a loop
  with more than one child helix;
- label source (`bpRNA` vs Rfam seed), encoded as a binary indicator.

Compute Spearman partial correlation between `DeltaF1_f` and `P1_f`, adjusting
for all covariates by rank-residualization: rank-transform each variable,
regress the `DeltaF1` ranks and `P1` ranks on the covariate ranks with an
intercept, and take the Pearson correlation of the two residual vectors. The
same procedure is used for P2 against `DeltaF1` (without the sign flip; the
reported P2 statistic is the partial correlation with pppl).

The family bootstrap resamples the 409 families 10,000 times with
`numpy.random.default_rng(20261113)`, recomputing the full partial correlation
in each resample.

## A12.5 — Success and falsification

- **E-C1 supported** if the partial `rho(DeltaF1, P1)` is at least 0.30 and
  its 95% family-bootstrap CI excludes 0, and P2 has the predicted negative
  direction (`rho <= -0.20`).
- **E-C1 falsified** if `rho < 0.10` or the CI includes 0. C is dropped and
  Paper A is written without it.
- A result between 0.10 and 0.30, or a positive P2 direction, is
  **inconclusive**: it is reported with the pre-registered thresholds and is
  not used to claim support.

## A12.6 — One-shot test_novel confirmation

If E-C1 passes, the same P1/P2 outcome and partial-correlation analysis is run
once on the frozen `test_novel` split. Before that run, the complete analysis
list is fixed in this amendment and the test_novel labels are unsealed once:

1. baselines: RNAfold, EternaFold, RNAstructure, MXfold2 retrained on the
   frozen train split, and BPfold with released weights (overlap flagged);
2. our greedy head and A-energy-fixed decoder;
3. E-C1 confirmation, and E-C2 if it passed;
4. audit summary statistics.

Every output is hashed immediately; no analysis is re-run after unsealing
except to fix a documented bug, recorded in `docs/corrections.md`.

## A12.7 — E-C2, only if E-C1 passes

Make the A-energy-fixed scale sequence-dependent:

```
s(x) = s_max * sigmoid(alpha - beta * pppl(x))
```

Fit `(s_max, alpha, beta)` by 5-fold clan-grouped CV on the frozen train split,
using the same objective as Amendment 7. Success: the gated decoder beats both
fixed-s decoder A and the best external baseline by at least 0.02 family-macro
F1 on withheld, with a family-bootstrap CI excluding 0, and confirms on
`test_novel`. No E-C2 run is allowed if E-C1 fails.

## A12.8 — What this amendment does not change

- D1-core remains falsified; H3 remains below its bar; H2.4 remains a
  controlled negative.
- The `(s,b)` calibration is frozen at `s=0.1875, b=0`.
- No D2/reactivity work is started; B is shelved.
- `test_novel` remains embargoed until §A12.6 is fixed and this amendment is
  committed.
