# Amendment 13 — E-C1b specificity control, coverage-gated decoder, and Paper A

Status: written 2026-11-13 **before any E-C1b or E-C2 outcome was computed**.
It carries PI Round 4m. E-C1 passed under Amendment 12; this amendment decides
what that pass is allowed to mean and freezes the E-C2 fit. `test_novel`
remains embargoed until Amendment 14 fixes the final one-shot analysis list.

## A13.1 — E-C1b: is the effect specific to pretraining?

The E-C1 outcome is `DeltaF1 = F1_FM - F1_RNAfold` per family. Abundant
families could advantage a learned model for reasons unrelated to pretraining
(easier/cleaner comparative labels), or could make RNAfold look worse. The
specificity test splits the correlation and adds a matched no-pretraining
control.

### Decomposition (descriptive)

On the same 409 families and with the same covariates as Amendment 12, compute
separately:

- partial `rho(F1_FM, P1)`;
- partial `rho(F1_RNAfold, P1)`.

This is descriptive; it is reported whichever way the primary control goes.

### Matched no-pretraining control (primary)

Train the same `PairMapHead` architecture and recipe with **learned one-hot
ribonucleotide embeddings** instead of FM embeddings:

- input dimension 4 (A,C,G,U; non-ACGU characters map to an all-zero row);
- same frozen train split and family-disjoint 5% validation split;
- same 12 epochs, learning rate 2e-3, weight decay .01, `pos_weight=20`,
  pair budget 640,000, threshold .5;
- seeds 17, 23, 41.

No pretrained weights enter this control. Evaluate it on the same withheld
records and v3 nested truth, form `DeltaF1_onehot = F1_onehot - F1_RNAfold` per
family, and compute the same partial correlation with P1.

**Criterion:** `Delta_rho = rho_FM - rho_onehot >= 0.15`, with a paired
family-bootstrap 95% CI excluding 0. If it passes, "pretraining coverage" is
earned. If it fails, the finding is reported as "learned-model advantage
tracks family abundance", and no pretraining attribution is made. The
secondary control is MXfold2 retrained on our train split, run through the
same analysis.

## A13.2 — E-C2: coverage-gated decoder

Only sequence-computable gate inputs are allowed. No Rfam family lookup and no
family median enters a deployable decoder.

### Gate

For a sequence `x` with per-sequence masked pseudo-perplexity `pppl(x)`:

```
z(x) = (pppl(x) - mean_train_pppl) / sd_train_pppl
G1: s(x) = s_max * sigmoid(alpha - beta * z(x))
G2: s(x) = s_max * sigmoid(alpha - beta * z(x) - gamma * log L)
```

`mean_train_pppl` and `sd_train_pppl` are computed once on the frozen train
split and frozen. G1 is primary; G2 is secondary because its length term was
motivated by withheld bins.

### Fit

- Train-split pseudo-perplexity uses the same four masks, 15% masking and
  seeds 17/23/41/53 as E-C1.
- Grid: `s_max in {0.1875, 0.5, 1.0}`, `alpha in {-1, 0, 1}`,
  `beta in {0, 1, 2}`; G2 adds `gamma in {0.25, 0.5}` with the same
  `(s_max, alpha, beta)` grid.
- Selection: 5-fold clan-grouped CV on the frozen train split, same objective
  and folds as Amendment 7. If the full 80-point grid exceeds 24 CPU-hours,
  use a family-stratified 3,000-sequence subsample with seed 17 and record its
  hash.
- Freeze the selected parameters and the train pppl mean/sd in
  `experiments/002_e1prime_runs/ec2_fit.json`, hash it, and record it in
  `state/amendments.lock.json` before scoring withheld.

### Success

On withheld, G1 must beat **both** fixed-s decoder A and the best clean
baseline (RNAfold, EternaFold, RNAstructure, or MXfold2 retrained on our
train split) by at least 0.02 family-macro F1, with a 95% family-bootstrap CI
excluding 0. BPfold may be reported but is flagged and cannot define "best"
because TR0 contains our withheld clans. The result then requires the
test_novel confirmation in Amendment 14.

Also report fixed-s A-energy-fixed against RNAfold and EternaFold on withheld.

## A13.3 — H2.4 length subgroup

The withheld bins are exploratory because they motivated this rule. The clean
test is `test_novel`:

**For sequences with L < 300, at the pair-budget operating point, decoder A
reduces long-range false positives and increases true positives relative to
greedy nested decoding.** The ≥300-nt failure is a known limitation and is
reported separately.

## A13.4 — Paper A structure

Write now: methods (splits, v3 projection, scorer convention and the corrected
internal error, pre-registration history), audit F1–F3, and the negative
mechanism results. P3 is stated as "for this head, full FM context gives no
gain over ±32 nt, in distribution (Δ −0.0095 [−0.025,+0.006]) or OOD"; H2.4
length bins are descriptive.

Hold the title, abstract and discussion until E-C1b, E-C2 and test_novel. The
headline is one of: (i) pretraining coverage plus a gated decoder; (ii)
pretraining coverage as a mechanism only; (iii) audit plus controlled
negatives. Keep a placeholder section for C with wordings for (ii) and (iii).

## A13.5 — Order and embargo

1. CPU: MXfold2 retraining on the frozen train split, with train hash, MXfold2
   commit and README-default hyperparameters in the manifest.
2. GPU: train-split pseudo-perplexity and one-hot heads.
3. CPU: E-C2 fit and freeze.
4. E-C1b decomposition and Δρ on withheld.
5. E-C2 evaluation on withheld.
6. Amendment 14 final test_novel analysis list, then one unblinding.

No E-C1b/E-C2 number may be quoted before this amendment is committed and its
hash is in the lockfile. `test_novel` remains unscored until Amendment 14.
