# Amendment 14 - One-shot `test_novel` analysis list

Status: written 2026-11-13, **before any `test_novel` prediction, score, or
per-sequence output exists**. At the time of writing the only `test_novel`
quantities that have been read are label counts and cohort composition, which
Amendment 8 and Amendment 9 permit. `test_novel` is spent by the first run
authorised here; there is no second unblinding.

## A14.1 - What this cohort can and cannot decide

| quantity | value |
|---|---:|
| sequences in the evaluation subset | 2,372 |
| families | 310 |
| clans | **6** |
| Amendment 1.6 clan gate | 150 |
| `gate_passed` | **false** |
| v3 true pairs, total | 57,303 |
| records excluded by the v3 `<4`-pair rule | 0 |
| records with an empty v1 label repaired by v3 | 463 (19.5%) |
| sequences with `L < 300` | 2,315 |
| sequences with `L >= 300` | 57 |

Consequences, fixed now:

1. **This is an external-validity check, not a second confirmatory cohort.**
   With 6 clans the cohort fails the pre-registered diversity gate, so a
   `test_novel` pass cannot convert a failed withheld endpoint into a confirmed
   claim, and a failure cannot by itself refute an in-distribution effect.
   What it can do is stop a withheld pass from being dismissed as
   in-distribution overfitting.
2. **Intervals.** The primary interval is a family bootstrap over 310
   families. A clan bootstrap over 6 clusters is reported for completeness but
   is labelled uninformative and may not be quoted as a diversity-controlled
   interval.
3. **Selection caveat.** The evaluation subset required `min_pairs = 4` on the
   consensus projection, so it contains no v3-sparse rows. The full
   `test_novel` seed split has 42,873 rows, of which 3,675 (8.57%) have zero v3
   pairs; those rows are not scored and this cohort therefore says nothing
   about weakly structured or unstructured RNA.
4. **Length caveat.** Only 57 sequences have `L >= 300`. Amendment 13 A13.3
   states the `>=300` nt failure from withheld/id_test as a known limitation;
   `test_novel` cannot test it.

## A14.2 - Frozen prediction list

Exactly these items are produced, once each, in this order. No item may be
added, dropped, retuned, or reseeded after the first output is read.

| id | predictor | frozen settings |
|---|---|---|
| P1 | greedy head | `p2fix_seed{17,23,41}`, threshold .5, nested projection, all pairs; per-seed values and their mean; no seed selection |
| P2 | decoder A, fixed scale | `s = .1875`, `b = 0` from `experiments/002_e1prime_runs/energy_fit_refined_p2fix.json` (sha256 `be54ddab...ad5ad5`) |
| P3 | decoder A, energy-fixed | ViennaRNA 2.8.0 energy term with the same frozen `(s, b)` |
| P4 | RNAfold | ViennaRNA 2.8.0 defaults, as in `docs/rnafold_baseline_2026-11-13.md` |
| P5 | EternaFold | the version and manifest recorded by the baseline job |
| P6 | RNAstructure | the version and manifest recorded by the baseline job |
| P7 | MXfold2 retrained | checkpoint from our frozen train split; manifest carries the train-split sha256, MXfold2 commit, and README-default hyperparameters |
| P8 | BPfold | reported and **flagged**: TR0 contains our withheld clans, so it can never define "best"; if the artifact is unavailable it is reported as unavailable, not substituted |
| P9 | E-C2 G1 | frozen `(s_max, alpha, beta)` from `ec2_fit.json` (sha256 `7008310f...`); primary gate variant |
| P10 | E-C2 G2 | frozen secondary variant; reported, cannot carry a confirmatory claim |
| P11 | one-hot control heads | `onehot_seed{17,23,41}`, same recipe and splits, for E-C1b |
| P12 | `P1` covariate | per-family Rfam full-region sequence count, analysis only, never a decoder input |

The E-C2 scale-control arm at `s = 0.09375` is **exploratory**, is run after
P1-P12, and may not enter any verdict.

## A14.3 - Frozen endpoints

Primary endpoint for every predictor: **family-macro F1 on v3 nested all-pairs
truth**, mean over the three seeds for learned predictors. Secondary: micro F1,
per-length-bin descriptive table, and the pre-registered subgroup tests below.
Every table reports macro and micro together and the dinucleotide-shuffle floor
of the same architecture on this cohort.

Pre-registered subgroup endpoints (unchanged from their amendments):

- **H2.4 / A13.3**: for `L < 300`, decoder A at the frozen pair-budget point
  must increase true positives and not increase long-range (>128 nt) false
  positives relative to greedy nested decoding, in each length bin
  `<80`, `80-149`, `150-299`. The `80-149`-only statement remains the fallback
  wording if the full subgroup fails.
- **E-C1**: partial `rho(DeltaF1, P1)` on the 310 families, `DeltaF1 = F1_FM -
  F1_RNAfold`, with the Amendment 12 covariates and family bootstrap.
- **E-C1b**: `Delta_rho = rho_FM - rho_onehot >= .15` with a paired
  family-bootstrap CI excluding 0. MXfold2 (P7) is the secondary control.
- **E-C2**: G1 must beat both P2 and the best clean baseline among P4-P7 by at
  least .02 family-macro F1 with the family-bootstrap CI excluding 0.

## A14.4 - One-shot execution rules

1. P1-P12 run once. Outputs go to `experiments/012_test_novel_oneshot/`, each
   JSON hashed, with a manifest recording the code commit, the label hash
   `b8fe71ab...`, the evaluation-subset hash `9d1ab36c...`, the fit hashes, and
   the scorer version.
2. After the first output is read, no threshold, scale, seed, subset, or
   decoding rule may change. If a code defect is found, stop, document the
   defect and fix in `docs/corrections.md`, rerun **every** affected item, and
   report explicitly that the unblinding was not clean.
3. `test_novel` is never used for model selection. After this run it is spent:
   later work may quote the frozen outputs but may not re-tune against them. A
   new claim that needs a fresh cohort needs a new held-out source, not a rerun
   here.

## A14.5 - Decision rules fixed in advance

| id | condition | licensed statement |
|---|---|---|
| D14.1 | G1 beats `max(P2, P4-P7)` by `>= .02` family-macro with CI excluding 0 on `test_novel`, **and** the withheld endpoint passed | coverage-gated decoding is supported out of distribution, with the 6-clan caveat stated in the abstract |
| D14.2 | G1 passes withheld but not `test_novel` | E-C2 is in-distribution only; the paper's headline is (ii) or (iii) |
| D14.3 | G1 fails withheld (the current frozen fit collapses to a constant scale) | E-C2 is reported as a negative/edge finding; no `test_novel` gating claim is made whatever its number |
| D14.4 | any learned predictor is below RNAfold on `test_novel` family-macro F1 | no "beats thermodynamics" wording anywhere in the paper |
| D14.5 | H2.4 fails any `L < 300` bin | report split-dependent decoding; only the `80-149` statement may be generalised |
| D14.6 | `Delta_rho < .15` or its CI includes 0 | wording stays "learned-model advantage tracks family abundance"; no pretraining attribution |
| D14.7 | any predictor's macro/micro gap exceeds .15, or its shuffle floor exceeds .05 | the corresponding number is reported as unstable and cannot support a claim |

## A14.6 - Not authorised by this amendment

- No architecture change, extra seed, ensemble, calibration, or threshold sweep
  after seeing `test_novel`.
- No per-family error analysis used to design the next method on this cohort;
  `test_novel` is spent at the end of the one-shot.
- No substitution of a missing baseline by a weaker one; a missing baseline is
  reported as missing.
- No claim that the cohort is clan-disjoint in the Amendment 1.6 sense.
