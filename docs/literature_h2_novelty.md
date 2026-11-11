# Novelty check for H2 ("spurious partners"), 2026-11-12

Purpose: decide what of H2 can be claimed as new before any result exists, so
the claim is fixed by the literature rather than by what the data happens to
show. Every statement below is either marked **verified** (read in the source
during this check) or **unverified lead**.

## What is already established

### 1. DL secondary-structure models fail to generalise across families — verified

Szikszai, Wise, Datta, Ward, Mathews, *Bioinformatics* (2022),
doi:10.1093/bioinformatics/btac415, "Deep learning models for RNA secondary
structure prediction (probably) do not generalize across families".

- They show it is "nearly trivial with convolutional neural networks to generate
  pseudo-free-energy changes ... that improve the accuracy of structure
  prediction for intra-family cases", and that intra-family performance is
  "insufficient proof of generalization".
- They propose a stricter inter-family cross-validation protocol.

**Consequence for us:** the family shortcut is *not* our discovery. H1's value is
that it is a *causal* test (distal-context swap) on a frozen foundation model,
not a correlation between train/test family overlap. Any claim text must cite
this paper and start from it.

### 2. The field has already reorganised around that failure — verified

*Machine Learning for RNA Secondary Structure Prediction: a review of current
methods and challenges*, arXiv:2511.02622 (2025).

- Describes the "generalization crisis": models that did well on TS0 failed on
  new families, which "mandated the development of more rigorous, homology-aware
  benchmarks" (bpRNA-new, PDB-derived TS2).

**Consequence for us:** S1′ (clan-held-out) is a *strengthening* of an accepted
protocol, not a new protocol. It cannot carry the paper; the mechanism and the
payoff must.

### 3. Models can beat physics-style folding and still not learn the physics — verified

Flamm, Wielach, Wolfinger, Badelt, Lorenz, Hofacker, *Frontiers in
Bioinformatics* 2:835422 (2022), "Caveats to Deep Learning Approaches to RNA
Secondary Structure Prediction".

- Argues models must learn a (complex) biophysical model rather than a
  phenomenological fit, highlights training-data biases, and proposes testing on
  synthetic data from a simplified biophysical model where biases are
  controllable.

**Consequence for us:** a "the model is not using the physics" claim is
publishable only if the diagnostic separates *what the model computed* from
*how it was thresholded* — which is exactly the job of K1 and the three-decoder
comparison. This is a point in favour of the current D1 design.

### 4. Shorter context can be better than global folding — verified for a different task

Lange et al., *Nucleic Acids Research* (2012), "Global or local? Predicting
secondary structure, accessibility and binding sites in mRNAs".

- Window-based local folding is reported to outperform global folding for
  accessibility and binding-site prediction in long mRNAs.

**Consequence for us:** "local beats global" is old wisdom in RNA
thermodynamics. Our version differs in three ways that must be stated
explicitly, or a reviewer will call it a re-discovery: (i) the object is a
learned 2D pair map from a frozen foundation model, not a thermodynamic
ensemble; (ii) the metric is base-pair F1 under family-held-out splits, not
accessibility; (iii) the manipulation is causal and per position (P2-block),
not a choice of algorithm.

## What H2 adds that is not covered above

| claim | status after this check | what has to be true to keep it |
|---|---|---|
| H2.1 long-range FP rate grows with length | **not verified as published**; the closest published statement is about generalisation, not length scaling | must be measured by us (D2), and reported as a finding even if it is only confirmatory |
| H2.2 those FPs are helix-shaped beyond a distance-matched null | **no published equivalent found** | the null must be the pre-registered one; a bare fraction is uninterpretable |
| P2-block: distal context is net-negative below ~300 nt | **no published equivalent found** | needs K1 to survive, and the length bands to be reported, not just the mean |
| K1: the e−a gap is an operating-point artefact | **no published equivalent found** | this is the methodological contribution if it fires; if it does not fire, it is the control that makes the e−a claim credible |

## Unverified leads (do not cite until checked)

- A frequently repeated claim in talks and reviews that DL predictors output a
  number of base pairs growing superlinearly with length. Plausible and worth
  checking in Flamm et al. §3, but **not** confirmed here, so it must not appear
  in any draft yet.
- Whether UFold-X (long-range U-Mamba, 2025) and similar long-range models
  actually improve *family-held-out* F1 rather than TS0/ArchiveII F1. If they
  do, the "distal context is net-negative" line needs to be scoped to
  frozen-backbone heads rather than stated generally.

## How this changes the paper's story

The story cannot be "we discovered that models memorise families" or "we
discovered that local context can beat global folding". Both are taken.

The defensible story, if K1 and H2 survive, is:

> On a family-held-out benchmark, a frozen RNA foundation model plus a 2D head
> loses accuracy because distal sequence supplies complementary partners that
> the head cannot reject; restricting each position to its own local context
> removes those partners and recovers accuracy below ~300 nt. The gain is not a
> threshold artefact: it survives an oracle operating point, and it disappears
> on a sequence-level control where the family shortcut is available.

That is a mechanism claim with a causal test, a calibration control, and a
positive control — none of which the three cited papers provide.
