# E1′ pre-registration — context-intervention factorial

Status: **frozen before any E1′ training run.** Amendments require a dated
entry at the bottom; results already produced under the older text stay reported
under it.

Source of the design: PI Round 1 (`docs/claude/transcripts/2026-11-11_round01_decision_and_E1prime.md`)
and Round 2 (`...round01b_e0_decisions.md`, Round 2 reply archived in this repo).

## Hypothesis

Distal sequence context in a frozen RNA foundation model acts as a family-identity
shortcut: it raises in-distribution accuracy and lowers out-of-distribution
accuracy. Removing it, and replacing the 2D decoder with exact nearest-neighbour
dynamic programming, recovers OOD accuracy without losing ID accuracy.

## Factors

| factor | levels |
|---|---|
| backbone | RiNALMo-giga (AGPL, isolated), RNA-FM (MIT), structRFM (MIT, contamination-flagged) |
| context half-width k | 8, 16, 32, 64, 128, 256, full (min(L, 1024)) |
| decoder | (A) learned energy corrections into nearest-neighbour DP; (B) matched-capacity 2D ResNet |
| distal-context condition (k fixed at 32) | a own, b dinucleotide-shuffled, c low-identity homolog (<60% id, same family), d other-clan, e pad/mask |

Decoders run for RiNALMo only; the other backbones use decoder A. The clan probe
(linear classifier, clan label from per-position embeddings) runs for every
backbone × k.

## Splits

- **Train:** bpRNA-1m(90) ∪ Rfam 15.0 seeds, minus S1′-novel clans, minus
  S1′-withheld clans, minus the ID-test families held out inside retained clans.
  Pair policy: canonical + GU.
- **ID test:** 10% of retained bpRNA families, sequence-deduped against train.
- **OOD-1 (k-curve fitted here):** S1′-withheld — clans present in bpRNA but
  deliberately removed from training. Same label provenance as train.
- **OOD-2 (confirmation and P5):** S1′-novel — clans with no bpRNA support at all;
  labels come from projected SS_cons, so it is confounded with Rfam's annotation
  conventions unless S4 agrees.
- **S4 (label-independent):** temporal PDB structures; P5 must hold here too.
- ArchiveII: reproduction (T4) and leave-one-clan-out (8 folds, tRNA+tmRNA merged)
  only. Not a primary OOD benchmark.

Split files are frozen by sha256 in `splits/s1prime/manifest.json` before training.

## Metrics

- Primary: **family-macro base-pair F1** on the held-out split, family as the
  bootstrap unit (10,000 resamples), 95% CI.
- Secondary: sequence-micro F1, precision, recall, Stem-F1, F1 by length bin,
  OOD/ID retention.
- Every ArchiveII number is reported both deduplicated and raw (313 duplicate
  sequences inside ArchiveII).

## Predictions and thresholds

| # | prediction | threshold |
|---|---|---|
| P1 | ID F1 roughly non-decreasing in k; OOD F1 lower at full context than at the best crop | OOD: full − crop(k*) ≤ −0.05, CI excludes 0, in ≥2 of 3 backbones |
| P2 | on ID data, homolog distal context keeps the full-context gain; other-clan context does not | (c) retains ≥70% of the (a)−(e) ID gain; (d) retains ≤30% |
| P3 | on OOD data the sequence's own distal context hurts | OOD F1: (b) ≥ (a) |
| P4 | how well the embedding encodes clan identity predicts the ID–OOD gap | Spearman ρ ≥ 0.7 across backbone × k |
| P5 | crop(k*)+DP beats BPfold on OOD, or beats it when added on top of BPfold's motif energies | Δ ≥ 0.03 family-macro F1 on S1′-novel **and** S4 |

## Go / no-go

- **Falsified → pivot to reactivity-only supervision (D2):** P1 fails in all
  backbones and P2 fails.
- **Mechanism-only (NAR route):** P1 and P2 hold, P5 fails; then one test of
  "BPfold + distal-ablated FM" before deciding.
- **Go:** P1, P2, P5 hold plus at least one of P3/P4.

## Compute plan and the measurement that gates it

The PI's estimate is ~120 GPU-h for embedding caches over 3 backbones × 7 k
values. That estimate is not yet measured, and crop embedding costs
`O(L · (2k+1))` residues per sequence, not `O(L)`, so it can be far larger.
`scripts/bench_backbone.py` measures seconds per residue for full context and for
k ∈ {8, 32, 128} on one 3090 **before** the grid is launched. If the extrapolated
cost exceeds 120 GPU-h, the grid is reduced in this order, with the reduction
recorded here before the run:

1. RiNALMo-giga → RiNALMo-micro for the sweep; giga only at k ∈ {32, full}.
2. Drop k = 16 and 64 (keep 8, 32, 128, 256, full).
3. Cap training sequences at length ≤ 512 and 8,000 sequences.

## What falsifies the mechanism independently of performance

- OOD F1 rising monotonically with k.
- Decoder A ≤ decoder B at matched k (the gain is the decoder, not the
  representation).
- Clan decodability flat in k while P1 still holds.

## Amendments

_None yet._
