# Claude Round 1 — Keep D1, recast as a causal paper

- Date: 2026-11-11
- Model: Claude Opus 5.5, High effort
- Session: https://claude.ai/chat/cc2ad40a-406e-481e-8e7b-89a64b78f2f4
- Executor input: independent verification of CHANRG, BPfold, structRFM,
  RiNALMo license/context limit; request for a decision among three options.

## Decision: option A, with one condition

Keep D1, recast as a causal paper about **why** foundation models fail out of
distribution, plus a fix. Reject B (reactivity-only supervision) and C
(calibration on a structRFM backbone) for now.

- C is incremental and builds on the most contaminated backbone available.
- B stays as the pivot if A's hypothesis is falsified; it costs about 3x the
  compute and carries identifiability risk.
- A is the cheapest path to a result nobody has: BPfold shows local energy +
  global inference works, but not why the FM approach fails or that removing
  context is what causes the OOD gain.

**Condition:** a diagnosis alone will not reach Nature Methods. The paper needs
two halves: (1) mechanism — distal context acts as a family-identity shortcut,
shown by intervention; (2) payoff — the fix (crop + exact DP) beats or
complements BPfold on clan-held-out data. If (1) holds and (2) does not, the
result is still a solid NAR / Genome Biology paper.

## Changes from Round 0

1. **Interventions replace correlation.** A k-sweep alone is confounded because
   smaller k also removes real long-range information. Add context-swap
   interventions with the local window fixed.
2. **Partner information stays.** Pair (i, j) is always scored from h_i and h_j,
   each embedded with its own local window. Crop/swap change the embeddings, not
   which partner the scorer sees.
3. **Length scope.** RiNALMo stops at 1024 positions, so "full context" means
   min(L, 1024). Primary analysis covers L <= 512. No claims beyond that.
4. **Backbones and licenses.** RiNALMo-giga (AGPL) behind a separate wrapper;
   RNA-FM (MIT); structRFM (MIT, frozen, contamination-flagged upper bound).
   Three backbones so the shortcut claim is not one model's quirk.

## E1' — context-intervention factorial

Training: bpRNA-1m(90) with every test clan removed (cmscan + identity dedup).
ID test: held-out families from training clans. OOD: S1' (clan-held-out from
Rfam 15.0) and S2 (ArchiveII leave-one-family-out, rebuilt).

**Arm 1 (k-sweep).** k in {8, 16, 32, 64, 128, 256, full} x 3 backbones x 2
decoders: (A) learned energy corrections into nearest-neighbour DP, (B) 2D ResNet
head of matched capacity. Decoders run for RiNALMo only; other backbones use
decoder A.

**Arm 2 (distal-context interventions).** Window fixed at +/-32; everything
outside replaced length-matched:

| condition | distal context | what it keeps |
|---|---|---|
| a. own | native | everything |
| b. dinucleotide-shuffled own | shuffled native | composition only |
| c. low-identity homolog (<60%, same family) | another family member | family identity, not instance sequence |
| d. other-clan sequence | from a different clan | neither |
| e. pad/mask | none | nothing |

**Arm 3 (clan probe).** Linear classifier predicts Rfam clan from per-position
embeddings for each backbone x k.

Baselines: BPfold (own and our splits), structRFM with our heads (flagged),
RNAstructure, EternaFold, MXfold2, Turner-only DP control.
Metrics: family-macro F1, Stem-F1, P/R, OOD/ID retention, 3 seeds,
family-level bootstrap CIs.
Budget: ~220 GPU-h (120 cache, 100 heads; DP on CPU).

## Pre-registered predictions

| # | prediction | threshold |
|---|---|---|
| P1 | ID F1 roughly non-decreasing in k; OOD F1 lower at full context than best crop | OOD full - crop(k*) <= -0.05, CI excludes 0, in >= 2 of 3 backbones |
| P2 | on ID, homolog distal context keeps the full-context gain; other-clan removes it | (c) retains >= 70% of (a)-(e) ID gain; (d) retains <= 30% |
| P3 | on OOD, the sequence's own distal context hurts | OOD F1: (b) >= (a) |
| P4 | embedding clan decodability predicts the ID-OOD gap | Spearman rho >= 0.7 across backbone x k |
| P5 (payoff) | crop(k*) + DP beats BPfold OOD on S1', or beats it when added on top of BPfold motif energies | delta >= 0.03 family-macro F1 |

## Go/no-go

- **Falsified, pivot to B:** P1 fails in all backbones and P2 fails.
- **Mechanism-only paper (NAR route):** P1 and P2 hold but P5 fails; then test
  "BPfold + distal-ablated FM" once before deciding.
- **Go (full paper):** P1, P2 and P5 hold, plus at least one of P3/P4. Next round
  is E2 (stem transplants) and S4 (temporal PDB).
- **Complication:** if P1 holds only with the ResNet decoder and not the DP, the
  effect is in the decoder rather than the representation.

## Splits

- **S1'** (new, primary): our own clan-held-out split from Rfam 15.0 with
  CHANRG-style methodology. Released with code. CHANRG numbers cited as context
  only, never as directly comparable.
- **S2:** ArchiveII leave-one-family-out rebuilt from raw ArchiveII, removing the
  test family's whole clan from the training pool. If BPfold's family-wise
  protocol did not do this, their 0.820 is inflated and that is worth reporting.
- **S3:** bpRNA-new from source, reported as defined and clan-filtered.
- **BPfold's splits:** only in T4 (baseline reproduction within +/-0.02), never
  for claims.

## Asks handed back to the executor

1. structRFM pretraining manifest and its overlap with every test set.
2. RiNALMo-giga positional handling: confirm crops shorter than 1024 need no
   special handling; check crop-vs-full embeddings agree at the sequence centre
   when k >= L/2.
3. Homolog pools for condition (c): count members below 60% identity per family;
   exclude families with fewer than 5 from P2.
4. Licenses: RNA-FM (verified MIT) and a repo-license decision.
5. BPfold training protocol: did their bpRNA training set exclude ArchiveII
   clans?

Evidence format unchanged; T1 gains an intervention column, plus two plots: OOD
F1 vs k (one line per backbone) and a bar chart of the five intervention
conditions split by ID and OOD.
