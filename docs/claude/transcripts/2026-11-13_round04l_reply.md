# Round 4l reply — Claude Opus 5.5 (High), PI

Source: https://claude.ai/chat/cc2ad40a-406e-481e-8e7b-89a64b78f2f4
Captured: 2026-11-13, after the D1-core, H3, and H2.4 report.

## Round 4l: decision. Write paper A now, and run one bounded test of a new hypothesis (C). B is shelved.

All three pre-registered outcomes came back negative, and the work behind them
was done carefully: common populations, fixed splits, frozen `(s, b)`, and the
fallback rule applied as written. That makes these negatives publishable. It
also means D1 can't carry a Nature Methods paper. Here's how I'd proceed.

### Decision

**A is the committed deliverable, starting now.** It's a real paper on its own:
a contamination and label audit plus well-controlled negative mechanism
results. The realistic venues are NAR, Bioinformatics or Genome Biology, not
Nature Methods.

**C gets exactly one cheap, pre-registered test** (about 2 GPU-h plus CPU time)
that could turn A into a mechanism-and-method paper. If it fails, we stop and
write A without it. No further hypotheses after that.

**B (reactivity-only supervision) is shelved, not just deferred.** Three
reasons:

- It builds on none of our findings.
- Its gap against EternaFold, RibonanzaNet and Krueger & Aviran is narrow.
- It's the most expensive option and has an unresolved identifiability risk.

Reopen it only if you specifically want a method attempt after A is written.

### Hypothesis C: pretraining coverage

Our results point here directly:

- The head doesn't take a family shortcut (P3 is null with a ±0.009 CI).
- Yet homologous distal context retains most of the performance
  (`R_c ≈ 0.75–0.9`) and unrelated context doesn't.
- Every withheld clan was unseen by the head but almost certainly seen by the
  FM during pretraining on RNAcentral.
- CHANRG's hardest split (GenF) is defined as families with sparse
  phylogenetic coverage.

**Hypothesis.** An FM-based predictor's advantage over thermodynamics on a
family scales with how well that family is represented in the FM's pretraining
corpus. The useful family signal lives in pretraining, not in head training.

**Novelty check needed [VERIFY]:**

- I believe ESMFold reported, for proteins, that structure accuracy tracks
  language-model perplexity, but I'm not certain of the exact claim. Check Lin
  et al., Science 2023, and related protein work on likelihood and performance.
- I know of no RNA analogue. If one exists, C becomes a replication: still
  useful for A, but weaker.

### E-C1: the decisive test, on withheld

**Unit:** the 409 common-population withheld families.

**Outcome:** family-macro FM advantage,
`ΔF1 = (p2fix head, greedy, condition a) − (RNAfold)`, on v3 nested labels.
Measuring the advantage over RNAfold rather than raw F1 controls for how
intrinsically hard each family is.

**Coverage proxies, both fixed in advance:**

- **(P1, primary):** log of Rfam 15.0's full-region membership count for the
  family, from `family.txt`. It costs nothing and is independent of the FM.
- **(P2, secondary):** RiNALMo-giga masked pseudo-perplexity: mask 15% of
  positions, average over 4 masks with fixed seeds, take the family median.
  Roughly 1–2 GPU-h. First check that the MLM head weights are released; if
  they aren't, P2 is dropped.

**Covariates for partial correlation:** median length, GC content, pair
density, multiloop presence, and label source (bpRNA vs seed).

**Success:** Spearman partial `ρ(ΔF1, P1) ≥ 0.30`, with a family-bootstrap 95%
CI excluding 0. P2 should show the same direction (`ρ ≤ −0.20` against
perplexity).

**Falsified:** `ρ < 0.10`, or a CI that includes 0. C is dropped, and A is
written without it.

**Pre-register first, in Amendment 12:** the proxies, outcome, covariates,
thresholds, and a confirmation run on `test_novel` with the identical
analysis. Pin it before computing either proxy, or ΔF1 on withheld, and state
in the amendment that nobody has looked at the coverage–F1 relationship before
pinning.

### E-C2: only if E-C1 passes, a coverage-gated decoder

The scale in A-energy-fixed becomes sequence-dependent:
`s(x) = s_max · σ(α − β · pppl(x))`. The decoder trusts the FM where coverage
is high and falls back toward Turner where it's low. Fit `(s_max, α, β)` by
the same clan-grouped 5-fold CV on train as Amendment 7.

**Success:** it beats both fixed-s decoder A and the best baseline (BPfold,
EternaFold, or MXfold2 retrained on our split) by ≥ 0.02 family-macro F1 on
withheld, with a CI excluding 0, and the result is then confirmed on
`test_novel`.

This is the only route left that offers both a mechanism and a method
together.

### What remains publishable from D1 and H2

| Claim | Status | Strength |
|---|---|---|
| ArchiveII is 61% exact / 88% near-duplicate of RNAStrAlign ∪ TR0. RNAStrAlign labels match ArchiveII's almost exactly (F1 0.994), so trained models can reproduce them; bpRNA labels cap agreement at 0.92 (telomerase 0.66) | Audit, solid | Core result of A |
| ArchiveII spans 8 clans; leave-one-family-out leaks through tRNA/tmRNA; Rfam-delta families leak through clan siblings (Cobalamin: 74 clan-sibling hits) | Audit, solid | Core |
| Projecting Rfam seed structures as strings corrupts labels (17.8% unparseable; 86k differing pairs vs pair-level projection); bpRNA labels contain crossing pairs in 23.7% of sequences | Audit, solid | Core. [VERIFY] which published Rfam-derived benchmarks used string-level projection; if any did, this becomes a headline finding |
| Removing or scrambling distal context in head training doesn't change OOD F1 (Δ = +0.0007 [−0.008, +0.009]) | Pre-registered negative, well powered | Core mechanism result: rules out head-level shortcut learning as the cause of FM OOD failure |
| At test time, homologous distal context mostly substitutes for native context; unrelated context doesn't (`R_c` 0.75–0.90, `R_d` 0.08–0.23) | Below the H3 bar | Report as directional and descriptive, with the pre-registered bar stated |
| Long-range false positives are stable helices: 4× enriched in stacks, isolated energy −1.43 vs −1.23 kcal/mol for true helices, only 26% in the MFE structure | Descriptive, solid | Explains why energy decoding fails to prune them |
| Energy-consistent decoding does not reduce long-range false positives (H2.4 negative), but at a matched pair budget it adds 3% (id_test) and 7.4% (withheld) true positives | H2.4 negative is confirmatory-grade; the TP gain is exploratory | A controlled answer to Flamm et al.'s open question, with a nuance |

**Not publishable:** anything from the original Round 4 table, H2 as a
mechanism, the 102% figure, or any number scored with `pair_prf`. The scorer
bug goes in the methods as a corrected internal error, not as a finding.

Two numbers still missing for the table:

- P3's id_test results. If e-block ≈ a there too, the finding broadens to "for
  this head, full FM context adds nothing over ±32 nt", which ties into
  CHANRG's wiring failure.
- H2.4's per-length-bin results, to see where the true-positive gain comes
  from.

### Unblinding test_novel

`test_novel` will be unblinded once, after Amendment 12 fixes the complete
analysis list:

- Baselines: RNAfold, EternaFold, RNAstructure, MXfold2 retrained on our train
  split, and BPfold with released weights, flagged for overlap.
- Our greedy head and decoder A-energy-fixed.
- The C confirmation, plus E-C2 if it passed.
- The audit-relevant summary statistics.

Every analysis gets run in one pass, the outputs are hashed immediately, and
nothing is re-run afterwards except to fix documented bugs, each recorded in
`corrections.md`.

### Next report

1. The Amendment 12 hash.
2. The E-C1 verdict for P1, P2 if available, with partial ρ and CI.
3. P3's id_test results and H2.4 by length bin.
4. MXfold2 retraining status on our train split.
5. A one-page outline of paper A: figure list, mapping each figure to a table
   row above.

If E-C1 fails, the round after that is writing.
