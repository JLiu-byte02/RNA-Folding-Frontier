# Paper A — one-page outline (2026-11-13)

**Working title:** *What RNA secondary-structure benchmarks actually measure:
leakage, label projection, and a controlled negative for foundation-model
context.*

**Thesis.** The apparent advantage of FM-based RNA folding on common benchmarks
is not explained by a head-level family shortcut. The benchmark and label layer
is itself a major source of apparent performance, and a carefully controlled
energy decoder does not remove the long-range false positives that motivate it.
The paper contributes an audit and well-powered negative mechanism results, not
a SOTA claim.

Target venues: NAR, Bioinformatics, or Genome Biology. Not Nature Methods.

**Current state (2026-11-15).** E-C1 passed (P1 rho .338, P2 rho -.242) but
E-C1b failed specificity (d_rho .085 < .15), so coverage is a shared predictor,
not evidence for a pretraining-specific mechanism. E-C2 failed at its boundary.
The four-method common cohort and retrained MXfold2 are frozen. Sections 1-4,
Tables T1-T7, and Figures 1/2/6/7 are drafted from committed artifacts; title,
abstract, discussion and conclusion remain held for the PI decision.

**PI framing (Round 4n, 2026-11-15).** Audit-led, with the clean clan-held-out
benchmark result at the centre. The paper is not a method claim, not a claim
that FMs fail, and not a coverage-specific mechanism. Position it against
CHANRG's out-of-distribution collapse. NAR/Genome Biology are the realistic
venues; Bioinformatics if compressed. Title and abstract stay held until H4
and the single `test_novel` unblinding are complete.

## Figures and claims

| figure | content | claim / table row |
|---|---|---|
| F1 | Benchmark contamination map: ArchiveII vs RNAStrAlign ∪ bpRNA TR0 exact/near-duplicate overlap; clan and family leakage (tRNA/tmRNA; Cobalamin clan siblings) | Audit: 61% exact / 88% near-duplicate; ArchiveII spans 8 clans; leave-one-family-out leaks |
| F2 | Label projection: string-level vs pair-level Rfam seed projection; unparseable rows; pair differences; bpRNA crossing-pair frequency | Audit: 17.8% unparseable; ~86k differing pairs; 23.7% crossing bpRNA sequences; verify which published benchmarks used string projection |
| F3 | Head-level intervention: P3 a/e-block/b-block design and withheld ΔF1 with family-bootstrap CIs | Mechanism negative: e Δ=+.0007 [−.0079,+.0093]; b Δ=−.0099 [−.0159,−.0036]; D1-core falsified |
| F4 | Test-time context substitution (P2-block a–e) and descriptive H3 ratios; common-population P2-central ratios | Secondary/descriptive: c retention ~.75–.90; d ~.08–.23; below the pre-registered H3 bar; the old 102% was a population mismatch |
| F5 | Long-range FP helix energetics and H2.4: greedy vs A-energy DP at matched pair budget, by split and length bin | Decoder negative: 4,391 FP helices; 25.6% in RNAfold MFE; DP adds TPs but also long-range FPs on both splits |
| F6 | E-C1: Rfam full-region membership and RiNALMo pseudo-perplexity vs family FM advantage; matched one-hot control | Coverage predicts both heads; the pretraining-specific control failed |

## Tables

- T1: dataset provenance, split sizes, exact/near-duplicate and clan-leakage
  counts.
- T2: label projection agreement and eligibility counts (v1/v2/v3; empty-truth
  and `<4` exclusions).
- T3: P3 common-population a/e/b family-macro and micro F1 with paired deltas
  and CIs.
- T4: P2-block a–e test-time table and H3 ratios with CIs.
- T5: H2.4 pooled counts, basis (`pair_budget`), and per-length-bin results.
- T6: publishable-claim table from the PI decision, with status/strength.

## Claims that are allowed

- ArchiveII is heavily duplicated against common trainers and its clan coverage
  is narrow; leave-one-family-out is not clan-disjoint.
- String-level Rfam projection corrupts labels; pair-level projection is the
  defensible construction.
- Removing or scrambling distal context during head training does not improve
  OOD F1; this rules out the head-level shortcut account, not pretraining
  contamination.
- Corpus coverage predicts per-family advantage, but the matched one-hot
  control shows a similar relation; call this a shared predictor, not a
  pretraining-specific mechanism.
- H4 is conditional: if the same PDB-clean chains score better against
  comparative labels than against structural labels, annotation convention is
  part of the learned-model advantage. Otherwise report the clean negative.
- Homologous distal context retains more than unrelated context at test time,
  but the pre-registered H3 bar was not met; report descriptively.
- Energy-consistent decoding does not reduce long-range FPs at matched budget;
  it does increase TPs, reported as exploratory.

## Claims that are not allowed

- Any Round 4 number produced with the old `pair_prf` zero-overlap convention.
- H2 as a mechanism, the 102% transfer figure, or a claim that masking distal
  context improves prediction.
- A SOTA or “first” claim not backed by a matched external baseline.
- Any test_novel number before the Amendment 14 one-shot unblinding.

## Remaining work before submission

1. H4: acquire the PDB-derived, cmscan-assigned, train-clean chain cohort; build
   fixed comparative and structural labels; run the three heads, RNAfold and
   EternaFold; report the family count and D.
2. Pre-register and run the gated shift refinement; report oracle recovery.
3. Run the clan-grouped `(s,b)` refit and A-fixed on withheld under old and
   refit parameters.
4. Add two predeclared head seeds (five total) and the default-weight MXfold2
   row; keep the short-schedule retrain clearly labelled.
5. Amendment 15 final analysis list, then the single `test_novel` unblinding;
   no test-set number may inform earlier choices.
6. Human review of the final narrative and citations; the current ten
   references were checked on 2026-11-15.
7. Refresh the T6 claim table after H4 and the unblinding.
