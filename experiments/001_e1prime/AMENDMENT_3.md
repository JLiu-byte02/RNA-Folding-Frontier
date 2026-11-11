# Amendment 3 — H2.4 restated, P>R as evidence against H2 on OOD, and the
# confirmatory/replication split

Date: 2026-11-12. Source: PI Round 4c ruling (Claude, 2026-11-12, reply to the
message archived as `docs/claude/transcripts/2026-11-12_round04c_reply.md`).
The frozen text and Amendments 1–2 stay as written; this amendment overrides them
where they differ. Committed and hashed **before any diagnostic dump or
P2-block cache existed**.

## A3.1 H2.4 restated — energetic consistency, not structural validity

The Round 4 numbers were produced by the production decoder, which **already**
enforces one partner per position and non-crossing (Amendment 2 §A2.3). So the
Round-4 "e helps" result cannot be explained by multi-partner false positives,
and the original wording of H2.4 — "bare threshold vs one-partner decoding" —
measures something much weaker than the claim it was meant to test. Restated
prediction:

> **H2.4.** Decoder A — exact nearest-neighbour DP over **loop energies plus
> learned pair corrections** — reduces long-range false positives (``|i-j| >
> 128``) relative to greedy nested matching **at matched recall**, and the
> reduction is larger than any accompanying loss of true positives.

Protocol:

1. Decode the same probability map with greedy nested matching and with decoder
   A. Choose the greedy threshold ``t_g`` by bisection so that greedy recall
   equals decoder-A recall to within 0.005; report ``t_g`` and the residual.
2. Count, per sequence, long-range false positives (LRFP) and true positives
   (TP). Statistics: ``drop_lrfp = LRFP_greedy(t_g) - LRFP_dp`` and
   ``drop_tp = TP_greedy(t_g) - TP_dp``.
3. Pass if ``drop_lrfp > drop_tp`` with a family-block bootstrap 95% CI on
   ``drop_lrfp - drop_tp`` that excludes 0.
4. Both splits: id_test and test_withheld; both reported, withheld is the
   evidential one.

Two decoder-A variants are pre-specified, and they are not interchangeable:

| variant | energies | role |
|---|---|---|
| A-struct | pairwise log-odds only, zero loop penalties | **exploratory ablation**: structural validity only |
| A-energy | ViennaRNA Turner energies + learned per-pair corrections (`sc_add_bp`) | **confirmatory**: this is the H2.4 test |

A-struct may be reported but cannot confirm or refute H2.4: the PI's point is
that greedy decoding already provides structural validity, so only energetic
consistency is new information.

## A3.2 Precision above recall on withheld is evidence against H2 on OOD data

H2 predicts over-prediction of long-range partners. CHANRG reports the opposite
for foundation models out of distribution: high precision, collapsing recall.
The two claims can coexist because our head uses ``pos_weight = 20``, but a
specific observation discriminates them:

> **A3.2.** If, on test_withheld under condition (a) at threshold 0.5 with the
> production decoder, family-macro precision exceeds family-macro recall
> (``P > R``), that is recorded as **evidence against H2 on OOD data**. The
> Round-4 (e) gain would then be read primarily as the logit shift K1 tests, not
> as removal of spurious partners.

This is a pre-committed interpretation rule, not a new test: it uses numbers D1
already produces.

## A3.3 Confirmatory claims on test_novel (frozen, embargoed)

test_novel stays unscored until the following are pre-registered **and** the
diagnostics have decided whether the mechanism is alive. When it is unsealed, the
confirmatory claims are:

- **H2.3** — P2-block rewriting changes the long-range error structure, measured
  with the native head;
- **H2.4** as restated in A3.1 (A-energy variant);
- **P5** — decoder A with FM features beats BPfold, EternaFold and RNAstructure
  on clan-novel data and on **S4** (temporal PDB structures).

Everything else that uses test_novel is secondary or exploratory and must be
labelled so.

## A3.4 H2.1 is replication, not novelty

Flamm, Wielach, Wolfinger, Badelt, Lorenz and Hofacker, "Caveats to deep learning
approaches to RNA secondary structure prediction", *Frontiers in Bioinformatics*
2:835422 (2022), doi 10.3389/fbinf.2022.835422, tests DL predictors against
synthetic data from a known generative model. Snippet-level reading suggests it
reports scaling failures and dependence on training families. **The full text
has not been read yet**, so:

- H2.1 (long-range FP rate grows with length) is reported on id_test, withheld
  and novel as **replication**, citing Flamm et al. only if the full text
  supports the specific claim;
- the exact claim and its figure/table number must be recorded in
  `docs/literature.md` before any citation is written;
- before any novelty statement about H2.4, check whether Flamm et al. or their
  follow-ups already test a **learned-score + exact-energy DP**; that is the
  closest possible overlap with our decoder-A contribution.

No novelty claim may rest on H2.1.

## A3.5 What this amendment does not change

- The metric (family-macro base-pair F1, family-block bootstrap) and the split
  definitions are unchanged.
- K1 keeps its role and its kill rule; its result is still reported first.
- The positive control (sequence-level split inside training clans) remains
  required, because a null mechanism result is only interpretable if the
  instrument can detect a shortcut when one exists.
