# Research State

Updated: 2026-11-15 14:58 +0800

This is the current authoritative snapshot. Earlier versions accumulated
superseded interim numbers; use `docs/decision_log.md`, `docs/corrections.md`,
and the experiment artifacts for history.

## Phase

Paper A is the committed deliverable: a benchmark/label audit plus
well-controlled negative mechanism results. The frozen four-method comparison
contains no head advantage over the thermodynamic tools. Round 4n was delivered
on retry and its reply is archived: the paper remains audit-led, no method claim
is made, and D2 is shelved. The authorized pre-unblinding queue is: H4
annotation-convention test, gated shift refinement, clan-grouped `(s,b)` refit,
two additional head seeds, default-weight MXfold2, then Amendment 15 and one
`test_novel` unblinding.

## Round 4n decision

- H4 starts the queue. It compares the same PDB-derived, train-clean chains
  against cmalign/SS_cons comparative labels and fixed PDB structural labels.
  Supported: D >= .03 with CI excluding 0 and same sign for the one-hot head.
  Killed: D < .01 or CI includes 0. Fewer than 40 two-label families is
  underpowered and ends the test without extension.
- The shift refinement is a secondary analysis with a fixed gate: withheld
  macro gain >= .008 with CI excluding 0 and id_test point estimate >= 0; kill
  below .004 or if the CI includes 0.
- The `(s,b)` refit must use Amendment 7's clan-grouped folds; a boundary fit
  gets one declared outward step. Only the refit values enter `test_novel`.
- Two additional seeds bring the head to five. The new seeds and the five-seed
  reporting rule must be frozen before training.
- MXfold2 default weights enter as one overlap-flagged row; no further
  retraining compute is authorized.
- Amendment 15 fixes the final `test_novel` list and supersedes, without
  editing, Amendment 14. Hash the one-pass outputs and log any rerun.

## Confirmed scientific verdicts

- **D1-core is falsified.** Replacing full foundation-model context with a
  +/-32 nt block does not lower family-macro F1 by the pre-registered margin;
  the id_test and withheld results both fail the decision rule.
- **H3 is not supported.** Homologous test-time context retains more signal
  than unrelated context, but the pre-registered retention bounds are not met.
  Report the direction descriptively, not as a mechanism claim.
- **H2 and H2.4 are retired.** Energy-consistent decoding does not remove
  long-range false positives at a matched pair budget. The exploratory TP gain
  is reported, not the old mechanism story.
- **E-C1 passes; E-C1b fails specificity.** On 409 common withheld families,
  Rfam full-region coverage has partial rho `.338 [.244,.429]` with the FM
  advantage, and RiNALMo pseudo-perplexity has rho `-.242 [-.332,-.147]`. The
  matched one-hot control has rho `.253`; `d_rho = .085 [.009,.160]` is below
  the pre-registered `.15`. Coverage is therefore a shared predictor, not a
  pretraining-specific mechanism.
- **E-C2 fails.** The gate collapsed to `s_max=.1875` with
  `alpha=beta=gamma=0`; the fitted arm is a constant scale bit-identical to its
  control, and the id_test comparison reverses sign.
- **No default-decoder advantage.** On the frozen common cohort the head's
  three-seed mean is `.5905` withheld and `.5901` id_test; every paired
  family-bootstrap interval against RNAfold, EternaFold or RNAstructure includes
  zero. Seed 17 is nominally above two tools, but the seed-level interval is
  wide and the study cannot resolve differences near `.01`.

## Frozen evaluation assets

- Common cohort: withheld **16,040 sequences / 408 families**, id_test
  **1,272 / 153**. The retrained MXfold2 600 nt cap binds at 597 / 572 nt.
- Baselines: RNAfold, EternaFold, RNAstructure and MXfold2 retrained on the
  frozen train split. The retrained MXfold2 scores `.5518/.5491` and is a
  control, not the published system.
- Labels: `labels_v3` pair-level projection parses 92,781/92,781 Rfam seed rows
  and yields 2,517,623 canonical pairs. On 37,797 exact-sequence comparisons
  against bpRNA-1m, mean pair F1 is `.9056` for v3 versus `.8148` for v2.
- Primary scoring is v3, nested truth, family-macro base-pair F1 with a
  `>=4` true-pair exclusion and family-block bootstrap. Never pool the old
  per-tool populations.
- `test_novel` remains embargoed. No scoring output from it exists. Only
  Amendment 14 can authorise the single unblinding pass.

## Audit results retained in Paper A

- ArchiveII is 61.2% exact and 88.4% near-duplicate of RNAStrAlign union bpRNA
  TR0; its clan coverage is narrow and leave-one-family-out leaks through
  tRNA/tmRNA and Rfam clan siblings.
- On the label mirror, ArchiveII labels agree with RNAStrAlign at pair F1
  `.9937` and with bpRNA TR0 at `.9195`; telomerase agreement is `.6606`.
- String-level Rfam projection is a label defect: v1 leaves 17.8% of seed rows
  unparseable, while pair-level v3 parses every row. Use v3 provenance only for
  v1/v2 comparisons.

## Exploratory diagnostics

- Corrected FP anatomy uses the deployed greedy head, not thresholded maps.
  On withheld it emits 188,158 FPs versus RNAfold's 311,768; 21.9% are within
  1-2 nt of truth versus 13.4%, and 51.6% are farther than 10 nt versus 64.3%.
  On id_test the near-miss shares are 31.1% versus 20.3%.
- A radius-2 shift oracle raises macro F1 by `.0251` on withheld and `.0263` on
  id_test. This is an upper bound that ignores crossing constraints,
  competing FPs, selection cost and model score. No shift decoder has been run
  or claimed.

## Paper A status

- Sections 1-4, Tables T1-T7, Figures 1/2/6/7, captions, availability and ten
  verified references are drafted from committed artifacts.
- Title, abstract, discussion and conclusion remain held for the PI.
- Claim 9 is a single-cohort comparison with no head advantage over the
  thermodynamic tools and a weaker retrained MXfold2.
- Every table number is checked against JSON by `tests/test_paper_numbers.py`;
  figure regeneration and family diagnostics are CI-gated.

## Operations

- Repository: `https://github.com/JLiu-byte02/RNA-Folding-Frontier`; current
  `fresh-clone` CI is green and the local suite is 187 passed / 2 skipped.
- SSH to the compute server is up. No GPU job is running.
- ZoogVPN reports `protected`, but its WireGuard peer retries handshakes every
  five seconds; system DNS is unusable while SSH by IP still works. See
  `docs/network_2026-11-15.md`.
- RiNALMo-giga is AGPL-3.0 and isolated behind `src/rff/models/backbone.py`;
  the repository licence is still undecided.

## Do not reuse

- `runs/p2_seed*` checkpoints selected under the buggy zero-TP `pair_prf`
  convention.
- Any family-macro number produced before the `pair_prf` fix and absent from
  `docs/corrections.md`.
- The old `102%` P2 transfer figure, H2 as a mechanism, or any test_novel
  performance number.
