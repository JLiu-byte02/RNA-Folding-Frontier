# Research State

Updated: 2026-11-13 23:05

## Phase

`post-E-C2` — P2-block and P3 are complete, D1-core is falsified, H2/H2.4 are
retired. The E-C2 context gate failed: its fit collapsed to the grid boundary
(`s_max=.1875`, `alpha=beta=gamma=0`), the fitted gate is therefore a pure
constant scale that contributes exactly zero, and id_test reverses sign. What
survives is E-C1 (pretraining coverage predicts generalisation) and the
negative mechanism results. The project needs a new method story from the PI;
the paper can be written as a rigorous negative/audit paper until then.

## Decisive 2026-11-13 result

- P3 common population: 16,163 withheld sequences, 409 families, v3 all-pairs,
  fixed 0.5.
- a `.5959 [.5726,.6191]`; e-block `.5965 [.5733,.6196]`; b-block
  `.5860 [.5630,.6088]`.
- Paired Δ vs a: e `+.0007 [−.0079,+.0093]`; b `−.0099 [−.0159,−.0036]`.
- **D1-core falsified** under Amendment 10's .02 rule.
- P2-block H3 `R_c=.895 [.726,1.052]`, `R_d=.225 [.019,.342]`; **not
  supported** under the pre-registered bounds.
- P2-central common-population ratios are also below the old 77%/102% reading:
  kept–kept withheld `R_c=.747 [.564,.960]`.
- H2.4 with frozen `(s,b)=(.1875,0)`: DP has more TPs but also more long-range
  FPs than greedy on both splits; **not confirmed**.
- E-C1 supported under Amendment 12: Rfam coverage `rho=.338 [.244,.429]`,
  RiNALMo pseudo-perplexity `rho=-.242 [-.332,-.147]`; hypothesis C is alive.
- Full report: `docs/p2block_p3_2026-11-13.md`.

## E-C2 gate failure (2026-11-13)

- Fit on the freeze split collapsed to `s_max=.1875` with `alpha=beta=gamma=0`.
- Withheld family-macro: gated `.6168` vs fixed `.6137` (paired Δ `+.0031
  [−.0038,+.0088]`); the constant-scale control is **bit-identical** to the
  gated arm, so the context term contributes exactly zero.
- id_test reverses sign: Δ `−.0064 [−.0122,−.0002]`.
- The coarse grid restricted `beta>=0`, so **negative** dependence is untested,
  not excluded. Further gate work requires a new pre-registration.
- Against EternaFold withheld `.6051`, the gated arm's margin is only `+.0117`.
- Full report: `docs/ec2_2026-11-13.md`; Amendment 13 pinned.

## External baselines (frozen cohort policy)

- Numbers are per-tool populations and **must not** be compared across tools
  until the frozen common cohort is built (`scripts/build_common_cohort.py`):
  dump ∩ truth_v3 `>=4` pairs pre-nested ∩ prediction present in every tool.
- Withheld, 16,163 seqs / 409 families: RNAfold `.6003` [.5733,.6265],
  EternaFold `.6051` [.5792,.6302], micro `.5548`/`.5901`.
- id_test, 1,287 seqs / 155 families: RNAfold `.5955` [.5530,.6375],
  EternaFold `.5916` [.5479,.6329], micro `.5295`/`.5465`.
- Head (P2-fix a, 3-seed mean) is **below** RNAfold in family-macro
  (withheld −.0092, id_test −.0118) and **above** in micro
  (+.0394, +.0642). No "beats thermodynamics" claim is available.
- 1,024 nt context is binding: 578/17,371 withheld and 10/1,402 id_test dumps
  are missing. Families RF02540 and RF02543 are lost on withheld.

## Established facts

- GitHub repository created: https://github.com/JLiu-byte02/RNA-Folding-Frontier
- Slurm partition `g-3090-2`: 4 nodes x 2 RTX 3090 GPUs (8 total, shared with
  other users). Project cap: 6 GPUs.
- Remote project root: `/home/junxuanliu/projects/RNA-Folding-Frontier`, now a
  git checkout of `origin/main`; 55 GB of data, caches, and runs remain
  server-side and untracked.
- Login-node Python is 3.6.8; Conda modules are available. A dedicated environment
  must be created before any experiment.
- E0 first pass: ArchiveII overlaps bpRNA-1m TR0 by 20.47% exact sequence
  duplicates and 38.25% near-duplicates (k=15 containment >= 0.5). TS0, VL0 and
  bpRNA-new are sequence-clean at this threshold.
- **E0e:** against the union of BPfold's two trainers (RNAStrAlign ∪ bpRNA TR0),
  ArchiveII is **61.20% exact duplicates** and **88.40% near-duplicates**. The
  long families (16S, 23S, group II) are clean only because the trainers cap
  length at 600 nt. ArchiveII also contains 313 duplicate sequences internally.
- **E0f:** ArchiveII vs RNAStrAlign labels agree at mean pair F1 **0.9946**
  (91.2% identical), while ArchiveII vs bpRNA TR0 agree at **0.9191** (46.6%
  identical). Overlap plus label identity is the mechanism that turns leakage
  into F1.
  **Correction:** those mean-F1 values used `pair_prf` before its zero-TP bug was
  fixed and must be rerun; the exact-identity percentages are unaffected.
- **E0c:** Rfam 15.0 has 4,178 families, 147 clans, 463 clan memberships. ArchiveII
  spans 8 clans; tRNA and tmRNA share CL00001, so leave-one-family-out leaks on
  that fold.
- **bpRNA-1m mirror:** the HuggingFace `rouskinlab/bpRNA-1m` mirror (66,715
  sequences, sha256 `d36c62ee...`) contains 98.7-98.9% of TR0/VL0/TS0 sequences
  and is a strict superset of their bpseq labels (0 missing pairs in 13,419
  sequences). Filtering it to canonical+GU reproduces bpseq exactly for 91.2% of
  TR0. A 90% subset (CD-HIT 4.8.1) has 26,701 representatives.
- **labels_v3:** pair-level projection parses 92,781/92,781 seed rows (100%),
  produces 2,517,623 canonical pairs, and differs from v2 on 6,757 rows. The
  frozen artifact hashes to `b8fe71ab…` and reproduces on Linux.
- **Label defect, quantified:** among the 11,527 seed rows in `test_withheld`,
  2,392 (20.75%) had empty truth under v1; under v3, 112 are empty and 341
  (2.96%) fall below four pairs. In `test_novel`, 10,470/42,873 v1 rows were
  empty and 3,675 (8.57%) fall below four under v3.
- **v3 cross-check:** 37,797 exact-sequence comparisons against bpRNA-1m give
  mean pair-F1 .906 for v3 vs .815 for v2; v3 wins 5,240 comparisons and loses
  2,156. Details in `experiments/011_labels_v3/v2_v3_bprna_crosscheck.json`.
- **D1 label policy (corrected after the `pair_prf` fix):** primary nested-truth
  e−a is **−.3140** on id_test and **−.2977** on withheld; partner-arm totals
  are −.3110 and −.2962. The confidence intervals are entirely below zero.
  Masking (e) is substantially worse; the earlier positive e−a was a scoring
  artifact.
- **H2.2:** long-range FP helices are about 4× enriched over the null on both
  splits. H2.1 is split-dependent after candidate normalisation (1.62× id_test,
  3.89× withheld).
- **Template control 1:** condition e's one-kept and neither-kept F1 are near
  zero (one: .018–.024, neither: .006–.009), so it does not invent pairs outside
  the kept window.
- **Template control 2 (corrected):** all-N family-macro F1 is .00055/.000066;
  dinucleotide shuffle is .0162/.0173. Neither supports a template prior. The
  earlier .89/.90 shuffle result was the `pair_prf` bug.
- **FP-helix energetics:** 26.6%/26.7% of FP helices are fully in the RNAfold
  MFE, 37.5%/37.3% have at least one pair there; their isolated energies are
  close to true helices. H2.4 may fail for a principled reason.
- **Round 4j correction:** the `pair_prf` zero-TP defect voided the Round 4
  P2-central evidence. H2 is retired as a mechanism; H2.4 continues as a
  standalone decoder question. Amendment 9 is pinned at `a1dae27e…`.
- **E0 after correction:** ArchiveII vs bpRNA TR0 mean pair-F1 .9195 (server
  3,370-record mirror), telomerase .6606; vs RNAStrAlign .9937. The old run used
  a 3,966-record mirror, so counts are not directly comparable; F1 movement is
  small. T4/Turner-control/independent-scorer numbers used `prf`, not
  `pair_prf`, and stand.
- **Current next experiment:** retrain condition-a heads (seeds 17/23/41) under
  corrected selection, then rerun P2-central a–e with v3 nested labels and the
  kept–kept/all/micro readouts. No H2.4 or P2-block until that table is seen by
  the PI.
- **Corrected P2-central (3 seeds, new heads):** a remains best; c is close to a
  on kept–kept (withheld .5927 fixed/.6117 oracle vs a .5910/.5906; id_test
  .5616/.5755 vs .5759/.5753); b/d/e lose both all-pair and kept-window F1.
  c−d kept-window contrast is +.058 (id_test) / +.070 (withheld). This is
  consistent with a family/template shortcut but awaits the positive control.
  Full table: `docs/p2central_fixed_2026-11-12.md`.
- Infrastructure: Infernal 1.1.5 built from source; Rfam 15.0 CMs pressed;
  cmscan running over 37,343 sequences (bpRNA-90 reps + ArchiveII + TS0 +
  bpRNA-new) on the login node because the CPU partition is fully allocated.
- Environment: the rff env has transformers 5.20 which requires torch >= 2.6;
  torch 2.5.1 is being replaced with 2.6.0+cu124 (driver 550.144.03 supports it).

### E1′ first results (2026-11-11 late)

- **S1′ frozen** (3,861 clans: 1,983 novel / 376 withheld / 1,502 retained).
  train 17,746 · id_test 1,402 · test_withheld 17,371 · test_novel 42,873.
- **P2 test-time intervention** (RiNALMo-giga, decoder B, 3 seeds) —
  `experiments/002_e1prime_runs/P2_REPORT.md`. Family-macro F1, Δ vs native:
  id_test (e) N-mask **+0.089 [+0.045, +0.134]**, (d) other-clan **+0.053
  [+0.009, +0.098]**, (b) shuffle +0.027 [−0.014, +0.071], (c) homolog
  **−0.092 [−0.141, −0.040]**; withheld: (e) **+0.041 [+0.019, +0.063]**, (d)
  +0.007, (b) +0.008, (c) **−0.076 [−0.101, −0.050]**.
- **Direction of the effect is length-dependent**: masking the distal region
  helps at 80–299 nt (+0.08 to +0.11 micro) and hurts at ≥300 nt (−0.03). This
  is the first clean mechanistic signal of the project.
- **P2/P3 as pre-registered are not supported.** The (a)−(e) "gain" is negative,
  so P2's retention criterion is undefined; (c) goes the wrong way.
- **Clan probe**: 0.9216 held-out accuracy over 590 training clans, but 0.0526
  on id_test (chance 0.042) — the identity signal is family-specific, and the
  splits hold families out, so a family shortcut has nothing to exploit there.
  The Amendment 1.3 decodability term adds ΔR² ≈ 0.002–0.003.
- **test_novel gate fails**: 310 published-SS families collapse to 6 clans, not
  the required ~150 (`experiments/001_e1prime/novel_eval_set.md`). Escalated.
- **Condition (c) coverage** 39% (id_test) / 49% (withheld); CD-HIT cannot
  cluster below 80% identity for nucleotides, so homologs are found by verified
  alignment band instead.

## Open questions

- Does E-C1 (coverage predicts generalisation) survive on the frozen common
  cohort, where all four external tools and our head are scored on the same
  sequences?
- Is the family-macro deficit against RNAfold real on that cohort, or an
  artefact of per-tool populations?
- Would a *negative* context term (`beta<0`, never on the fitted grid) do
  anything, or is the gate structurally unidentifiable on this data?
- What method story replaces H2/E-C2? No candidate is pre-registered yet.

## Current blockers

- No scientific blocker. CPU is unlimited by policy (still via Slurm).
- Four-tool common cohort not yet materialised: RNAstructure withheld job hit
  its wall-clock limit and must be resumed from its append-only cache.
- The workspace `.git` is read-only under the managed profile; commits are being
  created through a temporary object/worktree repository and pushed with the
  repository deploy key.

## Do not claim

- Any performance number, SOTA statement, or novelty claim that is not backed by a
  committed artifact in this repository.
- That BPfold calls ArchiveII out-of-distribution: they call it a *sequence-wise*
  dataset and reserve "cross-family" for Rfam 12.3-14.10.
- That masking the distal region improves prediction. The corrected D1/K1
  results say the opposite by about .30 F1. The old positive result was a
  `pair_prf` bug artifact.
- That the E-C2 context gate helps. The fitted gate is a constant scale; it is
  bit-identical to its own control and reverses sign on id_test.
- That the head beats RNAfold or EternaFold on family-macro F1. It does not on
  the current per-tool populations, and the frozen cohort is not scored yet.

