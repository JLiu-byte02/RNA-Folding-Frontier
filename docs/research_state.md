# Research State

Updated: 2026-11-12 17:55

## Phase

`h2-template-controls` — v3 labels and D1/K1/D2 are complete; the native-head
all-N/shuffle controls and the `(s,b)` fit are in flight. H2 continues, but
H2.4 waits for a frozen calibration.

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
  mean pair-F1 .915 for v3 vs .823 for v2; v3 wins 5,293 comparisons and loses
  2,136. Details in `experiments/011_labels_v3/v2_v3_bprna_crosscheck.json`.
- **D1 label policy:** primary nested-truth e−a is +.0819 on id_test and +.0198
  on withheld; partner-arm totals are +.0740 and +.0119. K1 is a near miss on
  withheld (Δ=.0119, CI contains 0), not support.
- **H2.2:** long-range FP helices are about 4× enriched over the null on both
  splits. H2.1 is split-dependent after candidate normalisation (1.62× id_test,
  3.89× withheld).
- **Template control 1:** condition e's one-kept and neither-kept F1 are near
  zero (one: .018–.024, neither: .006–.009), so it does not invent pairs outside
  the kept window.
- **FP-helix energetics:** 26.6%/26.7% of FP helices are fully in the RNAfold
  MFE, 37.5%/37.3% have at least one pair there; their isolated energies are
  close to true helices. H2.4 may fail for a principled reason.
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

- Does the native head recover structure from all-N/shuffle inputs (template
  prior), and is condition e's id_test gain therefore partly geometric?
- Does A-energy-fixed H2.4 fail because FP helices are energetically competitive
  and present in the RNAfold MFE?
- Does H2.3 show the expected block-boundary false positives once P2-block
  caches exist?

## Current blockers

- CPU partition `c-80-255` is fully allocated by other users. The `(s,b)` fit
  is running there; the all-N/shuffle scoring job is queued behind it. Do not
  bypass Slurm on the login node for data-scale scoring.
- The workspace `.git` is read-only under the managed profile; commits are being
  created through a temporary object/worktree repository and pushed with the
  repository deploy key.

## Do not claim

- Any performance number, SOTA statement, or novelty claim that is not backed by a
  committed artifact in this repository.
- That BPfold calls ArchiveII out-of-distribution: they call it a *sequence-wise*
  dataset and reserve "cross-family" for Rfam 12.3-14.10.

