# Decision Log

Append-only. Each entry: date, decision, evidence, alternatives, consequences.

## 2026-11-11 — Project bootstrap

- **Decision:** Start a clean, independent project (`RNA-Folding-Frontier`) rather
  than extending prior RNA_ThermoFusion work.
- **Rationale:** User requirement; avoids inheriting prior design choices and
  conclusions.
- **Consequences:** All baselines, splits, and claims must be rebuilt from scratch.

## 2026-11-11 — Compute envelope

- **Decision:** Use Slurm partition `g-3090-2`, max 6 GPUs concurrently, remote
  root `/home/junxuanliu/projects/RNA-Folding-Frontier`.
- **Rationale:** User-specified resource limit and server policy.
- **Consequences:** Favor frozen backbones, small adapters, and parallel short runs.

## 2026-11-11 — Round 0 direction: D1 selected

- **Decision:** Pursue **D1 — local evidence, global physics**: frozen FM
  embeddings on local crops -> additive corrections to a nearest-neighbor energy
  model -> Zuker/McCaskill DP. D3's evidence-stratified evaluation is adopted as
  standard protocol; D2 is the planned pivot.
- **Rationale:** Only D1 makes a causal, falsifiable claim about why FMs fail
  out of distribution (context length as a shortcut variable) and fixes it with a
  structural change rather than a new head.
- **Evidence status:** PI's literature claims are provisional and must be
  independently verified before E1 (see `docs/literature.md`).
- **Consequences:** E0 leakage audit is a hard prerequisite; no E1 number is
  trusted before T3 exists.

## 2026-11-11 — Round 1: D1 retained as a causal paper (option A)

- **Decision:** Keep D1, recast as "why do foundation models fail OOD, and what
  fixes it". Reject B (reactivity-only supervision, kept as pivot) and C
  (calibration on structRFM, rejected as incremental on a contaminated backbone).
- **Condition:** mechanism (distal context as family-identity shortcut, shown by
  intervention) plus payoff (crop + exact DP beats or complements BPfold on
  clan-held-out data). Mechanism alone caps the venue at NAR/Genome Biology.
- **Design:** E1′ = k-sweep {8..256, full} × 3 backbones × {DP, ResNet}, plus five
  distal-context swap conditions at fixed ±32, plus a clan probe. Pre-registered
  P1–P5 with thresholds; falsification and go/no-go rules recorded in
  `docs/claude/transcripts/2026-11-11_round01_decision_and_E1prime.md`.
- **Consequences:** backbones are RiNALMo-giga (AGPL, isolated module), RNA-FM
  (MIT), structRFM (MIT, contamination-flagged). Primary length scope L ≤ 512.

## 2026-11-11 — Round 1b: audit scoping, split construction, cmscan required

- **Decision (A):** the ArchiveII overlap is one panel of a broader evaluation
  contamination section, not a standalone contribution. Before wording: check
  BPfold's own framing, re-score BPfold on clan/near-dup-cleaned ArchiveII, and
  report label discordance on the exact duplicates.
- **Decision (B):** S1′ is pooled (bpRNA-1m(90) + Rfam 15.0 seeds from families
  absent from bpRNA-1m) with cmscan as the *single* family assigner
  (`--rfam --cut_ga --nohmmonly`), clans from the release table, whole clans held
  out (singletons treated as their own clan), 15–25% of clans held out, split
  hash frozen before training.
- **Decision (C):** cmscan is required before E1′. Meanwhile: RiNALMo positional
  check, embedding caches for conditions a/b/e, T4 reproduction, label
  discordance, BPfold re-scoring.
- **Consequences:** GPU work is not blocked on EBI; clan-dependent work is.

## 2026-11-11 — E0b result: identical sequences have different labels

- **Finding:** of the 812 ArchiveII sequences that are exact duplicates of bpRNA
  TR0 sequences, only 378 (46.6%) carry identical annotations. Pair-level F1
  between the two label sets: mean 0.919, median 0.982, p10 0.710. By family:
  tRNA 0.997, srp 0.970, grp1 0.965, RNaseP 0.934, 5s 0.877, tmRNA 0.796,
  telomerase 0.661. bpRNA has fewer pairs in 394/812 cases (mean precision 0.966
  vs recall 0.886).
- **Interpretation:** label disagreement, not model error, caps F1 near 0.92 on
  the overlapping subset and near 0.66 for telomerase. This is an independent
  argument for evidence-stratified evaluation (D3) and it changes how the
  contamination story is told: the duplicated subset is not merely leaked, it is
  ambiguously labelled.
- **Consequences:** every reported F1 on ArchiveII must be paired with the label
  ceiling for the corresponding subset.

## 2026-11-11 — E0c result: ArchiveII's nine families are eight clans

- **Finding:** Rfam 15.0 has 4,178 families, 147 clans, 463 clan memberships.
  ArchiveII maps to 5s→CL00113, 16s→CL00111, 23s→CL00112, RNaseP→CL00002,
  srp→CL00003, telomerase→CL00004, tRNA→CL00001, tmRNA→CL00001, grp1/grp2 →
  singletons.
- **Consequences:** leave-one-family-out over ArchiveII leaks tRNA↔tmRNA. S2 must
  remove the whole clan from the training pool; the tRNA/tmRNA pair is the
  specific fold where a naive implementation leaks.


## 2026-11-11 — E1′ P2 result: the prediction is inverted (PI Round 4)

- **Finding (three seeds, RiNALMo-giga, decoder B, family-macro F1, paired family
  bootstrap 10k):** masking the whole distal region (condition e) *improves*
  family-macro F1 by +0.089 [+0.045, +0.134] on id_test and +0.041 [+0.019,
  +0.063] on test_withheld. A low-identity same-family homolog (condition c) is
  the worst condition on both splits: −0.092 [−0.141, −0.040] and −0.076
  [−0.101, −0.050] against its own subset baseline. Other-clan (d) is +0.053
  [+0.009, +0.098] on id_test and +0.007 [−0.014, +0.028] on withheld;
  dinucleotide-shuffled (b) is +0.027 [−0.014, +0.071] and +0.008 [−0.012, +0.029].
- **P2 is not met.** The retention criterion is undefined because the (a)−(e)
  gain is negative.
- **P3 is not met.** Family-macro +0.008 with a CI that includes 0; micro −0.062.
- **Amendment 1.3 is not supported**, with the caveat that the clan-decodability
  covariate has no signal on these splits: the probe scores 0.9216 on held-out
  training sequences over 590 clans but 0.0526 on id_test (chance 0.042) and is
  undefined on withheld, because whole families inside a clan are held out.
- **Status:** the family-identity shortcut, as framed, is not supported by the
  first causal test. Everything that follows is post hoc and exploratory until
  confirmed on data not yet touched (test_novel, S4).
- **Not yet claimable:** that distal context is a "net-negative prior below
  300 nt". Two artifacts must be ruled out first — a threshold/calibration shift
  (the head was trained on native inputs and tested on swapped ones at a fixed
  0.5 cut) and the scoring geometry (one central ±32 window, not a window per
  position). The replacement candidate is H2, "spurious partners": distal
  sequence that is locally complementary to the kept region yields confident
  false long-range pairs that the 2D head cannot reject. H2 is pre-registered
  only when its diagnostics are in.
