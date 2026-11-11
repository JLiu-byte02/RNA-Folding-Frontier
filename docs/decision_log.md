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

## 2026-11-12 — P2-central recorded as a protocol deviation; H2 pre-registered (PI Round 4b)

- **Protocol deviation.** Round 1 specified a per-position ±32 nt window. The
  implementation kept one central window `[L//2-32, L//2+32]` and rewrote
  everything outside it. The Round-4 table is therefore renamed
  **P2-central (exploratory)**; only its kept–kept cell is a valid readout, and
  cross-band comparisons are not comparable (`1 - 65/L`: 19% rewritten at L=80,
  78% at L=300). No P2-central number may be reported as a test of P2.
- **P2-block replaces it** (pre-registered, Amendment 2 §A2.2): 16-nt blocks,
  each embedded with its own ±32 nt native flank, one pass per block, ~20 GPU-h.
- **Decoder fixed in writing before dumps** (§A2.3): symmetrise → sigmoid →
  `prob >= 0.5` and `j-i > 3` → greedy accept in decreasing probability under
  one-partner and non-crossing. D1 additionally reports a bare-threshold arm so
  the value of the global constraint is measured.
- **H2 "spurious partners" pre-registered** (§A2.4) with a kill criterion K1:
  if oracle-threshold plus one-partner decoding gives `e - a < 0.01` with a CI
  covering 0 in the 80–299 band, the "distal context is a net-negative prior"
  reading is dropped and H2 is not pursued.
- **test_novel stays embargoed** until H2's status is settled (§A2.5).
- **Order of work (PI Round 4 Q4):** D1/D2 diagnostics → train-and-test heads
  (giga a/e/crop-32, 3 seeds; micro paused) → sequence-level positive control
  (~80% identity inside training clans) → external baselines (BPfold,
  RNAstructure, EternaFold, MXfold2) → decoder A. Relaxed scan audit non-gating.
- **Amendment 2 hash:** file `experiments/001_e1prime/AMENDMENT_2.md`,
  sha256 `4a9d2f482cf12ad0ce28327bcc7b9d3b521ca24185953ddea35ef9987029de78`,
  written before any D1/D2 dump existed.

## 2026-11-12 — SSH failure is the aTrust SDP line, not sshd

- **Root cause (evidence in `docs/network_diagnosis_2026-11-12.md`).** Every SSH
  attempt is intercepted by aTrust (`flowSrc: wfpRedirect`) and accepted locally
  in ~0 ms, then fails upstream with `ec:1 "no proxy line"` and `tcpConnTime =
  int64 max`. The 0-byte banner is the local proxy, not sshd.
- **Why the line is gone:** the single WAN line `218.18.235.189:441` answered in
  46 ms at 18:34, 1055 ms at 21:42, and returned `i/o timeout` from 23:59 on;
  the client then marked all lines out and tore down the two live sessions at
  00:01:18.
- **Not our path:** `check-host.net` probes from Indonesia, Iceland and two US
  nodes all time out on `:441` while all four reach `:443` on the same host in
  236–288 ms. The host is up; the line port is closed internet-wide.
- **Not fixable locally without a credential.** The client caches the line list
  from login; only a re-login re-fetches it. No route, firewall or VNIC change
  can create a line the WAN cannot reach.
- **Tooling:** `remote.py doctor` now reports `suspected_local_interception`
  when the handshake completes in <5 ms, so this failure mode is no longer
  misreported as an sshd problem.

## 2026-11-12 — offline work while the tunnel is down: H2.3 and decoder A

Tunnel still closed (`Test-NetConnection 218.18.235.189 -Port 441` = False at
02:00). No GPU work is possible, so the two remaining pieces of the H2 measurement
that do not need the server were built and tested locally.

- **H2.3 (`scripts/d2b_block_boundary.py`).** Implements the pre-registered rule
  unchanged — extra FPs under (c) vs (a) that are not mutually native under the
  block geometry — and adds what the rule alone cannot say: a distance-stratified
  view with a distance-matched null. Reason, recorded before any P2-block number
  exists: `|W(i)| <= block + 2*flank`, so the 60% rule can pass purely because the
  extra FPs are long-range, which is H2.1's claim. The 60% rule stays the verdict;
  the stratified ratio is the mechanism check.
- **Decoder A (`src/rff/decoders/nn_dp.py`).** Exact nearest-neighbour DP with
  traceback, verified against exhaustive enumeration of all structures for
  `n <= 12` (value and traceback) under random scores with and without loop
  penalties. Three real bugs were caught this way: interior offsets with
  `dk = 0`/`dl = 0` (inner pair not strictly inside) inflated the optimum;
  `Wtb` defaulted to "paired" and looped in traceback; and skipping the cell when
  `pair_score` was disallowed also skipped the multiloop/exterior recursions,
  which made any sparse map decode to the empty structure.
- **H2.4 (`scripts/h24_decoder_compare.py`).** Compares bare threshold, greedy,
  and DP at fixed 0.5 with a family-block bootstrap on
  `drop_lrfp - drop_tp`. Zero loop penalties are the default on purpose: any
  penalty is a knob that could be tuned until the contrast appears.
- **Expectation to check, not a claim:** a long-range pair only conflicts with
  pairs that interleave its endpoints; nested pairs coexist, so DP may keep the
  FP and recover the local TPs. A null H2.4 result is informative and should be
  reported as such.
- Suite: 68 tests pass. Nothing in this entry is a result — no dump, no model
  output, and no number from id_test/withheld exists for H2.3 or H2.4 yet.

## 2026-11-12 02:12 — SSH: both WAN SDP addresses closed, portal healthy

- The 01:20 re-login *did* fetch a new line list. It has three lines, one of
  which the earlier diagnosis never saw: `vpn.suat-sz.edu.cn` →
  `218.18.235.188:441`, alongside `218.18.235.189:441` and a LAN line. So
  "controller publishes zero lines" was wrong; the correct statement is that all
  three fail detection within six seconds of the reload.
- `:441` is closed on both WAN addresses (15 s timeouts), while `:443` answers on
  both and `https://vpn.suat-sz.edu.cn/` returns HTTP 200 with the aTrust login
  page. Gateway healthy, line service unreachable — still school-side.
- The client is sticky: after `all line out` at 01:21:36 no further detection ran
  (healthy-period detection was triggered by network change/latency, not a
  timer). Recovery therefore needs a tray disconnect/reconnect once the school
  restores the port, not patience.
- Added `scripts/tunnel_watch.py --once` (probe + `remote.py doctor`, exits 3 while
  closed, 2 if the line is back but SSH fails, 0 on full recovery). Registration
  with Task Scheduler is left as an explicit step.

## 2026-11-12 03:05 — PI Round 4c: H2.4 restated, novelty split, no hand-written Turner

- **Round 4c rulings** (Claude Opus 5.5 High; archived at
  `docs/claude/transcripts/2026-11-12_round04c_reply.md`): the production
  decoder already enforces one partner and no crossings, so H2.4 as written was
  a minor contrast; it is restated as an **energetic** claim at **matched
  recall**. P > R on withheld under (a) is pre-committed as evidence against H2
  on OOD data. test_novel is reserved for H2.3, restated H2.4 and P5; H2.1 is
  replication only. Hand-written Turner is rejected in favour of ViennaRNA soft
  constraints (`sc_add_bp`), with our verified DP kept as the oracle.
- **Amendment 3** committed with those changes: sha256
  `7e035123b9bf43abbdf133a77e3bf8eed31064eabaa754185b15805a04a7962c`.
- **Flamm et al. 2022 read in full** (Europe PMC XML). Verified: abstract says
  pair counts scale quadratically with length; **Figure 3** shows it on
  matrix entries > 0.5 before postprocessing, 2,000 sequences per bin over
  30–250 nt. Their postprocessing is highest-entry-per-row + round + remove
  pseudoknots — no energy term. Discussion says the behaviour "can easily be
  addressed during postprocessing" but "it is not clear whether that would
  correct or merely hide the underlying problem", and they never run an
  energy-based decoder. Recorded in `docs/literature.md`. Consequence: H2.1 is
  anticipated (replication), H2.4 answers a question they pose but leave open.

## 2026-11-12 03:40 — H2.4 implemented at matched recall; one design caveat found

- `scripts/h24_decoder_compare.py` now emits two readouts: `at_0_5`
  (exploratory, old wording) and `matched_recall` (confirmatory, A3.1). The
  greedy threshold is chosen by grid search, not bisection.
- **Why not bisection:** greedy nested decoding is monotone in the threshold
  after all (adding lower-probability candidates cannot evict accepted pairs),
  so a bisection would work — but the *achievability* question is separate and
  matters more. A pair that blocks true positives necessarily outranks them, so
  raising the threshold to remove it removes the true pairs too. When no
  threshold reaches decoder A's recall, matching is impossible and the script
  reports `H2_4.status = "inconclusive_no_recall_match"` instead of a verdict.
  This is recorded in the module docstring and covered by a test. It is an
  implementation note, not a change to A3.1; the PI must say what to do if it
  fires on real data (candidate fallbacks: matched pair budget or matched
  precision).
- Tests: 79 pass, including the new matched-recall and unmatchable cases.
