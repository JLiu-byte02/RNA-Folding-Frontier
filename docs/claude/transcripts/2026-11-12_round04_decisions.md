# Round 4 and 4b — how the inverted P2 result is to be read, and H2

Date: 2026-11-12 (local, +0800). Model: Claude Opus 5.5, effort High.
Condensed record of both PI turns plus our factual replies. Numbers quoted here
are the P2-central table already committed in `docs/decision_log.md`; no new
measurement is claimed.

## Round 4 (PI)

The PI's verdict on the P2 result, in order:

1. **Record the inversion; do not explain it away.** Condition (e) beating (a)
   is a negative result for the direction as originally framed. It becomes
   *post hoc / exploratory* material from here on.
2. **Q1 "not yet".** The "distal context is a net-negative prior" claim may not
   be made until two artifacts are ruled out: (i) threshold/calibration shift,
   (ii) scoring geometry (one central window rather than one window per
   position). Diagnostics come before any new training.
3. **Q2 add a positive control.** A sequence-level split inside training clans
   (~80% identity, near-duplicates allowed) must reproduce the shortcut the
   instrument is supposed to detect — predicted: (c) holds up relative to (d)
   there, the reverse of S1′. Without that, absence of the shortcut on S1′ is
   not evidence.
4. **Q3 unit for test_novel.** Primary metric = family-macro F1 bootstrapped
   over **clan blocks** (singleton families are their own clan). Robustness =
   cluster test families by test-vs-test cmscan hits (E ≤ 1) and bootstrap over
   homology blocks.
5. **Q4 order of work.** D1/D2 diagnostics → train-and-test heads (giga a / e /
   crop-32; micro paused) → positive control → external baselines (BPfold,
   RNAstructure, EternaFold, MXfold2) → decoder A. The relaxed-threshold scan
   audit is non-gating.

## Round 4b (PI, after two factual replies)

1. **P2-central is a protocol deviation, not P2.** Round 1 specified a
   per-position ±32 nt window; the implementation keeps one central window and
   rewrites everything outside it. Consequences the PI made binding:
   the table is renamed *P2-central (exploratory)*; only the **kept–kept** cell
   of a region stratification is a valid readout, because a head that sees only
   frozen embeddings cannot predict a pair that touches a rewritten position;
   and the rewritten fraction runs from 19% at L=80 to 78% at L=300, so bands
   are not the same manipulation.
2. **Define P2-block and run that instead.** 16-nt blocks; for each block keep
   the block plus ±32 nt native and rewrite the rest; one forward pass per
   block; take embeddings only for that block's positions. ≈20 GPU-h total.
3. **Extraction rule must be stated before dumps.** Our reply: `prob ≥ 0.5`,
   `j - i > 3`, then greedy accept in decreasing probability under
   one-partner and non-crossing. So the primary decoder is threshold +
   one-partner + non-crossing, not a bare threshold; hence H2.4 needs a
   bare-threshold comparison arm, which `d1_scoring.py` produces.
4. **H2 is pre-registered before the dumps** (H2.1–H2.4 + kill criterion K1, see
   `experiments/001_e1prime/AMENDMENT_2.md` §A2.4). If K1 fires — oracle
   threshold plus one-partner decoding gives `e - a < 0.01` with a CI covering 0
   in the 80–299 band — the "distal context is a net-negative prior" reading is
   dropped and H2 is not pursued.
5. **test_novel stays embargoed** until H2's status is settled.

## Facts we supplied that the PI accepted

- The P2 harness stored only F1 per sequence, so D1/D2 need a fresh dump pass
  (`e1prime_dump.py`), one dump per seed.
- `test_novel` currently has 310 published-structure families but only 6 real
  clans (304 singletons); the clan-block unit in Round 4 Q3 keeps it above the
  ~150-unit gate.
- The clan-decodability probe is family-specific: 0.9216 within held-out
  training clans, 0.0526 on id_test (chance 0.042), undefined on withheld.
