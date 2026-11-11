# Next Actions

Updated: 2026-11-11

## Immediate (E0)

- [x] Get PI decision on direction and first experiment batch (D1 selected).
- [x] Verify PI's literature claims independently (CHANRG, BPfold, MXfold2 DP,
      RiNALMo weights).
- [x] Acquire BPfold release (bpRNA TR0/VL0/TS0, ArchiveII, bpRNA-new,
      RNAStrAlign) and ArchiveII mirror.
- [x] Archive Round 1 (E1′) and Round 1b (E0 decisions) transcripts.
- [x] Label-discordance analysis on the 812 exact duplicates (E0b).
- [x] Rfam 15.0 metadata: 4,178 families / 147 clans / 463 memberships; ArchiveII
      clan map (E0c). SEED dump parsed to 92,781 ungapped seed sequences.
- [ ] Download Rfam 15.0 CMs (EBI intermittent) and run cmscan to map every
      sequence to an Rfam family/clan.
- [ ] Build Infernal from source (in progress; no conda, no cmsearch on PATH).
- [ ] Check BPfold's own framing of ArchiveII (OOD vs standard benchmark).
- [ ] Re-score BPfold on clan/near-dup-cleaned ArchiveII (T4/A.2).
- [ ] RiNALMo-giga crop-and-positional sanity check (PI Round-1 ask 2).
- [ ] Homolog pool counts below 60% identity per test family (condition c).
- [ ] structRFM pretraining manifest and overlap with every test set.
- [ ] Build splits S1-S4 and run the leakage audit (identity, clan overlap,
      pretraining-corpus hits).
- [ ] Record cluster facts: scratch quota, CPU cores per node, preemption.
- [ ] Decide repository license with the PI (RiNALMo-giga is AGPL; rest is MIT).

## Queued

- [ ] E1′ context-intervention factorial (pre-registered P1–P5).
- [ ] E2 clan-shortcut probe.
- [ ] E3 evidence-stratified evaluation.

