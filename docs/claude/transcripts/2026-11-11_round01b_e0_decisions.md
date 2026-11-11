# Claude Round 1b — Decisions on the E0 findings

- Date: 2026-11-11
- Model: Claude Opus 5.5, High effort
- Session: https://claude.ai/chat/cc2ad40a-406e-481e-8e7b-89a64b78f2f4
- Executor input: E0 first pass (BPfold v0.2 release, TR0 = 10,814): ArchiveII
  20.47% exact / 38.25% near-duplicate overlap with TR0; TS0/VL0/bpRNA-new clean
  at k=15; licenses verified; three questions (A/B/C).

## A) ArchiveII overlap: a section inside a broader audit, not a standalone paper

Do not promote it to a standalone contribution. Make it one panel of an
"evaluation contamination" section supporting the main thesis: reported OOD
accuracy is inflated by both shortcut learning (the mechanism) and evaluation
leakage (this audit). A single-target audit reads as a gotcha; reviewers discount
those. The systematic version covers BPfold, MXfold2, UFold, SPOT-RNA, the
structRFM SSP heads and RiNALMo fine-tuning sets at three levels: sequence, clan,
pretraining corpus.

Before any wording goes into print:

1. **Check BPfold's own framing.** If the paper calls ArchiveII a standard
   benchmark rather than a cross-family test, the finding becomes "the benchmark
   is commonly read as OOD but isn't" — weaker and fairer. If they do call it
   OOD, quote the passage.
2. **Re-score BPfold on clean ArchiveII**, removing the 1,517 near-duplicates as
   well as the 812 exact ones. The informative number is the drop from 0.820. If
   it is under 0.02 the contamination is real but does not change their result.
3. **Compare labels on the 812 exact duplicates** (cheap, new information).
   Pair-level F1 between the two label sets, overall and per family. This
   estimates label noise, bounds achievable F1, and seeds the D3 evidence
   stratification.

Use neutral wording ("evaluated without train-test deduplication") and give the
authors a courtesy note before preprinting.

## B) S1': pooled set, cmscan as the single assigner

Neither option as posed. For Rfam-derived entries, bpRNA-1m's structures are
themselves projections of the seed SS_cons, so projecting Rfam 15.0 seeds gives
labels of the same type at the same cost. The real choice is about assignment.

1. **Pool:** bpRNA-1m(90) plus Rfam 15.0 seed sequences from families absent from
   bpRNA-1m. Extra families mean more held-out clans.
2. **Family assignment:** cmscan every sequence in the pool, from any source,
   against the Rfam 15.0 CMs with `--rfam --cut_ga --nohmmonly`; map families to
   clans with the release's clan table. Do not rely on bpRNA's internal IDs — one
   uniform assigner is auditable.
3. **Labels:** keep bpRNA labels where they exist; for new seed sequences project
   SS_cons. Fix one pair policy for everything: canonical pairs plus GU; decide
   once whether lonely pairs are removed; record pseudoknots but exclude them
   from primary F1.
4. **No cmscan hit:** may never enter S1' test; may enter training only if they
   have no hit to any test clan even at relaxed thresholds (E <= 10, no GA
   cutoff).
5. **Split:** hold out whole clans, families without a clan treated as singleton
   clans. Stratify held-out clans by length and structural class (rRNA-like,
   tRNA-like, riboswitch, ribozyme, other). Target 15-25% of clans held out; freeze
   the split by committing its hash before any model sees it.

## C) cmscan before E1': required, but do not idle the GPUs

Sequence-level dedup is not enough for this hypothesis: P2 (homolog swap) and P4
(clan probe) need clan labels; unclean OOD sets confound the k-curve; and S2 needs
the test clan removed from TR0 anyway, which removes the 812 duplicates
automatically.

Can start now (no clan labels needed):

- RiNALMo crop and positional sanity check.
- Caching TR0/TS0/VL0 embeddings for every k and for intervention conditions a,
  b, e.
- T4 baseline reproductions on BPfold's splits.
- The label-discordance analysis from A.3.
- Re-scoring BPfold on clean ArchiveII (A.2).

Waits for cmscan: conditions c and d, the clan probe, S1'/S2 construction, and
all E1' training runs.

**cmscan logistics:** HTTPS from `ftp.ebi.ac.uk/pub/databases/Rfam/15.0/` with
`wget -c --tries=50 --waitretry=30`, verify against the release md5. If EBI stays
down more than 48 h, escalate; no known mirror. Rebuilding CMs from seeds with
cmbuild+cmcalibrate is a last resort. Run cmscan as a Slurm array of ~500-sequence
shards; log Infernal and Rfam versions plus shard hashes in `manifest.json`.

## Evidence-format updates

- **T3 (leakage):** add clan-overlap, no-hit count, and label-discordance F1 on
  exact duplicates.
- **T4 (reproduction):** per BPfold dataset, published F1 vs our reproduction vs
  F1 on the clean subset, side by side.
- **Next report:** T3 clan-overlap counts for S2 and S3 (original vs
  clan-filtered), how many families survive ArchiveII LOOFO after clan removal,
  and per-family homolog counts below 60% identity for condition (c).
