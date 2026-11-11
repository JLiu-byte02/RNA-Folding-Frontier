# Round 3 — PI decisions (Claude Opus 5.5, effort High)

Session: https://claude.ai/chat/cc2ad40a-406e-481e-8e7b-89a64b78f2f4
Executor message: `docs/claude/outbox_round03.md` (sent as an attachment).
Date: 2026-11-11.

## Verification of our work

> "The assignment gate, the reverse-homology screen and the positional check are
> what I needed, and the self-test catching the sharder bug before the full run
> is how this should work. Your finding 4(b) is the most important item here:
> cropping degrades the representation for reasons that have nothing to do with
> family identity. I'm changing the design around it."

## Decision 1 — swaps become the primary causal evidence

> "The reduced sweep is approved: RiNALMo-micro at k in {8, 32, 128, full}, and
> giga at {32, full}. P1 (the k-curve) drops to secondary, descriptive evidence.
> Write it as a dated amendment; don't edit the frozen text."

The missing design rule, to be written into the amendment:

> "P2 (does the trained predictor rely on family identity?) is a test-time
> intervention. Train the head once on native inputs (condition a). Then, at
> evaluation only, swap the distal context to conditions b-e. Same model, only
> the input changes, which is the clean causal design."
>
> "P3/P5 (does removing distal content cause the OOD gain?) is train-and-test per
> condition. Train a separate head on each condition's embeddings and evaluate it
> on the same condition."

Swap-specific confound and its control:

> "Chimeric inputs are also out of distribution for the FM, and condition c (a
> homolog's distal context) will naturally distort the local embedding less than
> condition d (another clan's). So for every condition, record two numbers per
> position: distortion (1 - cos(swap(i), native(i))) and clan decodability (the
> probe's accuracy on the swapped embeddings). Then fit F1 ~ distortion + clan
> decodability across conditions and positions. P2 counts as supported only if
> clan decodability predicts F1 beyond what distortion explains. Treat condition
> e (pad/mask) as a lower anchor for distortion, not as a natural-input
> condition."

> "RNA-FM may use learned absolute positional embeddings; please check, because
> I'm not sure. If it does, crops must pass the original position_ids, or RNA-FM
> runs only the swap conditions, which keep positions fixed by construction."

## Decision 2 — NOHIT stays excluded; relaxed scan as audit only

> "The split is frozen, and a 6.4% gain in training data isn't worth re-freezing
> and re-auditing it. Run the 1 h relaxed scan (E <= 10, no GA cutoff) anyway,
> but only to report two counts: how many NOHIT pool sequences hit a test clan at
> the relaxed threshold; how many of the training sequences that do have hits
> also get secondary hits to test clans at E <= 10, a stricter version of your
> E <= 1 screen. If the second count is non-trivial (above about 0.5% of
> training), tell me before E1' results are interpreted."

## Decision 3 — stratify, and treat label provenance as the bigger issue

> "Two confounds need handling before anything is evaluated. (i) Label
> provenance. Training labels are all bpRNA; test_novel labels are all Rfam 15.0
> seed SS_cons projections. A gap between withheld and novel could just be a
> label-source difference. test_withheld contains both sources (5,844 bpRNA +
> 11,527 seeds), so use it to calibrate: report F1 on the bpRNA and seed subsets
> of the same withheld clans, per model."

Per-family structure source, from the SEED ``#=GF SS`` and ``#=GF TP`` lines, for
all 1,983 novel clans, before committing the tables.

Pre-registered tables (to be fixed in the amendment before any model sees
test_novel):

- Primary: clan-macro F1, bootstrapped by clan, on novel families with a
  published structure source and at least 4 base pairs in SS_cons; within each
  family deduplicate at 90% identity and keep at most 10 sequences, sampled with
  a fixed seed.
- Stratified rows: by type (cis-regulatory/riboswitch, ribozyme, snoRNA, miRNA,
  sRNA, other) and by complexity (has a multiloop or not).
- Supplementary: the full 42,873 sequences, all families, no caps.

> "If the published-structure subset turns out to be under about 150 clans,
> report back and we re-decide; don't silently relax the rule."

> "No test_novel model numbers until the amendment containing the primary-table
> definition is committed and its hash is recorded."

## Decision 4 — GPUs run the context factorial with decoder B

> "Decoder A as pre-registered is the nearest-neighbor DP with learned
> corrections, and it doesn't exist yet. A Nussinov decoder with learned scores
> isn't the matched comparison, so A vs B has to wait."

Starting now on GPUs:

1. Embedding cache, giga and micro, conditions a-e at a fixed +/-32 local window,
   on train, id_test, test_withheld and the capped test_novel sample.
2. The P2 test-time intervention: one decoder-B head trained on condition a,
   3 seeds, evaluated under a-e on id_test and withheld; the clan probe runs on
   the same cache.
3. P3 train-and-test heads for conditions a, b and d on micro first.
4. In parallel on CPU: connect MXfold2's DP so it accepts external score tensors,
   check it against Turner-only output on 100 sequences, then build decoder A.

> "Within the factorial, giga takes priority, because P2 has to hold on the
> strongest backbone."

## Next report

- P2 table: the condition-a head evaluated under a-e, for id_test and withheld,
  3 seeds, clan-macro F1 with CIs.
- Distortion and clan decodability per condition, plus the regression of F1 on
  both.
- Withheld calibration: F1 on bpRNA vs seed labels within the same clans.
- Novel-set label metadata: counts by SS source (published vs predicted) and by
  type, and the size of the primary subset.
- The relaxed-scan audit counts.
- Decoder A status, including the Turner-only check against an existing tool.
