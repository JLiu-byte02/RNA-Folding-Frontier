# Next Actions

Updated: 2026-11-15 15:04 (A15A pinned; H4 cohort construction active)

Governing documents: `experiments/001_e1prime/PREREGISTRATION.md`, Amendments
1-14, the Round 4n ruling in
`docs/claude/transcripts/2026-11-15_round04n_reply.md`, and
`state/amendments.lock.json`. Amendment 8 fixes `labels_v3` as primary and
secondary labels; 10/11 pre-register P2-block/P3/D1-core; 12 freezes E-C1; 13
freezes E-C2; 14 pins the one-shot `test_novel` list.

`AMENDMENT_15A.md` is the pre-unblinding protocol and must be pinned before the
new seeds, H4, shift refinement, or refit are run. `AMENDMENT_15.md` remains
reserved for the final `test_novel` list, as the PI specified; A15A must not edit
or supersede A14.

## Immediate queue

1. **[done] Pin A15A before any new result.** It froze:
   - new head seeds **43 and 47**, giving five seeds total;
   - the five-seed reporting rule (family-bootstrap CI and across-seed t
     interval; no post-hoc seed selection);
   - H4's population, two label constructions, predictors, D statistic,
     supported/killed/underpowered gates, and minimum family count;
   - the exact shift-refinement rule, `.008` pass / `.004` kill thresholds,
     oracle-recovery reporting, and three-seed evaluation;
   - the clan-grouped `(s,b)` refit, one-step outward rule at the boundary,
     and the requirement to report old and refit parameters;
   - the default-weight MXfold2 row and its overlap annotation.
2. **H4 cohort (CPU/data, current task).** Obtain bpRNA-1m PDB-sourced chains
   outside train; exclude k=15 containment >= 0.5 against train; cmscan-assign
   to an Rfam family; cmalign each chain to its family CM; build pair-level
   comparative labels and fixed PDB structural labels; apply canonical/GU,
   maximum-nested, >=4-pair policy to both; run p2fix 3-seed, one-hot,
   RNAfold and EternaFold; report family count and bootstrap D.
3. **Shift refinement (CPU; frozen cohort).** Implement the exact helix-level
   register-shift rule from Round 4n, evaluate all three seeds on withheld and
   id_test, and report the gate verdict, oracle recovery, and near/far FP
   changes. No looser threshold may be substituted after seeing results.
4. **Clan-grouped `(s,b)` refit (CPU; frozen train).** Reproduce Amendment 7's
   clan-grouped folds, extend one step outward only if at the boundary, then
   rescore A-fixed on withheld under old and refit parameters. Only the refit
   enters `test_novel`.
5. **Additional seeds 43 and 47 (GPU; frozen recipe).** Train and select under
   the same corrected metric and split, then produce the five-seed head table.
   No seed may be dropped after looking at its test result.
6. **Default-weight MXfold2 (CPU; external row).** Run the released/default
   weights on the common cohort, flag training-data overlap, and compare with
   the short-schedule retrain. No more MXfold2 training.
7. **Amendment 15 and one unblinding.** After 1-6 are frozen, write A15 with the
   final list from Round 4n, pin its hash, then read `test_novel` exactly once.
   Hash every output; log any bug repair in `docs/corrections.md`.

## Closed this cycle

- [x] Round 4n delivered on retry; full reply archived.
- [x] SSH and Claude DNS recovered; the first send was confirmed absent from
      history before retry.
- [x] MXfold2 short-schedule retrain and the four-method common cohort frozen.
- [x] FP anatomy corrected to the deployed head; regression test added.
- [x] Paper A Introduction/C1 coverage specificity wording corrected.
- [x] Research state, Paper A outline, decision records and task queue refreshed.
- [x] Full suite currently 187 passed, 2 skipped; the latest fresh-clone CI is
      green.

## Standing

- `test_novel` remains embargoed until Amendment 15 is pinned. Do not read or
  score it, even for a “quick check”, before that one-pass unblinding.
- No method claim: the paper is audit-led; report parity, not FM failure.
- One common cohort only; never compare per-method subsets.
- The default-weight MXfold2 row is descriptive and must not be used to rescue
  the short-schedule retrain.
- RiNALMo-giga is AGPL-3.0; the repository licence remains undecided.
- No credentials enter the repository.
