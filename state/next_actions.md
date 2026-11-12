# Next Actions

Updated: 2026-11-12 19:15 (after Claude Round 4j correction)

The governing documents are `experiments/001_e1prime/PREREGISTRATION.md`,
Amendments 1–8, and `state/amendments.lock.json`. Amendment 8 supersedes the
label clause of Amendment 7: Rfam seed primary **and** secondary labels are
`labels_v3`; v1/v2 are provenance only.

## Completed this round

- [x] `labels_v3` pair-level projection: 92,781/92,781 rows parse, sha256
      `b8fe71ab…`, full map in `experiments/011_labels_v3/`.
- [x] Eval-subset label audit: v1 463 empty-truth / v3 0 empty / v3 0 `<4`;
      Amendment 8 exclusion drops 0/2,372 records.
- [x] Amendment 8 written, hashed `5411a050…`, committed before any primary F1.
- [x] Fresh-clone CI added; local fresh clone verified all locks against blob
      bytes and passed 151 tests (1 skipped).
- [x] All 87 old-greedy failures stored as a replayable regression fixture.
- [x] Round 4h reply archived; Round 4 numbers relabelled secondary-style.
- [x] SSH restored and authenticated; the server project is now a git checkout
      of `origin/main`, with a 395 MB non-data backup and no data loss.
- [x] Gzip artifacts are binary-flagged and byte-reproducible on Windows and
      Linux; both frozen label hashes reproduce on the server.
- [x] Server-side per-split counts: withheld has 2,392/11,527 v1 empty seed rows
      (20.75%); v3 primary excludes 341 (2.96%). `test_novel` excludes 8.57%.
- [x] P2 probability dumps for conditions a/e, 3 seeds, both splits completed;
      each dump carries v1 and v3 truth.

## Running now

- [x] `pair_prf` zero-TP bug fixed and regression-tested.
- [x] Amendment 9 pinned, scorer audit table written, E0/T4/Turner control
      audited.
- [x] condition-a head seed17 retrained under corrected selection.
- [ ] condition-a head seed23 retraining (job 14792), then seed41.
- [ ] corrected P2-central a–e with the new heads; wait for PI before H2.4 or
      P2-block.

## Immediate local

- [x] Round 4j sent and archived; Claude ruled H2 retired and Round 4 P2
      withdrawn.
- [x] Corrected D1/K1/template/crosscheck artifacts pushed; the workspace
      `.git` remains read-only, so commits use the temp-object workaround.
- [ ] Keep probing with `scripts/remote.py doctor` plus an authenticated command;
      the public `:441` probe is not a liveness test.

## Next after retraining

- [ ] Dump native-head probabilities for P2-central a–e on id_test/withheld.
- [ ] Score v3 nested primary at fixed 0.5 and oracle thresholds; report
      kept–kept, all-pairs and pooled micro F1 for every condition.
- [ ] Enforce macro/micro ±.15 and shuffle-floor ≤.05 reporting invariants.
- [ ] After the PI sees corrected P2-central: P2-block + positive control, then
      H2.4 as a standalone decoder question.
- [ ] Restart `(s,b)` fit only with the new heads; freeze/hash before H2.4.

## Standing

- [ ] Keep `test_novel` embargoed: labels/counts only until H2/P5 are
      pre-registered and the PI authorises confirmatory scoring.
- [ ] RiNALMo-giga is AGPL-3.0; repository licence remains undecided.
- [ ] Baselines and relaxed-scan audit when there is slack.
