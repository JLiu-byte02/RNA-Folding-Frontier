# Next Actions

Updated: 2026-11-12 21:20 (after Claude Round 4k ruling)

The governing documents are `experiments/001_e1prime/PREREGISTRATION.md`,
Amendments 1–10, and `state/amendments.lock.json`. Amendment 8 supersedes the
label clause of Amendment 7: Rfam seed primary **and** secondary labels are
`labels_v3`; v1/v2 are provenance only. Amendment 10 pre-registers P2-block,
P3, H3, and the D1-core rule.

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
- [x] All three corrected condition-a heads retrained (`runs/p2fix_seed*`).
- [x] Corrected P2-central a–e table completed; invariants pass.
- [x] Amendment 10 written and hashed `5ec10662…` before any P2-block cache
      was read.

## Running now

- [x] `pair_prf` zero-TP bug fixed and regression-tested.
- [x] Amendment 9 pinned, scorer audit table written, E0/T4/Turner control
      audited.
- [ ] Build P2-block test caches for b/c/d/e on id_test/withheld, then score
      a–e at fixed 0.5 and compute H3.
- [ ] Build train caches for b-block and e-block, train P3 heads (3 seeds each),
      and compute the D1-core verdict on withheld.
- [ ] CPU in parallel: dump p2fix_seed17 train logits, fit/freeze `(s,b)`,
      record D2C FP-helix expectation, then H2.4.

## Immediate local

- [x] Round 4j sent and archived; Claude ruled H2 retired and Round 4 P2
      withdrawn.
- [x] Corrected D1/K1/template/crosscheck artifacts pushed; the workspace
      `.git` remains read-only, so commits use the temp-object workaround.
- [x] Round 4k reply archived; P2-block and P3 authorized.
- [ ] Keep probing with `scripts/remote.py doctor` plus an authenticated command;
      the public `:441` probe is not a liveness test.

## Next

- [ ] P2-block test-time table a–e on id_test/withheld, with H3 ratios and CIs.
- [ ] P3 train-and-test table a/e-block/b-block, with the D1-core verdict.
- [ ] Bootstrap CIs for the corrected P2-central retention ratios.
- [ ] H3 dose-response by homolog identity bin.
- [ ] Enforce macro/micro ±.15 and shuffle-floor ≤.05 reporting invariants.
- [ ] H2.4 only after the `(s,b)` freeze is committed.

## Standing

- [ ] Keep `test_novel` embargoed: labels/counts only until H2/P5 are
      pre-registered and the PI authorises confirmatory scoring.
- [ ] RiNALMo-giga is AGPL-3.0; repository licence remains undecided.
- [ ] Baselines and relaxed-scan audit when there is slack.
