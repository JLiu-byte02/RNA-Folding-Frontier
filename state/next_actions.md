# Next Actions

Updated: 2026-11-12 (after PI Round 4h and Amendment 8)

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

## Immediate local

- [ ] Send Round 4i to the PI, leading with the v3 parse rate and eval-subset
      counts; attach commit `a9de0a2` and the CI/lock evidence.
- [ ] Record the PI's Round 4i ruling in a new transcript + `docs/decision_log.md`.
- [ ] Keep probing the tunnel with `scripts/tunnel_watch.py --once`; do not
      change network configuration or stop aTrust.

## When the cluster returns

- [ ] A6.1 provenance: verify the remote worktree, membership-list hashes and
      `git status`; push any missing frozen split metadata.
- [ ] Run `scripts/build_labels_v3.py --split-dir splits/s1prime` (or the
      actual split directory) and record per-split v1-empty / v3-`<4` counts.
- [ ] D1/D2 on the existing P2-central dumps with v3 labels swapped at scoring
      time, in both versions: empty-truth rows included and excluded.
- [ ] K1 with oracle threshold + one-partner decoding on the same dumps, both
      label versions; if K1 fires it goes on the first line of the report.
- [ ] A6.1 split provenance must precede any quoted split number.
- [ ] Recompute the `(s, b)` grid cost on the real train split; if below 24
      CPU-h, run the A7.4 coarse grid, otherwise use the pre-registered 3,000
      sequence seed-17 subsample. Freeze and hash `(s, b)` before scoring.

## Standing

- [ ] Keep `test_novel` embargoed: labels/counts only until H2/P5 are
      pre-registered and the PI authorises confirmatory scoring.
- [ ] RiNALMo-giga is AGPL-3.0; repository licence remains undecided.
- [ ] Baselines and relaxed-scan audit when there is slack.
