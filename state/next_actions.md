# Next Actions

Updated: 2026-11-12 16:35 (after D1/K1/D2 completion)

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

- [x] D1 seeds 17/23/41, four label policies plus nested-truth primary.
- [x] K1 full-truth and nested-truth; neither fires.
- [x] D2 H2.1/H2.2 on seed17; H2.2 passes strongly, H2.1 confounded.
- [ ] D2 length-normalised H2.1 sensitivity (job submitted; SSH transiently
      reset, result not yet collected).

## Immediate local

- [ ] Commit/push the uncommitted D1/K1/D2 JSONs and
      `docs/d1_k1_d2_2026-11-12.md`. The workspace `.git` is read-only and the
      escalation approval review is timing out; retry when the approval channel
      recovers.
- [ ] Reopen Claude once ZoogVPN has a non-China IP. The current tab redirects
      to `claude.com/app-unavailable-in-region`; do not change VPN settings
      without user approval.
- [ ] Prepare Round 4i with the v3 parse rate, eval-subset counts, A6.1
      provenance, D1/K1 and D2 in that order.
- [ ] Record the PI's Round 4i ruling in a transcript and `docs/decision_log.md`.
- [ ] Keep probing with `scripts/remote.py doctor` plus an authenticated command;
      the public `:441` probe is not a liveness test.

## Next after diagnostics

- [ ] Recompute the `(s, b)` grid cost on the real train split; if below 24
      CPU-h, run the A7.4 coarse grid, otherwise use the pre-registered 3,000
      sequence seed-17 subsample. Freeze and hash `(s, b)` before scoring.
- [ ] P2-block cache only if the D1/K1 evidence still supports H2.

## Standing

- [ ] Keep `test_novel` embargoed: labels/counts only until H2/P5 are
      pre-registered and the PI authorises confirmatory scoring.
- [ ] RiNALMo-giga is AGPL-3.0; repository licence remains undecided.
- [ ] Baselines and relaxed-scan audit when there is slack.
