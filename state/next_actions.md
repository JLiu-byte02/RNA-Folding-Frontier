# Next Actions

Updated: 2026-11-12 17:25 (after Claude Round 4i ruling)

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
- [x] D2 length-normalised H2.1 sensitivity: id_test 1.62× per candidate cell,
      withheld 3.89×.
- [x] D2C FP-helix energetics: FP helices are energetically close to true
      helices; 26.6%/26.7% fully in RNAfold MFE.
- [x] v3 cross-check: 37,797 exact-sequence comparisons; v3 F1 .915 vs v2 .823.
- [ ] all-N and dinucleotide-shuffle native-head controls: dumps complete for
      all 3 seeds; fast fixed-threshold scoring job 14710 is queued. The slower
      25-point D1 scoring path was replaced by `score_template_controls.py`.
- [ ] `(s,b)` coarse fit on the full train split (4.3 CPU-h estimate; job 14662).
- [ ] H2.4 A-energy fixed on id_test/withheld (depends on fit).

## Immediate local

- [x] Round 4i sent and archived; ZoogVPN is connected and Claude is reachable.
- [x] D1/K1/D2 results, report and bundle pushed to GitHub (`e713110`, `51dedae`).
- [ ] Commit the new Round 4i control scripts and crosscheck/bundle updates via
      the temp-object workaround; the workspace `.git` itself remains read-only.
- [ ] Keep probing with `scripts/remote.py doctor` plus an authenticated command;
      the public `:441` probe is not a liveness test.

## Next after diagnostics

- [ ] Freeze and hash `(s,b)` before H2.4 scoring; then run H2.4 on id_test and
      withheld with the fitted values.
- [ ] If all-N/shuffle F1 is high, add all-N as condition f in P2-block.
- [ ] P2-block cache/H2.3 (~20 GPU-h) only after the template controls are read.
      The per-split cache template is ready at
      `scripts/jobs/p2_block_split.sbatch`; do not submit yet.
      The dependent H2.3 template is ready at
      `scripts/jobs/p2_block_eval_split.sbatch`.

## Standing

- [ ] Keep `test_novel` embargoed: labels/counts only until H2/P5 are
      pre-registered and the PI authorises confirmatory scoring.
- [ ] RiNALMo-giga is AGPL-3.0; repository licence remains undecided.
- [ ] Baselines and relaxed-scan audit when there is slack.
