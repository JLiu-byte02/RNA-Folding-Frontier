# Round 4j reply — Claude Opus 5.5 (High), PI

Source: https://claude.ai/chat/cc2ad40a-406e-481e-8e7b-89a64b78f2f4
Captured: 2026-11-12, after the `pair_prf` correction.

## Round 4j: the bug voids Round 4’s evidence rather than settling it. Don’t pivot yet

Finding this bug through the gap between macro and micro F1 was the right instinct.
I should have asked for that consistency check myself before building H2 on the
Round 4 table. The spec I gave in Round 4f listed edge cases but never stated the
zero-true-positive convention, so part of this is mine.

### What the bug actually invalidates

Every per-sequence family-macro number computed with `pair_prf` is void. That
covers more than D1/K1:

- **The whole Round 4 P2-central table.** All five conditions (a–e), the “P2
  inverts the prediction” verdict, (c) being worst, (d) ≥ (a), and the
  length-sign table if its “micro” column was a mean of per-sequence F1. The
  original family-shortcut test has not been answered; its Round 4 verdict is
  withdrawn, not reversed.
- **H2.** It was inferred from that table, and its premise (“e helps”) is false.
  As a mechanism, H2 is retired.
- **Checkpoint selection.** The bug gave F1 = 1.0 to sequences with zero true
  positives, so selection rewarded models that fail completely on hard
  sequences. The current checkpoints shouldn’t be reused for anything.
- **Possibly the E0 audit numbers.** Check whether the label-discordance F1s
  (mean 0.919, per-family values, telomerase 0.66) used `pair_prf`. Zero-overlap
  duplicate pairs would have scored 1.0 and inflated those means.
- **The independent-scorer cross-check.** It agreed with ours exactly on 3,370
  RNAfold predictions, and ArchiveII almost certainly includes sequences where
  RNAfold finds no true pairs. If both scorers gave those 1.0, the two
  implementations shared the convention error.

What stands: D2’s count-based results (H2.2 about 4×, H2.1 split-dependent), the
region-stratified counts, the corrected v3 cross-check (0.906 vs 0.815), the
controls showing no template prior, and all label/provenance work.

### Rulings

1. **Cancel the `(s, b)` fit** unless its log shows it started from a commit
   with the fix and uses the new checkpoints. Neither is true. Restart it after
   step 3.
2. **Scorer convention, frozen in Amendment 9:**
   - both truth and prediction empty: F1 = 1;
   - zero true positives with either set non-empty: F1 = 0;
   - primary F1 excludes sequences with fewer than 4 true pairs, so the
     empty-truth case never reaches the primary table.
   Add a test showing that the old convention fails it.
3. **Audit everything the scorer touched, before any new GPU work.** Add one
   table to `docs/corrections.md`: every artifact or number ever reported, which
   scoring function produced it, the old value, and the corrected value. This
   covers E0, T4, the Turner control, the smoke test, P2-central, the K1
   history and the v3 cross-check. Every number in the reports is either in that
   table as unaffected (with the reason) or recomputed.
4. **Retrain the native heads** (condition a, 3 seeds) with corrected
   selection. At about 5 minutes per head, retraining is cheaper than
   re-evaluating old checkpoints.
5. **Re-run P2-central for all five conditions with the new heads.** Regenerate
   the b, c and d probability dumps from the caches you already have. Score with
   v3 nested labels at both a fixed 0.5 threshold and the oracle threshold,
   three ways:
   - kept–kept pairs only (the valid readout);
   - all pairs;
   - pooled micro F1.
   The e result, much worse, is a sanity check rather than a finding. The
   informative contrasts are b, c and d on kept–kept pairs.
6. **Standing invariants, enforced automatically in every report:**
   - family-macro and pooled-micro F1 printed side by side; report fails if they
     differ by more than 0.15 without a stated reason;
   - a shuffle-control floor on every split; if it scores above 0.05, the report
     fails.
   These are reporting checks, not new pre-registration.

### What we’re testing now

We’re back to the original D1 question, with a corrected instrument:

- **P2-block** (the pre-registered geometry) plus the **positive control**, after
  step 5 shows the instrument behaves sensibly. Roughly 20–25 GPU-h.
- **H2.4 stays as a standalone question:** on fixed FM pair scores, does
  energy-consistent decoding cut long-range false positives versus greedy
  nested matching? It never depended on “e helps.” H2.2 is the descriptive
  observation that motivates it. The helix-energetics check still applies.
- **H2.3 folds into P2-block.**

Decision point after corrected P2-central, P2-block and the positive control:

- shortcut evidence appears (c retains performance vs d in the positive control,
  and S1′ contrasts show it): D1 is alive;
- positive control detects identity use but S1′ shows none: a well-controlled
  negative result; the paper rests on audit + H2.4 + P5; D2 becomes the method
  attempt;
- instrument fails the positive control: redesign the probe before drawing any
  conclusion.

The pivot decision waits for those results, not for numbers from a voided table.

### Records

`docs/decisions.md` gets three entries:

1. the Round 4 P2 verdict is withdrawn because of the scorer bug;
2. H2 is retired as a mechanism because its premise is false;
3. H2.4 continues as a standalone question.

Amendment 9 covers the scorer convention, the head retraining, the recomputation
plan and the new invariants, pinned before step 4. The embargo is unchanged.

### Next report

1. The scorer-scope table, including whether E0’s label-discordance numbers move.
2. Whether the independent scorer shared the convention, and its corrected
   agreement.
3. The P2-central table from the new heads: a–e, kept–kept / all-pairs / micro,
   at both thresholds.
4. The Amendment 9 hash.

H2.4 and P2-block are not to be launched until I’ve seen item 3.
