# Round 4h reply — Claude Opus 5.5 (High), PI

Source: https://claude.ai/chat/cc2ad40a-406e-481e-8e7b-89a64b78f2f4
Captured: 2026-11-12 (Asia/Shanghai), after the Round 4h report.

## Round 4h: the empty-label rows matter more than anything else here

The hash basis is clean, and the provenance audit has paid for itself twice: the
`.gitignore` hole and the nested-subset reference bug. Make the fresh-clone check
permanent with a CI job (GitHub Actions or equivalent) that, on every push,
clones the repo, installs it, runs the tests, and verifies the lockfile hashes
against blob bytes. Then “the repo is the record” is tested continuously rather
than claimed.

### 1. The 17.8% empty-truth seed rows: a label defect that may also explain part of Round 4

A sequence whose true structure was silently replaced by `pairs = []` because of
a parser failure has a corrupted label, not an unstructured RNA. Two
consequences:

- **Primary F1 can’t stay on v1.** Amendment 7 put primary F1 on v1, and v1 is
  defective for these rows.
- **Round 4’s withheld numbers are suspect in a direction that matters.**
  Withheld includes 11,527 seed rows. On a row with empty truth, predicting fewer
  pairs is never penalised and sometimes rewarded. Masking (condition e)
  plausibly suppresses predictions, so part of e’s +0.041 on withheld could come
  from these rows. I need to know how the scorer treats empty truth: does
  both-empty score 1, 0, or get excluded? Whatever the answer, K1 and D1 must
  also be run with empty-truth rows excluded, and both versions reported.

Why 8,870 rows still fail under v2: I’d bet the projection works on the string
level. If `SS_cons` is projected onto the ungapped sequence by deleting gap
columns, a bracket whose partner sat in a gap column becomes unbalanced. The
correct method works on pairs, not strings:

1. Parse `SS_cons` on alignment columns (all bracket types plus letter pairs).
   `SS_cons` is well-formed, so this always succeeds.
2. Map each column pair `(c1, c2)` to sequence coordinates, and drop the pair if
   either column is a gap in this sequence.
3. Apply the canonical/GU filter.

This should parse essentially 100% of rows. Build it as `labels_v3`:

- Report the parse rate.
- Report agreement with v2 on rows where v2 parses. It should be identical; any
  difference is a bug in one of them.
- Report how many rows still have zero or very few canonical pairs after correct
  projection. Those are real sparse labels, not parser failures.

**Amendment 8, committed before any primary scoring:**

- **Primary labels are v3, and secondary F1 also uses v3.** v1 and v2 are kept
  only for provenance. This is a correction made without outcome knowledge,
  since no primary F1 has ever been computed. Say so explicitly.
- **Scoring-time exclusion:** sequences with fewer than 4 true pairs under v3
  are excluded from primary F1. Their predicted-pair counts are reported
  separately as a false-positive diagnostic. Membership and hashes stay
  unchanged. This also replaces the “≥4 pairs in `SS_cons`” eval-subset
  condition with a per-sequence version, at scoring time.
- **Round 4 relabel:** all Round 4 numbers become “secondary-style, v1 labels
  including empty-truth rows”. They are never shown in the same table as primary
  results without that label.

Can run now, because the file is committed: count v1 empty-truth rows and rows
with fewer than 4 pairs in `test_novel_eval.jsonl`, and recompute both counts
under v3. Report the counts only, with no predictions involved; the embargo
isn’t touched. If more than about 10% of the eval subset drops out under the
fewer-than-4 rule, tell me before Amendment 8 is pinned.

### 2. The nested subset

Agreed on all points. Note that “pseudoknot components are tiny” was my
assumption, not yours; the 326-arc component settles it. Keep the 87 suboptimal
cases as a regression test.

### 3. Revised outage order

1. `labels_v3` with pair-level projection, the agreement checks against v2, and
   the eval-subset label counts (§1).
2. **Amendment 8** (§1).
3. **CI job** for fresh clone, tests and lock verification.
4. **Confirmatory-path rehearsal** on local ArchiveII plus the per-length timing
   table, as you planned.

### Cluster queue: one change

Step 1 (provenance under A6.1) now includes printing per-split counts of v1
empty-truth and v3 fewer-than-4-pair rows, using `--split-dir`. Everything after
it runs on v3 labels. That includes K1 and D1, which run on the existing
P2-central dumps with labels swapped at scoring time, in both versions
(empty-truth rows included and excluded).

Lead the next report with the v3 parse rate and the eval-subset counts.
