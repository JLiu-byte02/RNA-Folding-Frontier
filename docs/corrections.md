# Dated corrections

Append-only. Each entry names the artifact that was wrong, what it said, what is
true, and how the error was found. Frozen files are never edited; the correction
is recorded here and, where it affects data, in a separate versioned artifact.

## 2026-11-12 — Split manifest says "pseudoknots dropped"; they are not

- **Where:** `scripts/build_s1prime.py`, the generated `manifest.json` field
  `rules.pairs`: *"canonical AU/UA/GC/CG plus GU/UG; non-canonical and pseudoknots
  dropped"*.
- **True:** `filter_pairs` removes only non-canonical pairs. Crossing
  (pseudoknotted) pairs are **retained** in the labels. Measured on the bpRNA-1m
  mirror: 15,794 of 66,715 sequences (23.67%) contain at least one crossing pair,
  and 369,951 of 5,990,889 pairs cross (6.18%)
  (`experiments/005_label_policy/pseudoknot_prevalence.json`).
- **Found by:** reading the label pipeline for Amendment 6 §A6.2 against the
  manifest text.
- **Action:** the string is corrected in `build_s1prime.py` for any future
  rebuild. Existing manifests stay as they are — the split files are frozen and
  are not rebuilt for a wording fix.

## 2026-11-12 — ~7,670 seed rows enter the pool with an empty truth

- **Where:** `scripts/build_s1prime.py`: `pairs_from_dbn(dbn)` inside
  `try/except ValueError: pairs = []`.
- **What happens:** the v1 projection keeps only `()[]{}<>` and turns the other
  consensus characters into `.`. When the dropped level leaves a bracket
  unbalanced for that particular seed row, parsing raises and the row keeps its
  place in the split with **no pairs at all**. Any predictor is charged a false
  positive on every pair it proposes for such a row, and the row cannot
  contribute a true positive.
- **Measured:** 16,540 of 92,781 Rfam seed rows (17.8%) have a v1 projection that
  does not parse; the v2 WUSS projection parses 7,670 of those
  (`experiments/008_labels_v2/labels_v2_summary.json`). Per-split counts need the
  membership lists, which are server-side (Amendment 6 §A6.1).
- **Found by:** building `labels_v2` and counting parse failures explicitly
  instead of swallowing them.
- **Action:** not fixed in place — that would change frozen split hashes. The
  corrected pairs ship as the `labels_v2` delta artifact, and no confirmatory
  number may be quoted on seed rows without stating which label version was used.

## 2026-11-12 — WUSS letter pairs were dropped by the seed projection

- **Where:** `scripts/build_rfam_index.py:project_ss`.
- **True:** Rfam marks many pseudoknots with letter pairs (`AAA…aaa`); v1 maps
  them to `.`. Across all 92,781 seed rows this loses 59,632 pairs on 13,280 rows
  that otherwise parse, on top of the empty-truth rows above.
- **Found by:** Amendment 6 §A6.2 audit.
- **Action:** `rff.data.wuss` implements the v2 parser;
  `scripts/build_labels_v2.py` writes the delta artifact. Membership untouched.
