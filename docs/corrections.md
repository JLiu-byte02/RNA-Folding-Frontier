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

## 2026-11-12 — `pair_prf` scored zero-overlap rows as F1=1.0

- **Where:** `src/rff/analysis/labels.py:pair_prf`.
- **What was wrong:** when `precision + recall == 0` the function returned
  `f1 = 1.0`. That is correct only for the fully vacuous case where truth and
  prediction are both empty (`tp + fp + fn == 0`). For non-empty truth and
  prediction with zero overlap it must return `0.0`.
- **Consequence:** every sequence with no true positives scored as perfect. The
  error inflated family-macro F1 and affected, at minimum:
  - D1 and therefore K1 (per-sequence F1);
  - the all-N/shuffle template controls;
  - `train_head_s1.py` / `train_pairhead.py` validation F1 used for model
    selection;
  - earlier label-discordance and v2/v3 cross-check F1 summaries.
- **Found by:** the shuffle control produced family-macro F1 ≈0.89 while its
  pooled micro F1 was ≈0.009, which is mathematically impossible if both are
  computed from the same pair counts. A row-level print showed zero TP with F1
  1.0.
- **Immediate corrected evidence (template controls):** all-N family-macro F1
  is ≈0.00055 (id_test) and ≈0.000066 (withheld); the dinucleotide shuffle is
  ≈0.0162 and ≈0.0173. The earlier 0.89/0.90 was entirely the bug.
- **Immediate corrected D1 evidence (seed 17, v3 primary):** id_test nested
  a≈0.58, e≈0.28; withheld v1-all a≈0.55, e≈0.28. The earlier claim that e
  improves F1 is reversed; masking the distal region hurts.
- **Action:** fix the edge case, add `test_pair_prf_zero_overlap_is_not_vacuous`,
  rerun D1/K1/template controls and the v3 cross-check, and mark every earlier
  number based on `pair_prf` as invalid until rerun. Frozen amendments are not
  edited.

## 2026-11-12 — Scorer-scope audit table

This table is the Round 4j gate. A number not listed here as unaffected, or
recomputed below, may not appear in a report.

| artifact / number | scorer | old value | corrected / status |
|---|---|---|---|
| E0 ArchiveII vs bpRNA TR0 mean pair-F1 | `pair_prf` | .9191 | **.9195** on the server's 3,370-record ArchiveII mirror; the old value used a 3,966-record mirror. Telomerase .6606 unchanged. Re-run on the frozen E0 archive before a paper number. |
| E0 ArchiveII vs RNAStrAlign mean pair-F1 | `pair_prf` | .9946 | **.9937** on the 3,370-record mirror; identity 91.2% → 90.1%. Input-version caveat as above. |
| T4 baseline reproduction, ArchiveII RNAfold | `prf` (not `pair_prf`) | .5198 family-macro / .5407 micro | **unaffected** |
| Turner-control / independent-scorer crosscheck | `prf` + independent `f1_from_counts` | .5198 / .5407 | **unaffected**; the independent scorer did not use `pair_prf` and did not share the bug. |
| Pipeline smoke test | `prf` | reported smoke values | **unaffected** (pipeline-only, not a protocol number) |
| Round 4 P2-central a–e | `pair_prf` | old table | **withdrawn**; corrected table pending new heads |
| K1 history | `pair_prf` via D1 per-sequence | old e−a positive | **corrected:** e−a ≈ −.31 id_test, ≈ −.30 withheld; CIs entirely below zero |
| all-N / dinucleotide-shuffle template controls | `pair_prf` | .89 / .90 shuffle | **corrected:** all-N .00055/.000066; shuffle .0162/.0173 |
| v2/v3 cross-check against bpRNA-1m | `pair_prf` | v2 .8226 / v3 .9152 | **corrected:** v2 .8148 / v3 .9056 |
| D2 H2.1/H2.2 | `prf` / count arithmetic | H2.2 ≈4×, H2.1 split-dependent | **unaffected** |
| D2C FP-helix energetics | `prf` counts + ViennaRNA | 26.6%/26.7% in MFE | **unaffected** |
| Region-stratified D1 | `prf` (not `pair_prf`) | one/neither ≈0 | **unaffected** |
| labels_v3 parse/agreement counts | no pair scorer | parse 100%, 6,757 v2/v3 differences | **unaffected** |

The P2-central table must be recomputed with the retrained heads (Amendment 9
§A9.4–A9.5) before it can be quoted.

## 2026-11-12 — P2-central retention percentages used mismatched populations

- **Where:** the first paragraph of “Immediate reading” in
  `docs/p2central_fixed_2026-11-12.md`.
- **What was wrong:** condition c exists only on records with a verified
  homolog (449 id_test, 7,572 withheld after the ≥4-pair rule), while the a and
  e point estimates in that table came from the full split. `(c−e)/(a−e)` then
  compared a numerator from the c subset with a denominator from a different
  population. The 77% / 102% kept-window figures are not a like-for-like
  retention estimate.
- **Corrected, common a/c/e population (family-bootstrap 95% CI):**
  - id_test kept–kept: `R_c = .612 [.413, .853]`, `R_d = −.173 [−.617, .080]`;
  - withheld kept–kept: `R_c = .747 [.564, .960]`, `R_d = .082 [−.040, .203]`;
  - withheld all-pairs: `R_c = .504 [.390, .606]`, `R_d = −.069 [−.103, −.036]`.
- **Consequence:** on P2-central, c retention is below the H3 90% bar on both
  splits, with wide intervals. P2-central remains exploratory; H3 is decided on
  P2-block, where every condition is measured on the same population.
- **Artifact:** `experiments/002_e1prime_runs/p2central_ratios_ci.json` and
  `docs/p2central_ratios_ci_2026-11-12.md`.
