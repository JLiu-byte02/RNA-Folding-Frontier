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

## 2026-11-12 — `(s,b)` fit objective pooled counts, not the primary metric

- **Where:** `scripts/fit_energy_calibration.py:family_macro_f1`.
- **What was wrong:** it pooled tp/fp/fn per family and took the F1 of the
  pooled counts. The project's family-macro primary metric is the mean of
  per-family **mean per-sequence** F1; pooling weights long sequences more
  heavily and is a different selector.
- **Consequence:** no H2.4 number has been produced, but the first fit job
  started with the old objective and is cancelled; its output may not be
  frozen or quoted.
- **Action:** compute per-sequence F1 with the Amendment 9 zero-TP convention,
  then mean within family and across families; added
  `tests/test_energy_fit.py`. Refit from the existing p2fix train dump before
  freezing `(s,b)`.

## 2026-11-12 — `(s,b)` fit worker table was copied once per grid point

- **Where:** `scripts/fit_energy_calibration.py` `ProcessPoolExecutor`.
- **What was wrong:** the full item table, including a dense `L x L` score
  matrix per sequence, was passed through `initargs` for every one of the 25
  grid points. Slurm reported `MaxRSS` above 200 GB; the job produced no output
  and left orphaned joblib processes. The main process was gone while the
  allocation still showed `RUNNING`.
- **Action:** store the sparse scored map per item, inherit the read-only table
  through Linux `fork`, and build the dense matrix inside each worker. The
  Windows/local spawn fallback keeps the initializer path. The objective and
  grid are unchanged.
- **Evidence:** a synthetic 8-sequence fit completes end to end; no H2.4 result
  was produced or quoted from the stuck job.

## 2026-11-13 — Template-control artifact and P2-block b were not a global shuffle

- **Where:** `experiments/002_e1prime_runs/template_controls.json` and the
  provisional `score_p2_block.py` invariants.
- **What was wrong:** the JSON still contained the pre-fix `pair_prf` values
  (`shuffle_*` ≈ .89), even though the scorer had been corrected. Separately,
  P2-block condition b preserves the native block plus its 32-nt flank; it is a
  local-context swap, not a global sequence shuffle, so its F1 (≈ .56) cannot be
  used as a sequence-free floor.
- **Corrected template controls (3 seeds, fixed .5, v3 primary):** all-N
  `.000549` id_test / `.000066` withheld; dinucleotide shuffle `.016185` /
  `.017319`.
- **Action:** rerun `scripts/score_template_controls.py`; P2-block reports the
  b condition as `local_context_b_f1` and uses the corrected global shuffle
  mean for the ≤.05 invariant.

## 2026-11-13 - server reset overwrote two tracked-modified files

- **Where:** `junxuanliu@10.22.19.39:/home/junxuanliu/projects/RNA-Folding-Frontier`,
  files `docs/ec1_2026-11-13.md` and
  `experiments/002_e1prime_runs/ec1_p1_p2.json`.
- **What happened:** the server checkout was behind `origin/main` by several
  commits (it was at `4bc399e`). Before syncing, a backup loop was supposed to
  copy the two tracked-but-modified files aside. That loop failed on a
  PowerShell variable-interpolation error and did not report the failure, so
  the following `git reset --hard origin/main` replaced both files with the
  committed versions.
- **Why no reported result was lost:** the server JSON after the reset has
  sha256 `d3b8910f296b2751...9df12b`, 123,965 bytes, `n_families=409`,
  `supported=true`, `falsified=false`, `p1.rho=.337652 [.243549,.42936]`,
  `p2.rho=-.241971 [-.331827,-.146982]`; it is byte-identical to the local
  committed copy and reproduces every number quoted in `docs/ec1_2026-11-13.md`.
  Untracked server outputs were untouched by the reset.
- **Residual uncertainty:** the pre-reset server bytes were never hashed, so an
  edit that existed only on the server cannot be ruled out; the documented
  numbers are reproducible from the surviving artifact, not from that edit.
- **Action:** never run `reset --hard` on a checkout with local modifications.
  Sync procedure is now: `git fetch`, print `git status --porcelain`, and if it
  is non-empty, copy the listed paths to a timestamped directory and verify
  their hashes *before* any reset; the copy step must fail loudly (set `-e`
  or check the exit status) instead of continuing.
- **Evidence:** `git status --porcelain` on the server after the reset lists no
  tracked modifications; the hash above was re-checked with `sha256sum`.

## 2026-11-13 - RNAstructure's missing withheld rows are IUPAC skips, not timeouts

- **Where:** `docs/external_baselines_2026-11-13.md` sections 1, 4 and 4a, and a
  working note that attributed the 104 missing withheld predictions to the
  `<out>.errors.jsonl` timeout rows.
- **What was wrong:** two separate claims. (i) "111 ambiguous rows are inside
  our dump, so they will shrink the cohort by up to 111." (ii) The 104
  predictions missing from the cohort came from the timeout rows.
- **Verified** on the server by joining
  `rnastructure_test_withheld.{skipped,errors}.jsonl` to
  `dump_p2fix_seed17/a_test_withheld.jsonl`, the label-eligible set
  (`truth_v3` not null and >= 4 pairs) and `cohort/test_withheld.jsonl`:
  - `skipped.jsonl`: 111 rows, all `reason = non-ACGU`. 105 are in the dump;
    104 are label-eligible; 1 is label-ineligible; 6 exceed 1,024 nt and are
    not in the dump. The 104 label-eligible skips are exactly the cohort's
    `missing_prediction = 104`.
  - `errors.jsonl`: 41 rows, all `err = timeout`, lengths 3,337-4,381 nt
    (median 3,720). None is in the dump: all lie outside the 1,024 nt context,
    so none can belong to `missing_prediction`.
  - label-eligible rows without a prediction = 104 = 104 skipped + 0 errors +
    0 unexplained.
- **Action:** sections 1, 4 and 4a of the baseline document now state 111
  skipped / 104 inside the cohort / 41 timeouts outside the context. Any
  statement that timeouts explain the missing rows is withdrawn.
- **Evidence:** the counts above are reproducible from the two `.jsonl` files
  and the dump; the id sets are small enough to inspect directly.

## 2026-11-13 - `lt4 = 630` is `512 + 118`, not a single exclusion class

- **Where:** the coverage tables in `docs/external_baselines_2026-11-13.md` and
  working notes that quoted "630 rows below four pairs".
- **What was wrong:** `scripts/build_common_cohort.py` records
  `empty_truth` (118) and `lt_min_true_pairs` (512) as separate reasons;
  `scripts/score_external_baseline.py` reports their union as `lt4 = 630`.
  Calling 630 "rows with 1-3 pairs" overstates that class by the 118
  empty-truth rows, which are a label-parse defect, not short labels.
- **Action:** tables now print 118 + 512 = 630. No score changes; both
  populations were already excluded from the cohort by the same rule.

## 2026-11-14 - Energy calibration used record-level folds, not clan groups

- **Where:** `scripts/fit_energy_calibration.py:59-60,175-176` assigns folds by
  `zlib.crc32(record_id) % k`. `experiments/001_e1prime/AMENDMENT_7.md:111-112`
  (A7.1) specified "5-fold cross-validation grouped by clan so no clan is split
  across folds".
- **True:** the frozen `(s, b)` fit used record-level folds. Clan-mates of a
  held-out record can therefore sit in the training folds, so the fold F1 that
  selected `(s, b)` is optimistic relative to the pre-registered design.
- **Scope:** this affects the decoder `(s, b)` fit only. The later E-C2 gate
  fit (`scripts/fit_ec2_gate.py`) implements the clan-grouped design and is not
  affected; its failure therefore does not stem from this deviation.
- **Impact on `(s, b)`:** the selected point `(.1875, 0)` sits at the refined
  grid's upper boundary `s_max`, and the gate E-C2 built on top of it collapsed
  to a constant scale, bit-identical to its own control, with a sign reversal on
  id_test. No reported conclusion is rescued or overturned by the deviation,
  but the fold F1 behind `(s, b)` is optimistic relative to the pre-registered
  design. A clan-grouped refit has not been run; it would cost one 5-fold grid
  fit on the train split.
- **Found by:** a methods-vs-code audit of the Paper A draft on 2026-11-14.

## 2026-11-14 - Paper A Methods section did not match the code

- **Where:** `docs/paper_A_draft_2026-11-13.md` section 3.
- **What was wrong, and what is true:**
  - "outer-product ... symmetric outer-sum" described the head wrongly. The head
    concatenates `[z_i, z_j, z_i * z_j]`, adds a learned relative-distance
    embedding, runs eight dilated 3x3 residual blocks, and symmetrises the
    logits (`src/rff/models/pairhead.py:59-78`).
  - "BCE on v1 pairs ... and on v3 pairs thereafter" was wrong: every head was
    trained on the frozen bpRNA-era `pairs` field; v3 changes evaluation labels
    only (`scripts/train_head_s1.py:81-90`, Amendment 10).
  - "linear in memory and quadratic in time" for decoder A was wrong: it is
    O(n^3) time and O(n^2) memory (`src/rff/decoders/nn_dp.py:88`). The 33,714
    component / 87 suboptimal check belongs to the nested-truth solver
    `max_nested_subset_of_arcs`, not to decoder A.
  - "adds `s * (stacking energy) + b` to each cell's log-odds" was wrong:
    ViennaRNA receives a per-pair pseudo-energy `-(s * score_ij + b)`; Turner
    stacking stays inside ViennaRNA (`src/rff/decoders/vienna_soft.py:91`).
  - "`(s, b)` selected at the refined grid's interior" was wrong: it sits at the
    upper boundary `s_max = .1875`.
  - `s1prime` was named as the evaluated split and the cmscan flags were
    incomplete: the evaluated files are `s1prime_filtered` (train 17,746 /
    id_test 1,402 after the reverse-homology and 15-mer screens) and the
    assigner runs `--cut_ga --rfam --nohmmonly`.
  - "identical decoding" held only within model-arm comparisons; external tools
    use their native decoders.
  - The `<4`-pair predicted-pair diagnostic is emitted by the D1 and
    external-baseline scorers, not by `block_eval`, which drops those rows
    silently.
  - E-C1's adjustment set was understated: median length, GC, pair density,
    multiloop presence, and label source.
  - "Fourteen amendments (`AMENDMENT_*.md`)" counted files; there are 13 files
    (Amendments 2-14) plus Amendment 1 inside `PREREGISTRATION.md`.
  - `state/amendments.lock.json` does not pin label or split hashes; those live
    in `experiments/011_labels_v3/README.md` and
    `docs/remote_provenance_2026-11-12.md`. The embargo guard compares locked
    files against `HEAD`; it does not reject unrelated worktree changes.
  - The CI job installs test dependencies (the repository has no package
    metadata), and the repository is the record for code, protocols and
    evaluation artifacts, not for weights, raw logs, or corpora.
- **Action:** section 3 of the draft is rewritten to each of these statements.
  No result number changes; this is a description-vs-code repair.
- **Found by:** a read-only methods audit against the code and amendments.

## 2026-11-14 - Paper A Results mislabelled six populations

- **Where:** `docs/paper_A_draft_2026-11-13.md` section 2.
- **What was wrong, and what the artifacts say:**
  - 2.1 quoted ArchiveII as 3,370 records while the overlap and 313-duplicate
    statistics were computed on the 3,966-record corpus ingested for the leakage
    audit. Both populations are now named, and the exact/near counts are given
    as 2,427/3,966 and 3,506/3,966.
  - 2.2 called 92,781 "seed rows in our splits"; it is the full Rfam seed-label
    corpus, of which 54,400 rows appear in the frozen split files.
  - 2.2 said 9,322 of the 16,540 v1 failures still fail under v2. 8,870 is the
    residual; the other 452 are rows v1 parses and v2 does not. The table now
    shows 7,670 repaired / 8,870 residual / 452 new / 9,322 total.
  - 2.3 described the D1-core table as nested truth. The table is the all-pairs
    truth (`unpaired_allpairs`: .5959 / .5965 / .5860); the nested means are
    .5951 / .5956 / .5847 and are now quoted in the text.
  - 2.4 presented a five-condition ladder as one 16,163-sequence population.
    Condition c is restricted to the 8,476 sequences that have a same-family
    homolog block; the text now separates the four full-population conditions
    from c.
  - 2.5's "favourable below 300 nt" held only on withheld. On id_test the
    150-299 nt bin adds 125 long-range FPs and the <80 nt bin loses 41 TPs.
  - Claim 3's "8.1% of the rest" was ambiguous; 8.1% is 6,757 / 83,459, i.e. of
    rows where both v2 and v3 parse.
  - 2.5 omitted the decoder-A constant-scale comparison, the only arm that beats
    a thermodynamic tool: .6137 vs RNAfold .6003 on withheld, paired
    +.0134 [+.0025, +.0241], but +.0204 [-.0020, +.0392] on id_test and below
    the pre-registered best-baseline bar (`beats_baseline_by_min_gain: false`).
    It is now reported, and 2.7's "no advantage" is narrowed to the default
    greedy decoder.
- **Action:** section 2 and claim 3 rewritten; no metric changes. The three
  E-C2 score JSONs are now committed so the deltas above are reproducible.
- **Found by:** a read-only results audit of each number against its artifact.

## 2026-11-14 - MXfold2 export manifest mislabels its label source

- **Where:** `experiments/013_mxfold2/export_manifest.json`, field
  `label_policy`: "vim3 pairs projected to maximum nested subset".
- **True:** `scripts/export_mxfold2_bpseq.py` reads the frozen `pairs` field of
  `splits/s1prime_filtered/train.jsonl`, which is the bpRNA-era projection, not
  v3. Only the string is wrong; the exported structures are the ones intended
  (the same labels every head was trained on). The manifest is regenerated by
  the export and its text is hard-coded in the script, so the script is fixed
  here and the existing manifest is left untouched.
- **Scale of the export:** 17,746 train rows -> 17,155 exported (577 longer than
  600 nt, 14 with fewer than one pair); 1,846 rows had crossing pairs, projected
  to the maximum nested subset (454,222 -> 439,299 pairs). Train-split sha256 is
  `b8af2335...f94ecd`, list sha256 `d2ba5c96...b5291c`.
- **Found by:** reading the export manifest while documenting the baseline.
