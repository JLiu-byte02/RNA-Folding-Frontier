# Amendment 4 — H2.4 comparison order, the two A-energy variants, and the reproducibility gap

Status: written and committed **before any D1/D2 dump exists and before any
H2.4 number has been computed**. It records the PI's Round 4d rulings
(`docs/claude/transcripts/2026-11-12_round04d_reply.md`). Amendment 3 is not
superseded; this file fixes the details Amendment 3 left open.

## A4.1 — Comparison order for H2.4 (frozen)

Decoder A is compared with **greedy nested matching** (threshold, then
one-partner, then non-crossing — the production decoder of Round 4). The
comparison is fixed in this order:

1. **Primary: recall-matched.** One greedy threshold per split and seed, a single
   global value, chosen so that pooled recall equals decoder A's; if no threshold
   matches exactly, take the closest point *at or above* A's recall. Choosing the
   threshold on the evaluation data favours greedy, so the test is conservative
   toward the null.
   **Verdict:** H2.4 is supported if A has fewer long-range false positives
   (`|i-j| > 128`, pooled) **and** at least as many true positives.
2. **Fallback, only when no threshold reaches A's recall: pair-budget matched.**
   Choose the greedy threshold whose total predicted pair count equals A's.
   **Verdict:** supported only if A has **both** more true positives and fewer
   long-range false positives at equal count. This is stricter than the recall
   ceiling alone.
3. **Always reported, descriptive only: greedy's recall ceiling.** Maximum pooled
   recall over all thresholds.
4. **If neither match is feasible:** report `inconclusive` for that split, length
   bin or seed, and report how often that happens. No further fallback may be
   added after seeing data.

All three readouts are reported per length bin, but the **verdict is the pooled
result**.

## A4.2 — A-energy-fixed: the confirmatory variant for H2.4

`A-energy-fixed` is what Amendment 3 calls the confirmatory A-energy variant.

- Pseudo-energy per candidate pair: `dG_ij = -(s * logit_ij + b)`.
- The offset `b` is required because `pos_weight = 20` shifts the logits by
  about `log 20`; without it the scale `s` would absorb calibration and the
  decoder would no longer be a monotone function of the ranking greedy sees.
- `(s, b)` are fitted by **5-fold cross-validation on the frozen train split**,
  maximising family-macro F1, then frozen. The fitted values are hashed into this
  amendment's companion record before id_test or withheld is scored again.
- Because `(s, b)` are fixed and the map is monotone, A-energy-fixed isolates
  "energy-consistent decoding vs greedy" on the same pair scores. That is the
  H2.4 claim.
- Exploratory, reported but not decisive: F1 as a function of `s` on id_test, so
  a reviewer can see the verdict does not hinge on one lucky value.

## A4.3 — A-energy-learned: the P5 variant, not evidence for H2.4

- Head trained with structured max-margin **through** the decoder, with the
  per-pair bonus `tanh`-bounded at a cap `C` kcal/mol.
- `C` is chosen by cross-validation on the frozen train split over
  `{1, 2, 3, 5}` and then frozen.
- Its advantage over A-energy-fixed is end-to-end calibration of the energy
  scale. It is the method entered against BPfold in P5, and it is reported as
  such. It is **not** used to adjudicate H2.4.

## A4.4 — Narrowed novelty framing for H2.4

The claim must be stated as a **controlled decoder comparison on fixed FM pair
scores**, not as "no one has combined learned scores with an exact energy DP".
That stronger statement is false for the field: MXfold2 combines learned scores
with Turner terms inside an exact DP, and BPfold places energy inside its model.
Both co-train the score and the decoder, so neither isolates the decoder.

The defensible question: *does energy-consistent decoding remove the
length-scaling false positives Flamm et al. identified, compared with greedy
nested matching on the same scores?* Flamm et al. leave that question open.

## A4.5 — Reproducibility gap, recorded rather than worked around

The frozen split files (`splits/s1prime_filtered/*.jsonl`) and their manifests
(`splits/s1prime/manifest.json`, `splits/s1prime_filtered/filter_report.json`)
exist **only on the cluster**, and no membership list or split hash is committed
to this repository. Until the cluster returns, no local run can verify that a
rebuilt split is the frozen one. This is a reproducibility gap and it is
recorded as one:

- Local rebuild requires `cmscan` (Infernal) against Rfam 15.0 plus the pinned
  bpRNA-1m mirror; the rebuild is only meaningful if its sha256 matches the
  frozen manifest, which is currently unreachable.
- Consequence: CPU baselines run locally before the cluster returns are
  **pipeline exercises on non-protocol data**, not protocol numbers, and must be
  labelled that way.
- First actions when the cluster returns: copy both manifests into the
  repository, then re-check every split hash.

## A4.6 — Provenance resolution procedure (added after PI Round 4e)

Added on request of the PI (Round 4e) to fix exactly how the A4.5 gap is closed.
This section was appended after the hash in `docs/decision_log.md` was recorded,
so that recorded value no longer matches this file; A4.1–A4.5 are unchanged and
the new hash is recorded alongside the old one. Nothing here changes a claim.

Order of operations when the cluster returns, before any new computation:

1. `git status` and `git log origin/main..HEAD` on the server, then push
   everything, including per-split sequence-ID membership lists and per-split
   sha256 hashes.
2. Record the modification times of the split files and of the test_novel freeze
   file, plus a statement that no scoring output exists for test_novel.
3. Determinism check: rebuild every split from committed code, seeds and the
   pinned inputs into a fresh directory and compare hashes with the on-disk
   originals. Match → the freeze is credible after the fact. Mismatch → stop and
   report before going further.

Partial local fact, recorded 2026-11-12: the *evaluation subset*
`experiments/001_e1prime/test_novel_eval.jsonl` is committed, and its sha256
recomputes to `9d1ab36ca8d07a47b67b5972bf69b2be6ebd34127b777eabbc552b0706ee209b`
over 2,372 records, exactly as declared in `novel_eval_set.md`. That is an
artifact for the subset, not for the 42,873-record `test_novel` split, whose
membership list and hash remain uncommitted.
