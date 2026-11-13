# Next Actions

Updated: 2026-11-13 23:40 (external three-tool cohort frozen; MXfold2 still training)

The governing documents are `experiments/001_e1prime/PREREGISTRATION.md`,
Amendments 1-14, and `state/amendments.lock.json`. Amendment 8 fixes
`labels_v3` as primary and secondary labels. Amendments 10/11 pre-register the
P2-block/P3/D1-core design. Amendment 12 freezes E-C1. Amendment 13 freezes
E-C2 and its closure conditions. Amendment 14 pins the one-shot `test_novel`
confirmatory family list.

## Closed this cycle

- [x] D1-core falsified; H3 not supported; H2 and H2.4 retired.
- [x] E-C1 supported (`rho=.338` coverage, `rho=-.242` pseudo-perplexity).
- [x] E-C1b failed (delta-rho `.0849 [.0093,.1597] < .15`).
- [x] E-C2 failed: fit at boundary, gate is a constant scale, id_test reverses.
- [x] `(s,b)` frozen at `s=.1875, b=0` (hash `be54ddab`...) before H2.4.
- [x] `labels_v3` 92,781/92,781 parse, sha `b8fe71ab`...; eval subset sha
      `9d1ab36c`...
- [x] Scorer `pair_prf` zero-TP bug fixed with regression tests; audit table
      in `docs/corrections.md`.
- [x] Frozen common cohort: withheld 16,059 seqs / 409 families; id_test
      1,287 / 155. Same ids for every method.
- [x] RNAfold, EternaFold, RNAstructure scored on that cohort; paired
      family-block bootstrap (10,000 resamples, seed 20261113).
- [x] Verdict: head 3-seed mean below all three tools on macro F1; no paired
      delta excludes zero. `docs/external_baselines_2026-11-13.md` section 5.
- [x] Regression test `tests/test_external_baseline.py` (common-cohort contract).

## Running now (server)

- [ ] `15631 rff_mx_train` - MXfold2 MixC retrain, 10 epochs. Still RUNNING at
      23:35 (6h31 elapsed, start 17:04). `ckpt/epoch-7` written 22:51, so
      ~50 min/epoch and 3 epochs left: training ETA ~01:20, plus export.
      Output `experiments/013_mxfold2/run/`, Slurm log `logs/mx_train_15631.out`.

## Immediate local

- [x] Script sync: `scripts/jobs/p2central_cohort_eval.sbatch` recovered from the
      server (the only job script that existed only there). The three tracked
      files the server shows as modified are its own stale checkout, not
      unrecovered work: current `main` already carries the `--max-length 1024`,
      `--tool rnastructure`, pairs-prediction and shortest-first changes.
- [x] MXfold2 inference path built and smoke-tested against the epoch-8
      checkpoint: `scripts/predict_mxfold2.py` (chunked, resumable, strict
      parse, embargo refusal) plus `scripts/jobs/mxfold2_fold.sbatch`.
      3/3 sequences folded; 16 sequences took 11.9 s on one CPU worker.
      `tests/test_predict_mxfold2.py` pins the parser and the embargo rule
      (178 passed, 1 skipped).
- [x] Paper A: new subsection 2.7 + Table T7 with the frozen single-cohort
      numbers, limitation 6 and claim 9 rewritten (no head advantage, not
      "blocked"); an independent numeric audit caught two misquoted ranges
      (seed17 delta lower bound, family win rate) and both are corrected.
- [x] Methods-vs-code audit of Paper A: 12 discrepancies corrected, including
      the head description (concatenation + residual conv stack, 32-d
      projection), decoder A complexity (O(n^3)/O(n^2)), the ViennaRNA
      pseudo-energy formula, the `<4`-pair diagnostic scope, and the amendment
      count. All recorded in `docs/corrections.md` (2026-11-14).
- [x] **Protocol deviation disclosed:** `fit_energy_calibration.py` assigned
      folds by `crc32(id) % 5`, while Amendment 7 required clan-grouped folds.
      The frozen `(s,b)=(.1875,0)` is at the grid boundary and E-C2 failed
      independently, so no conclusion depends on it; a clan-grouped refit is
      optional and would cost one 5-fold grid fit on train.
- [x] Committed the frozen cohort id lists and per-row exclusion files, so the
      section-5 table is reproducible without the prediction files. Verified
      md5 against the server copies.
- [ ] On 15631 completion: verify `model.pth`, `model.conf`,
      `train_manifest.json`; then `scripts/jobs/mxfold2_fold.sbatch` folds both
      splits (CPU, 16 workers).
- [ ] Rebuild the common cohort with `mxfold2_<split>.jsonl` added; expect the
      600 nt cap to shrink it. Rerun every score and
      `scripts/cohort_paired_eval.py`; mark section 5 of
      `docs/external_baselines_2026-11-13.md` superseded.
- [ ] Run the full local suite (>170 tests) and commit.
- [ ] Sync server-only job scripts (`external_baseline_rnastructure_idtest.sbatch`
      etc.) into the repo so runs are reproducible from GitHub.

## Next

- [ ] Round 4n to Claude: lead with E-C2 failure + constant-scale proof, the
      single-cohort table, the paired verdict, the RNAstructure IUPAC-vs-timeout
      correction, and H2.4 `pair_budget/not_confirmed`. Ask (a) when/whether to
      invoke the Amendment 14 one-shot `test_novel` unblinding and (b) Paper A
      title/abstract/discussion direction.
- [ ] Revise Paper A claim 9: a valid single-cohort comparison now exists but
      shows no head superiority; keep "no SOTA claim".
- [ ] H3 dose-response by homolog identity bin (exploratory, still unfiled).
- [ ] Enforce macro/micro +-0.15 and shuffle-floor <=0.05 invariants in every table.

## Standing

- [ ] Keep `test_novel` embargoed: labels/counts only until the PI authorises
      the Amendment 14 one-shot unblinding.
- [ ] One cohort only; never compare per-method subsets.
- [ ] RiNALMo-giga is AGPL-3.0; repository licence remains undecided.
