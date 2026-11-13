# Next Actions

Updated: 2026-11-13 24:00 (external three-tool cohort frozen; MXfold2 still training)

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

- [ ] `15631 rff_mx_train` - MXfold2 MixC retrain, 10 epochs; epoch-7 present at
      23:40, ~48 min/epoch, ETA ~02:00. Output `experiments/013_mxfold2/run/`.

## Immediate local

- [ ] On 15631 completion: verify `model.pth`, `model.conf`,
      `train_manifest.json`; `scripts/export_mxfold2_bpseq.py` folds both splits.
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
