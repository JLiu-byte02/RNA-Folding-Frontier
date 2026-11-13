# Next Actions

Updated: 2026-11-13 23:10 (E-C2 failed; external-baseline cohort in progress)

The governing documents are `experiments/001_e1prime/PREREGISTRATION.md`,
Amendments 1–14, and `state/amendments.lock.json`. Amendment 8 fixes
`labels_v3` as primary and secondary labels. Amendments 10/11 pre-register the
P2-block/P3/D1-core design. Amendment 12 freezes E-C1. Amendment 13 freezes
E-C2 and its closure conditions. Amendment 14 pins the one-shot `test_novel`
confirmatory family list.

## Closed this cycle

- [x] D1-core falsified; H3 not supported; H2 and H2.4 retired.
- [x] E-C1 supported (`rho=.338` coverage, `rho=-.242` pseudo-perplexity).
- [x] E-C1b failed (Δρ `.0849 [.0093,.1597] < .15`).
- [x] E-C2 failed: fit at boundary, gate is a constant scale, id_test reverses.
- [x] `(s,b)` frozen at `s=.1875, b=0` (hash `be54ddab…`) before H2.4.
- [x] `labels_v3` 92,781/92,781 parse, sha `b8fe71ab…`; eval subset sha
      `9d1ab36c…`.
- [x] Scorer `pair_prf` zero-TP bug fixed with regression tests; audit table
      in `docs/corrections.md`.
- [x] Per-tool external baselines scored (RNAfold, EternaFold) on id_test and
      withheld; policy in `docs/external_baselines_2026-11-13.md`.

## Running now (server)

- [ ] `15631 rff_mx_train` — MXfold2 MixC retrain, 10 epochs, ~48 min/epoch,
      ETA ~03:00. Output `experiments/013_mxfold2/run/`.
- [ ] `15649 rff_rs_base` — RNAstructure 6.6 withheld, will hit the 8h wall at
      ~9,700/17,060; **resume needed from `rnastructure_test_withheld.partial.jsonl`**.
- [ ] `15820 rff_rs_idtest` — RNAstructure id_test, folds then scores.

## Immediate local

- [ ] Sync the server-only `scripts/jobs/external_baseline_rnastructure_idtest.sbatch`
      into the repo so the run is reproducible from GitHub.
- [ ] Re-run the full local suite (>170 tests) and commit.
- [ ] Draft Paper A methods + audit sections; hold title/abstract/discussion.

## Next

- [ ] When 15649 ends: resubmit resume+score with the same `--cache`,
      `--workers 16`, withheld scoring block only.
- [ ] When 15631 ends: verify `model.pth`/`train_manifest.json`, export with
      `scripts/export_mxfold2_bpseq.py`, fold both splits, score.
- [ ] Build the one frozen cohort over four tools + `dump_p2fix_seed17`; assert
      `n_missing_prediction == 0`; rescore all with `--restrict-ids`; append to
      `docs/external_baselines_2026-11-13.md` §5.
- [ ] Round 4n to Claude: lead with E-C2 failure, the constant-scale proof, the
      O(n²) incident, and the cohort status. Ask explicitly when/whether to
      invoke the Amendment 14 one-shot `test_novel` unblinding.
- [ ] H3 dose-response by homolog identity bin (exploratory, still unfiled).
- [ ] Enforce macro/micro ±.15 and shuffle-floor ≤.05 invariants in every table.

## Standing

- [ ] Keep `test_novel` embargoed: labels/counts only until the PI authorises
      the Amendment 14 one-shot unblinding.
- [ ] One cohort only; never compare per-method subsets.
- [ ] RiNALMo-giga is AGPL-3.0; repository licence remains undecided.
