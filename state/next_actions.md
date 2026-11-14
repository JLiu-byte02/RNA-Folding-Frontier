# Next Actions

Updated: 2026-11-14 16:45 (four-method cohort frozen; Round 4n drafted, unsent)

Governing documents: `experiments/001_e1prime/PREREGISTRATION.md`, Amendments
1-14, and `state/amendments.lock.json`. Amendment 8 fixes `labels_v3` as primary
and secondary labels; 10/11 pre-register P2-block/P3/D1-core; 12 freezes E-C1;
13 freezes E-C2; 14 pins the one-shot `test_novel` list.

## Closed this cycle

- [x] Cluster SSH outage 00:33-14:23 diagnosed and recovered; a per-user banner
      watchdog (5 min) logs every probe. Detailed routing diagnosis is kept
      outside the public repository; the elevated split-tunnel route is **not**
      required while SSH works.
- [x] MXfold2 retrain `15631` completed 01:18 (8 h 13 m, exit 0); checkpoint,
      config and manifest hashes are recorded in `experiments/013_mxfold2/`.
- [x] Fold job `16459` (16,774 + 1,377 sequences) and cohort rebuild `16489`
      finished. Four-method cohort frozen: withheld 16,040 / 408 families,
      id_test 1,272 / 153; the 600 nt MXfold2 cap binds (597 / 572 nt).
- [x] Verdict unchanged: the head has no advantage over RNAfold, EternaFold or
      RNAstructure; no paired interval excludes zero. The retrained MXfold2 is
      weak (.5518 / .5491) and is reported as a retraining control only.
- [x] Length bands are descriptive; the within-family interaction is
      underpowered (24 withheld / 5 id_test families), so no length mechanism is
      claimed.
- [x] H3 homolog-identity dose-response filed as exploratory (R_c .876 vs .897;
      intervals overlap, not evidence).
- [x] Paper A tables and README are generated/checked against artifacts; CI now
      fails on number drift, on credentials, and on lock-hash disagreement.
- [x] Full suite: 184 passed, 1 skipped; fresh-clone CI green.

## Open

- [ ] **Send Round 4n** to the PI when a permitted browser channel exists. Draft:
      `docs/claude/transcripts/2026-11-14_round04n_prompt.md`.
- [ ] PI decisions: (a) whether/when to spend the Amendment 14 `test_novel`
      unblinding; (b) Paper A framing (audit-led vs negative-led vs coverage +
      decoder); (c) clan-grouped refit of `(s,b)`; (d) whether any method
      contribution is still defensible or the next cycle tests a new hypothesis.
- [ ] Optional compute decision: a longer MXfold2 schedule if it is to serve as
      a strong baseline rather than a control.

## Standing

- `test_novel` stays embargoed: labels and counts only until the PI authorises
  the one-shot unblinding.
- One cohort only; never compare per-method subsets.
- Invariants enforced by `scripts/jobs/rebuild_cohort_4tools.sbatch`: macro/micro
  gap <= .15 and shuffle floor <= .05.
- RiNALMo-giga is AGPL-3.0; the repository licence remains undecided.
- No credentials in the repository; `tests/test_public_repo_hygiene.py` enforces
  this in CI.
- If SSH drops again, read the local routing note before changing anything; the
  watchdog log makes the failure mode visible.
