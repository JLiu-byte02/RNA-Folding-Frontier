# Next Actions

Updated: 2026-11-15 14:25 (Round 4n send attempted; PI reply pending network recovery)

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
- [x] FP anatomy lead corrected: the first pass used thresholded probability
      maps, not the deployed greedy head. Corrected artifacts now match the
      frozen head-17 counts exactly; the radius-2 shift oracle bounds the
      possible macro gain at +.0251/+.0263, not a reachable result.
- [x] Paper A tables and README are generated/checked against artifacts; CI now
      fails on number drift, on credentials, and on lock-hash disagreement.
- [x] Full suite after the correction: 187 passed, 2 skipped.
- [x] SSH doctor at 13:52 +0800: TCP reachable, 21-byte OpenSSH 7.4 banner;
      no retry storm.
- [x] GitHub push succeeded through a direct-IP one-shot SSH route; fresh-clone
      CI run 125 completed successfully at 06:02 UTC.

- [x] SSH recovered again at 12:43 (watchdog banners). `family_diagnostics.py`
      reproduces the three committed diagnostics on the server with identical
      MD5s; runtime scripts uploaded individually. No `reset --hard` on the
      server worktree (untracked dumps at risk).

## Open

- [ ] **Recover the Round 4n exchange.** The message was sent into the existing
      Claude chat at 14:13 +0800, but the generation stream died and the page
      then failed DNS. Reopen the chat after ZoogVPN data-plane recovery, check
      whether the message is in history, and do not duplicate it blindly. The
      exact body is `docs/claude/transcripts/2026-11-14_round04n_prompt.md`.
- [ ] **Human/VPN action:** reconnect or switch the ZoogVPN node from its UI with
      administrator rights. System DNS and Claude access are down while SSH to
      `10.22.19.39` remains up; see `docs/network_2026-11-15.md`.
- [ ] PI decisions: (a) whether/when to spend the Amendment 14 `test_novel`
      unblinding; (b) Paper A framing (audit-led vs negative-led vs coverage +
      decoder); (c) clan-grouped refit of `(s,b)`; (d) whether any method
      contribution is still defensible or the next cycle tests a new hypothesis;
      (e) whether to run the shift-corrected decoder, with its success bar and
      kill criterion pre-registered before any result.
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
- Do not retry Claude or edit VPN/DNS configuration repeatedly while
  `docs/network_2026-11-15.md` describes an unconfirmed data-plane recovery.
