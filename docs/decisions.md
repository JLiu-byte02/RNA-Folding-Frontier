# Standing decision rules

Rules that constrain how the project is run, as opposed to what it studies.
`docs/decision_log.md` records individual decisions; this file records the rules
those decisions must obey. Rules are added, never rewritten.

## Amendments are append-only as files

A frozen amendment is a document whose sha256 we have quoted to each other. If it
can be edited afterwards — including by appending — the hash stops proving
anything about what was pre-registered.

- Once an amendment's hash has been reported (to the PI, in a commit message, or
  in `state/amendments.lock.json`), **the file is immutable**.
- Any addition, correction, or clarification goes into a **new** amendment file.
- Superseded wording stays visible in git history. Never rewrite a quoted file to
  make the record look cleaner.
- `state/amendments.lock.json` holds the authoritative hash per amendment; the
  embargo guard (`src/rff/analysis/embargo.py`) refuses confirmatory work when a
  file no longer matches.

**Origin:** 2026-11-12. Amendment 4 was appended to (new §A4.6) after its hash
`08d6771a…` had been quoted, which moved it to `c6471002…`. PI ruling (Round 4f,
Fix 1): restore the original bytes, move §A4.6 into Amendment 6, keep the history.
Done — `AMENDMENT_4.md` hashes to `08d6771a…` again and `AMENDMENT_6.md` carries
the provenance procedure as §A6.1.

## 2026-11-12 — Round 4 P2 verdict withdrawn

The Round 4 P2-central table and the “P2 inverts the prediction” verdict are
withdrawn. They were computed with `pair_prf`, which scored every zero-TP
sequence as F1=1.0. The original family-shortcut test is unanswered; it is not
reversed. Evidence: `docs/corrections.md` and the corrected D1/K1 rerun.

## 2026-11-12 — H2 retired as a mechanism

H2 ("distal context is a net-negative prior") was inferred from the voided P2
table. After the scorer fix, masking the distal region (condition e) is about
0.30 F1 worse than native a, not better. The premise is false; H2 is retired as
a mechanism.

## 2026-11-12 — H2.4 continues as a standalone decoder question

H2.4 — does energy-consistent decoding cut long-range false positives versus
greedy nested matching on fixed FM pair scores? — never depended on “e helps.”
It remains a standalone question, motivated by H2.2 and the helix-energetics
check. It is not launched until the corrected P2-central table has been seen by
the PI.

## 2026-11-12 — Round 4k: P2-block becomes the decisive D1 test

P2-central is exploratory because one forward pass cannot give every position
its own native flank. The decisive test is **P3 train-and-test**: train separate
heads on condition a, an e-block cache, and a b-block cache, then evaluate each
on the same condition in withheld. The D1-core rule, the H3 test-time estimand,
the c dose-response bins, and the `(s,b)`/H2.4 order are frozen in
`experiments/001_e1prime/AMENDMENT_10.md` before any block cache is read.

The near-duplicate positive control is deferred until after P3, not cancelled.
`test_novel` remains embargoed until Amendment 10's outcomes are recorded.
