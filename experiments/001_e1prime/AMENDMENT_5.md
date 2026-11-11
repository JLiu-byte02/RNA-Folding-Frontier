# Amendment 5 — How baselines are trained, scored and compared (P5)

Status: written before any baseline has been run or scored, on the PI's Round 4e
ruling. It exists because the released parameters of the strong baselines make
the naive comparison uninterpretable on `withheld`.

## A5.1 — The problem this fixes

Released MXfold2 and BPfold weights were trained on TR0 and/or RNAStrAlign. Our
`S1'-withheld` split consists of clans that *are* present in bpRNA-1m and that we
removed from our own training. So on `withheld`, a released-weight baseline has
seen the held-out clans and our model has not. That direction is conservative
for us, but it makes the comparison uninterpretable as a claim about method
quality, which is what P5 needs.

## A5.2 — Rules, frozen

1. **Trainable baselines** — MXfold2 certainly; BPfold if its training code
   permits — are **retrained on our frozen train split** and scored that way for
   P5. A retrained baseline is the only comparison that can carry the claim.
2. **Released-weight scores** are reported separately, and each number is
   labelled with the overlap between the baseline's training data and our test
   clans. They are context, never the P5 comparison.
3. **Fixed-parameter methods** — RNAfold, RNAstructure, EternaFold — use
   parameters fit on other data and cannot be retrained on ours. They are
   reported as they are, with that caveat stated next to the number.

## A5.3 — What must be recorded for every baseline, before it is scored

1. Version or commit, the exact command line, and the licence.
2. What data the released parameters were trained on, with a citation.
3. Whether training code exists, and whether retraining on our train split is
   feasible (with the cost if it is not cheap).
4. A clan-overlap check: once `cmscan` is reachable, the fraction of each
   baseline's training sequences that fall in our withheld and novel clans.

Recorded in `docs/baselines.md`. A baseline without items 1–3 does not get a
number in the paper.

## A5.4 — What this does not change

H2.1–H2.4 involve no external baseline and are unaffected. The Turner-only
control on the local ArchiveII mirror (`experiments/004_turner_control/`) stays a
pipeline exercise and is not a P5 entry.
