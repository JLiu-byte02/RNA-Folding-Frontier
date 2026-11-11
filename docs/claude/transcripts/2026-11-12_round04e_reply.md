# Round 4e — PI rulings (Claude Opus 5.5 High, 2026-11-12)

Condensed to binding content. Full exchange: conversation
`cc2ad40a-406e-481e-8e7b-89a64b78f2f4`.

## 1. A provenance contradiction that comes first

"Your Round 3 report said the frozen split hashes were in the repository, and
Round 4 quoted a test_novel hash (9d1ab36c…). Neither is in the pushed repo.
Either they're in unpushed commits or untracked files on the server, or the
Round 3 statement was wrong. Without them, the embargo on test_novel currently
rests on your word, not on an artifact."

Resolution, to be added to A4.5 and done **before any new computation** when the
cluster returns:

1. `git status` and `git log origin/main..HEAD` on the server, then push
   everything, including membership lists (sequence IDs per split) and per-split
   sha256 hashes.
2. Record file modification times for the split files and the test_novel freeze
   file, plus confirmation that no scoring output exists for test_novel.
3. Determinism check: rebuild every split from committed code, seeds and pinned
   inputs into a fresh directory, compare hashes with the on-disk originals.
   Match → the freeze is credible after the fact. Mismatch → stop and report.

## 2. ArchiveII mirror (3,370 vs 3,966)

"Don't quote anything from it." Diagnose with a per-family count table and a
length distribution against the cmscan input when there is slack. The gap is too
large for the 600-nt cap (136 records) and does not match the duplicate count
(3,966 − 515 = 3,451), so it is probably a differently processed distribution;
record which one.

## 3. (a) RNAfold on the local mirror — yes, and cross-check the scorer

"A scorer bug would silently corrupt every number from here on, and this is the
cheapest point to catch one." Score the same predictions with our scorer and one
independent implementation and require exact agreement per sequence. Include
these as unit tests: empty truth / empty prediction / both; pseudoknotted true
pairs under the pre-registered pseudoknot policy; non-canonical true pairs; the
±1 slip variant, reported separately and never as the primary.

## 4. (b) Baseline documentation — yes, and record what each model was trained on

Released MXfold2 and BPfold weights were (per the PI's premise) trained on TR0
and/or RNAStrAlign, while S1′-withheld consists of clans present in bpRNA-1m that
we removed from our own training. So released-weight baselines have seen the
held-out clans and our models have not: conservative for us, but the P5
comparison becomes uninterpretable.

For each baseline record: (1) version or commit, exact command line, licence;
(2) what the released parameters were trained on, with a citation; (3) whether
training code exists and retraining on our train split is feasible; (4) a
clan-overlap check once cmscan is reachable.

**Amendment 5 rule, committed before any baseline is scored:**

- Trainable baselines (MXfold2 certainly; BPfold if its code allows) are retrained
  on our frozen train split and scored that way for P5.
- Released-weight scores are reported separately, labelled with their overlap
  with our test clans.
- RNAfold, RNAstructure and EternaFold use fixed parameters and are reported as
  is, with that caveat.

## 5. Queue when the cluster returns

Provenance resolution → D1/D2 dumps, then K1 → P2-block cache → A-energy-fixed
(fit `(s, b)`, then H2.4 on id_test and withheld) → train-and-test heads →
positive control → baseline retraining → A-energy-learned.

"Lead the next report with the provenance check, then K1."
