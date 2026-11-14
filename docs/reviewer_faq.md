# Reviewer-objection checklist

Prepared 2026-11-15 from committed artifacts, for the PI's framing discussion.
Each answer points at the file a reviewer can inspect; none of it is a new
claim.

**1. "The head loses to RNAfold and EternaFold — why is this publishable?"**
The contribution is the audit plus negative mechanism results, not a
state-of-the-art claim. The four-method cohort is frozen, every row shares one id
set, and the paired family bootstrap shows no advantage (Table T7). The one
positive mechanism, coverage, is reported as a shared predictor because the
one-hot control fails the specificity bar.

**2. "The cohort stops at 597 nt; long RNAs are missing."**
The binding cap is the retrained MXfold2's 600 nt; 19 withheld and 15 id_test
sequences are lost, and RNAstructure's 104 IUPAC skips are also inside the shared
exclusion set. Every method is scored on the same ids, and the limit is stated in
section 4 and section 7 of the baseline document.

**3. "Your retrained MXfold2 is far worse than published MXfold2, so the
comparison is meaningless."**
It is a retraining control on our split (MixC, 10 epochs, nested labels, 600 nt),
not the published system, and the paper says so. Its margin is never used as a
competitiveness claim; its role is to show the comparison is not vacuous. A
longer schedule is an open compute decision for the PI.

**4. "Three seeds cannot support a no-advantage claim."**
Correct, and the draft says it: the seed-level SD is .014, larger than the ~.01
effect, with a t interval of [-.045, +.025]; Table T7's family CIs condition on
the three seeds and understate total uncertainty (Limitation 8). We claim neither
an advantage nor a tie.

**5. "The length bands suggest the head fails on long RNA."**
We do not claim that. Only 24 withheld families (5 on id_test) have sequences in
both a short and a long band, and the within-family interaction crosses zero
(head-RNAfold +.011 [-.046, +.067]; head-EternaFold +.040 [-.007, +.087]). The
bands are descriptive; family identity and length are confounded.

**6. "v3 labels are your own construction; why trust them?"**
v3 parses 92,781 of 92,781 rows. On the 83,459 rows where both v2 and v3 parse,
the pair sets disagree on 8.1%. Against bpRNA-1m pair sets, mean pair F1 on
37,797 exact-sequence comparisons is .9056 for v3 versus .8148 for v2.
Artifacts:
`experiments/011_labels_v3/labels_v3_summary.json` and
`v2_v3_bprna_crosscheck.json`.

**7. "The overlap percentage uses 3,966 records, but the label mirror has
3,370."**
Both populations are named wherever they are used, and neither number is mixed
into a rate with the other denominator. Source:
`experiments/000_e0_leakage/union_overlap/summary.json` and the
`label_discordance_*` summaries.

**8. "A scorer bug was found — are the numbers trustworthy?"**
Every number in the draft postdates the fix (Amendment 9). An independent scorer
cross-checks the counts, `tests/test_scoring.py` pins the zero-TP convention, and
`docs/corrections.md` lists each affected number with its replacement.

**9. "The energy decoder failed only because the grid was constrained."**
The coarse grid imposed `beta >= 0`, so negative dependence is untested rather
than excluded, and the draft states that. Inside the tested region the gate
collapses to a constant scale that is bit-identical to its control and reverses
sign on id_test.

**10. "Why not unblind test_novel?"**
Amendment 14 fixes a one-shot list; the PI decides when to spend it. Until then
only labels and counts are inspected, no scoring output exists anywhere, and none
of the reported results depend on test_novel.

## Reproducibility

`tests/test_paper_numbers.py` re-parses the paper tables and fails on drift;
`tests/test_family_diagnostics.py` regenerates the family diagnostics
byte-for-byte; `tests/test_paper_figures.py` re-renders the figures;
`tests/test_public_repo_hygiene.py` blocks credentials. The fresh-clone CI runs
all of them against exactly what a reader sees.
