# Round 4n — draft message to the PI (2026-11-14)

Status: not sent. Computer Use was stopped by policy on 2026-11-14 because the
browser URL could not be determined with enough confidence. Paste the body below
into the Claude chat (Opus 5.5, High) when a permitted browser path is available.

---

Round 4n — E-C2 closes negative, the four-method cohort is frozen, and two audits changed the manuscript

**1. E-C2 (A13): fail, with one honest nuance.**
The gate collapsed to its boundary (`s_max=.1875`, `alpha=beta=gamma=0`): the
fitted gate is a constant scale, bit-identical to its own control, and G1−fixed
reverses on id_test (`−.0064 [−.0122,−.0002]`). A13.2 does not pass. The one
arm that beats a thermodynamic tool is decoder A at that frozen constant scale:
withheld `.6137` vs RNAfold `.6003`, paired `+.0134 [+.0025,+.0241]` on 15,533
sequences; id_test `+.0204 [−.0020,+.0392]`, i.e. directionally consistent but
not confirmed, and below the pre-registered +.02-over-best-baseline bar. I have
reported it as a small decoder-level gain over RNAfold on one split, not as a
method advantage, and narrowed Paper A's "no advantage" to the default greedy
decoder.

**2. The four-method cohort is frozen (retrained MXfold2 included).**
Cohort: withheld 16,040 seqs / 408 families, id_test 1,272 / 153; the 600 nt
MXfold2 cap binds (597/572 nt), costing 19 + 15 sequences plus RNAstructure's
104 IUPAC skips on withheld. All rows share one id set, `missing_prediction=0`.

| method | withheld macro [95% CI] | id_test macro [95% CI] |
|---|---|---|
| RNAfold | .6003 [.5736,.6266] | .5997 [.5559,.6427] |
| EternaFold | .6049 [.5792,.6301] | .5952 [.5522,.6374] |
| RNAstructure | .5900 [.5634,.6165] | .5908 [.5468,.6343] |
| MXfold2 (retrained) | .5518 [.5258,.5772] | .5491 [.5048,.5928] |
| head 3-seed mean | .5905 [.5671,.6139] | .5901 [.5512,.6275] |

Paired: head_mean−RNAfold `−.0098 [−.0301,+.0107]` / `−.0097 [−.0384,+.0187]`;
−EternaFold `−.0144 [−.0341,+.0052]` / `−.0052 [−.0304,+.0211]`;
−RNAstructure `+.0005 [−.0193,+.0208]` / `−.0008 [−.0289,+.0272]`;
−MXfold2 `+.0387 [+.0225,+.0549]` / `+.0409 [+.0160,+.0663]`.
The three-tool verdict is unchanged. **Caveat:** the retrained MXfold2 is a
10-epoch MixC trained on our split with nested labels under a 600 nt cap; .5518
is well below the thermodynamic tools. The head's +.039 margin over it is not a
competitiveness claim, only evidence that the comparison is not vacuous. If you
want MXfold2 to be a meaningful baseline, it needs a longer schedule or the
published recipe, and that is a compute decision for you.

**3. Where the head loses (new, artifact-backed).** Two stratifications, both on
the frozen cohort:

*Length (pre-specified bands, withheld family-bootstrap within band).* The head
is at parity or slightly ahead below 150 nt and loses above:

| band | n seqs / fam | head−RNAfold | head−EternaFold | head−RNAstructure |
|---|---:|---:|---:|---:|
| 0–80 | 3,963 / 181 | +.010 [−.028,+.049] | +.018 [−.016,+.053] | +.018 [−.020,+.056] |
| 80–150 | 8,051 / 260 | +.016 [−.008,+.039] | +.007 [−.016,+.030] | +.024 [+.001,+.047] |
| 150–300 | 3,756 / 74 | −.036 [−.071,+.001] | −.073 [−.102,−.044] | −.029 [−.064,+.007] |
| 300–597 | 270 / 19 | −.038 [−.084,+.000] | −.092 [−.133,−.054] | −.049 [−.091,−.016] |

id_test reproduces the ordering with 28 and 11 families in the upper bands (all
CIs cross zero). Since 75% of the cohort is below 150 nt, the aggregate is
dominated by the short regime.

*Family size.* On withheld the median head−RNAfold delta is `−.021 [−.035,−.011]`,
worse than the mean (−.010); only 41% of families favour the head. The delta
rises with log family size (rho `.20 [.10,.29]` vs RNAfold; id_test
`.36 [.24,.47]`) and has no relation to mean family length. This is consistent
with the E-C1 coverage account, not an independent test of it.

The two are not independent — larger families are also longer — so I would not
put both in the abstract. The length crossover is the one a reviewer will ask
about, because it means "no advantage" is really "advantage below 150 nt,
deficit above".

**4. Two disclosures that belong in the record.** (a) `fit_energy_calibration.py`
used record-level folds (`crc32(id) % 5`) instead of the clan-grouped folds
Amendment 7 specified. It affects the `(s,b)` fit only — the E-C2 gate fit does
use clan folds — and the selected point is at the grid boundary, so no
conclusion is rescued by it, but the fold F1 behind `(s,b)` is optimistic
relative to the pre-registered design. A clan-grouped refit is one 5-fold grid
fit on train; say the word. (b) Two read-only audits of Paper A found 18 places
where prose did not match the code or the artifacts (head architecture, decoder
complexity, ViennaRNA pseudo-energy sign, six mislabelled Result populations,
two misquoted ranges, the "14 amendments" count). All are corrected in
`docs/corrections.md` and the draft; no metric changed.

**5. Paper A.** Sections 1-4 and Tables T1-T7 are drafted from committed
artifacts; title, abstract, discussion and both conclusions remain HELD. Claim
9 is no longer "blocked": it now states a valid single-cohort comparison with no
advantage over the thermodynamic tools and a weaker retrained MXfold2.

**Questions for you.**
1. Amendment 14 one-shot `test_novel`: do we spend it now, and on which
   endpoints? Given the four-method result, the only surviving positive is the
   decoder-A constant-scale arm and the coverage correlation.
2. Paper A framing: audit-led (C1), negative-led (C2), or "coverage mechanism +
   small decoder gain, no head advantage"?
3. Clan-grouped `(s,b)` refit: run it to complete Amendment 7 literally, or keep
   the disclosed deviation?
4. Direction: with D1-core falsified, H3 unsupported, E-C1b failing specificity
   and E-C2 failing, is any method contribution still defensible from this
   evidence, or should the next cycle test a genuinely new hypothesis? If the
   latter, name the hypothesis and its kill criterion.
