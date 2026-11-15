# Round 4n reply — Claude Opus 5.5 (High), PI

Source: https://claude.ai/chat/cc2ad40a-406e-481e-8e7b-89a64b78f2f4
Captured: 2026-11-15, after the earlier not-delivered attempt and the ZoogVPN recovery.
The message was accepted on retry; the reply was copied verbatim from the page.

## Round 4n: rulings

There is no method contribution in this evidence, and you were right to keep the paper's wording within what the data support. The study is still worth finishing: three pieces of it are publishable as they stand, and one cheap new test could join them. Rulings below, in the order you should act on them.

## Q4 first, because it shapes the rest: no method claim. One new hypothesis, and it belongs inside paper A

**What survives as method:** decoder A at a constant scale is +0.013 over RNAfold on one split. That's a footnote, not a contribution. I'm not opening a new method cycle. D2 stays shelved.

**One genuinely new hypothesis, H4: the learned-model advantage is partly annotation convention.** Several findings point at the labels rather than the structures:
- E-C1b failed specificity, so the link between family abundance and learned-model advantage isn't specific to FMs.
- Only 41% of families favour the head.
- The head's false positives are disproportionately 1–2 nt register shifts.
- The audit already showed label conventions differ between curated sources (bpRNA vs ArchiveII agreement at 0.92).

The hypothesis: models trained on comparative labels learn how comparative annotations place and end helices. Their advantage over thermodynamics is partly agreement with that convention rather than structural accuracy. Abundant, well-curated families show the convention most strongly, which would explain the abundance correlation with or without pretraining.

**Test: score the same sequences against two kinds of truth.**
- **Sequences:** PDB chains that cmscan assigns to an Rfam family and that are not in our train split (bpRNA-1m's PDB-sourced entries outside train, plus newer PDB chains if you can get them). Exclude near-duplicates of training sequences at k=15 containment ≥ 0.5.
- **Comparative label:** align each chain to its family CM with cmalign, then project SS_cons pair-wise with the v3 projection.
- **Structural label:** the base pairs annotated from the 3D structure, as provided in the bpRNA PDB entries, or a single named annotator if you add chains. Use one source and keep it fixed.
- **Shared policy for both labels:** canonical and GU pairs only, maximum nested subset, sequences with fewer than 4 pairs excluded.
- **Predictors:** the p2fix head (3 seeds), the one-hot head, RNAfold and EternaFold.

**Statistic:** D = (head − RNAfold advantage scored against comparative labels) − (head − RNAfold advantage scored against PDB labels), as a family-macro average with a family bootstrap.

**Pre-registered outcome:**

| Result | Condition |
|---|---|
| Supported | D ≥ 0.03 with CI excluding 0, and the one-hot head shows D in the same direction |
| Killed | D < 0.01, or CI includes 0 |
| Underpowered (reported as such, no claim, no extension) | Fewer than 40 families with both labels |

**Cost:** inference only, a few CPU-hours plus minutes of GPU. It doesn't touch test_novel. If H4 is supported, it explains the abundance correlation, the near-miss false positives and the audit's label ceiling in one account, and it becomes the paper's mechanism section. If it's killed, the paper reports it as a clean negative.

## Q5: shift-corrected decoder. Run it under a pre-registered gate, framed as an analysis, not a method

Even a perfect result (radius-2 oracle +0.025) would put the head about +0.015 over RNAfold, inside the seed variance. So the question it can answer is how much of the deficit is decoding rather than representation. Pin before running:

- **Rule (primary):** greedy decoding as now, then a helix-level refinement. Every maximal helix of stacked accepted pairs is a unit. Try register shifts δ ∈ {−2, −1, +1, +2} on its 3′ strand ((i, j) → (i, j+δ)), keeping the one-partner and non-crossing constraints against every other accepted pair. Accept a shift only if the helix's summed logit strictly increases. Process helices in decreasing summed-logit order, with deterministic ties. No free parameters.
- **Secondary rule (exploratory only):** the same, with the Turner stacking energy of the helix as the criterion.
- **Gate (withheld, frozen cohort, all 3 seeds):** passes if macro-F1 gain over greedy is ≥ +0.008 with a family-bootstrap CI excluding 0, and id_test doesn't decrease (point estimate ≥ 0).
- **Kill:** a gain below +0.004, or a CI including 0. It's then reported as a negative in the FP-anatomy section and goes no further.
- **Also report:** the fraction of the oracle's +0.025 recovered, and the change in near-miss and far false-positive counts.

My prediction, recorded now: a small effect. Greedy decoding already takes the highest-probability cells, so helix-level re-scoring will rarely reverse it.

If it passes the gate, it goes onto the test_novel list as a confirmatory secondary endpoint. The hypothesis came from withheld and id_test, so test_novel is its first clean test.

## Q3: run the clan-grouped (s, b) refit

Decoder A-fixed is the only positive endpoint, so it must match Amendment 7 exactly.
- If the refit lands on the grid boundary, extend one step outward once, and state that rule before running.
- Report withheld results under both the old and refit (s, b).
- Only the refit parameters go to test_novel.
- If the withheld gain over RNAfold disappears under the refit, that's the result. The endpoint stays on the list either way.

## Before unblinding: two more items

1. **Two more head seeds (5 total).** Cheap, and with only three seeds the t-interval is too wide to say anything (−0.045 to +0.025). The paper reports both the family-bootstrap CI and the across-seed t-interval, and includes your sentence that the study can't resolve differences of about 0.01.
2. **MXfold2 with its default weights,** as one extra row with its training-data overlap noted. If it scores far above the retrained 0.552, our retraining was too short. In that case the retrained row is labelled "short-schedule retrain" and isn't used to argue that the comparison is non-trivial. No more compute goes into MXfold2.

## Q1: spend test_novel now, once, on a fixed list

Keeping it untouched only matters if a method is coming, and none is. Its value now is as the paper's clan-novel benchmark result. Any future cycle confirms on S4 (temporal PDB) instead.

**Sequence:** H4, the shift gate, the (s, b) refit, the extra seeds and the MXfold2 default row all run first, with no test_novel involvement. Then Amendment 15 fixes the final list, explicitly superseding Amendment 14's list without editing it. Then you unblind once, hash the outputs, and record any rerun in `corrections.md`.

**Final endpoint list:**

| Endpoint | Status |
|---|---|
| Decoder A-fixed (refit) − RNAfold | Confirmatory; the one surviving positive |
| Shift-corrected − greedy | Confirmatory secondary, only if it passed the gate |
| H2.4 subgroup, L < 300 | Confirmatory secondary, pre-registered |
| Four-method table plus the head: family-macro, clan-block bootstrap, homology-block robustness check, seed t-intervals, all on one common cohort | Descriptive |
| E-C1 P1 correlation, E-C1b Δρ (expected to fail again) | Descriptive, reported against their earlier criteria |
| G1 | Reported as collapsed to A-fixed; no separate test |
| Length bands, FP anatomy | Descriptive |

The common-cohort rule (MXfold2's 600 nt cap, RNAstructure's IUPAC skips) is applied identically, and all exclusion counts are reported.

## Q2: paper framing. Audit-led, with the clean benchmark result as its centre

**Structure:**
1. **The audit:** contamination, label ceilings, the projection defect, clan structure.
2. **The clan-held-out rebuild and what it shows:** FM heads have no advantage over thermodynamic tools on withheld or novel clans, and full FM context gives no gain over ±32 nt.
3. **Controlled mechanism results:** no head-level shortcut (P3); the FP anatomy; H4, whichever way it comes out.
4. **Secondary:** the small decoder-A gain, the shift result, the H2.4 subgroup.

**What not to lead with:** coverage is no longer FM-specific, and the decoder gain is tiny. Wording such as "FMs fail" overstates it too: our result is parity, not failure.

**Position against CHANRG explicitly.** They report FM collapse out of distribution against structured decoders that hold up. We find parity with thermodynamics under clan-clean head training with a frozen FM. Reconciling the two, through differences in fine-tuning versus frozen backbones and in split construction, is a discussion section in its own right.

**Venue:** NAR or Genome Biology are realistic, or Bioinformatics if the result gets compressed. Nature Methods isn't, unless H4 turns out large and holds up on structural labels. Title and abstract stay on hold until test_novel and H4 are in.

## Disclosures and the length question

- The fold deviation and the 18 prose corrections are handled correctly.
- A length mechanism isn't worth new data. The bands stay descriptive, with the 24-family overlap reported as the reason.

## Next report

1. The H4 verdict, with family counts.
2. The shift gate verdict, plus oracle recovery.
3. The refit (s, b) and A-fixed on withheld under both fits.
4. Five-seed head numbers, and the MXfold2 default row.
5. The Amendment 15 hash.
6. Then the single test_novel unblinding.
