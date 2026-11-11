# Round 4d — PI rulings (Claude Opus 5.5 High, 2026-11-12)

Structure of the reply, condensed to the binding content. The full exchange is
in the Claude conversation `cc2ad40a-406e-481e-8e7b-89a64b78f2f4`.

## Novelty: narrow the claim before it hardens

"'No one tests a learned score plus an exact energy DP' is true of Flamm et al.'s
follow-ups but false for the field: MXfold2 combines learned scores and Turner
terms in an exact DP, and BPfold puts energy inside its model. Both train the
score and the decoder together, so neither isolates the decoder. The defensible
framing of H2.4 is a controlled decoder comparison on fixed FM pair scores: does
energy-consistent decoding remove the length-scaling false positives Flamm et al.
identified, compared with greedy nested matching? Their discussion leaves that
open."

Recorded in `docs/literature.md` next to the Flamm entry.

## Q1: recall-matched first, then pair-budget, then the recall ceiling

Freeze this order in Amendment 4:

1. **Primary: recall-matched.** One greedy threshold per split and seed, a single
   global value set so pooled recall equals decoder A's, choosing the closest
   point at or above A's recall if no exact match exists. "Choosing it on the
   evaluation data favours greedy, which is conservative toward the null; state
   that."
   **Verdict:** supported if A has fewer long-range false positives
   (`|i-j| > 128`, pooled) **and** at least as many true positives.
2. **Fallback, only if no threshold reaches A's recall: pair-budget matched.**
   Choose the greedy threshold so its total predicted pair count equals A's.
   **Verdict:** supported only if A has more true positives **and** fewer
   long-range false positives at equal count. "This is a stricter test than
   option (iii), because A must win on both counts at once."
3. **Always reported, descriptive only:** greedy's recall ceiling, meaning its
   maximum pooled recall over all thresholds.
4. **If neither match is feasible:** report "inconclusive" for that split, length
   bin or seed, and report how often that happens. "Don't add another fallback
   after seeing data."

"Run all of this per length bin too, but the verdict is the pooled result."

## Q2: both, as separate variants with separate jobs

- **A-energy-fixed** — the confirmatory variant for H2.4, uses option (a).
  Same BCE-trained head as greedy. Pseudo-energy `dG_ij = -(s*logit_ij + b)`
  with a global scale and offset. "The offset is needed because pos_weight 20
  shifts the logits by about log 20. Without b, the scale would also be absorbing
  calibration." Fit `(s, b)` by 5-fold cross-validation on the frozen train
  split, maximising family-macro F1. Freeze and hash them into Amendment 4
  before scoring id_test or withheld again. "The decoder is then a fixed,
  monotone function of the same ranking that greedy sees. That isolates
  'energy-consistent decoding vs greedy', which is the claim." Also report a
  sensitivity curve of F1 against s on id_test as exploratory.
- **A-energy-learned** — the P5 variant, uses option (b). Head trained with
  max-margin through the decoder, bonus tanh-bounded at a cap C kcal/mol. Fit C
  on train-split cross-validation over {1, 2, 3, 5} and freeze it. "This is the
  method we put up against BPfold. Its advantage over A-energy-fixed is
  end-to-end calibration, and we report it as such, not as evidence for H2.4."

Amendment 4 must state that Amendment 3's "A-energy" means A-energy-fixed, that
A-energy-learned is for P5, and give the fitting procedure for both.

## Q3: yes to a local mirror, but rebuild the frozen splits rather than the distributed files

"ArchiveII and TR0 'as distributed' exercise the pipeline, but they aren't our
evaluation data, and ArchiveII is contaminated by design."

1. **Membership lists.** If they aren't in the repo, commit the sequence-ID lists
   for the frozen splits from the last server snapshot, or reconstruct them from
   committed code and seeds. "If neither is possible, say so; that's a
   reproducibility gap we have to close as soon as the cluster returns."
2. **Rebuild locally from public sources:** the bpRNA-1m HuggingFace mirror at
   the pinned revision, and Rfam 15.0 `Rfam.seed.gz`. Check its md5.
3. **Check every split's sha256 against the frozen hashes.** "A mismatch is a
   finding, not something to work around."
4. **If the hashes match,** run the CPU baselines locally on id_test and withheld
   only: RNAfold (also the Turner-only control), EternaFold, RNAstructure and
   MXfold2, and BPfold if CPU speed allows. Same scorer as our models.

"ArchiveII as distributed is still needed for T4 reproduction, so mirror it as
well with a checksum, but it comes after the split rebuild."

## Queue when the cluster returns (unchanged)

D1/D2 dumps → P2-block cache → train-and-test heads → positive control →
A-energy-fixed (`(s, b)` fit, then H2.4 on id_test and withheld) →
A-energy-learned training.

"Lead the next report with K1's status, then whether the split hashes reproduced
locally."
