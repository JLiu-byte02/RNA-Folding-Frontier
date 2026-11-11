# Literature Watch

Living document. Add papers with: citation, venue/year, core idea, relevance,
and a concrete implication for this project.

## Seed areas

- RNA foundation models: RiNALMo, RNA-FM, RNA-MSM, ERNIE-RNA.
- Structure prediction: bpRNA / ArchiveII benchmarks, UFold, SPOT-RNA,
  E2Efold, RNAfold / LinearFold baselines.
- ML mechanisms worth transplanting: structured decoding, energy-based models,
  diffusion on graphs, equivariant/geometric priors, test-time adaptation,
  uncertainty calibration, retrieval-augmented prediction.

## Entries

### BPfold — base-pair motif energy learned from a thermodynamic library

- Zhu et al., *Nature Communications* 16:5856 (2025). Code MIT.
- **What they claim.** A base-pair motif energy (BPM) map predicted from sequence
  by a CNN, combined with thermodynamic energy minimisation, generalises better
  than end-to-end structure predictors on unseen families.
- **What their training set actually is (verified from the paper text).**
  Trainable models are trained on **RNAStralign (19,313 sequences, <= 600 nt) and
  bpRNA-1m** (MXfold2's TR0, 12,114 sequences), jointly. RNAStralign covers eight
  RNA families: 5S, 16S, 23S, group I intron, group II intron, RNase P, SRP,
  tRNA — i.e. eight of ArchiveII's nine families.
- **Their evaluation framing (verified verbatim).** "ArchiveII (3966 RNAs) and
  bpRNA-TS0 (1305 RNAs) are sequence-wise datasets while Rfam12.3–14.10 (10,791
  RNAs) contains cross-family RNA sequences". ArchiveII is called "the most
  widely used benchmark dataset for evaluation of RNA secondary structures", not
  an OOD set. Their cross-family claim rests on the Rfam 12.3→14.10 delta, which
  is a **temporal** holdout (entries added between releases), not a clan-held-out
  split.
  **Correction (same day, after re-reading the text).** The *dataset* is the
  Rfam 12.3→14.10 delta, but the *protocol* is described as "family-wise
  cross-validation" and the headline results are on named unseen families
  (Skip-rope n=240, Twister n=236, C-di-G n=197, RAGATH n=160, Cobalamin n=286).
  So it is a family-level holdout over the delta, not merely a temporal one. Two
  flaws remain: it is not clan-aware (a new delta family that is a new member of
  an old clan leaks), and the accessible main text does not state how homology
  between the held-out family and the training pool was screened. Our S1′ is
  therefore a clan-level strengthening of BPfold's protocol, not a different
  protocol — the novelty must come from the mechanism and the payoff, not from
  the split itself.
- **Implication for us.** (i) Do not attribute an OOD claim about ArchiveII to
  BPfold; the honest statement is that ArchiveII is routinely *reported* as a
  generalisation benchmark while being in-distribution under BPfold's protocol.
  (ii) ArchiveII's 20.47% exact overlap with bpRNA TR0 is only half of the
  contamination; RNAStralign is the other trainer and must be audited too.
  (iii) The Rfam-delta cross-family set is a fair comparison target for our
  payoff experiment, but it must itself be checked for clan overlap — new Rfam
  families are frequently new members of old clans.
- **Release/paper discrepancy.** The released TR0 has 10,814 sequences; the paper
  reports 12,114. The released pool is split into TR0 (10,814, train) and VL0
  (1,300, validation), and 10,814 + 1,300 = 12,114. So the release re-splits the
  paper's TR0. Any reproduction table must say which convention it uses.

