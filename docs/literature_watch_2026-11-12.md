# Literature watch — 2026-11-12

Scope: a targeted web scan while the GPU server is unreachable. **These are
snippet-level reads.** No full text, abstract page, or supplementary file has
been downloaded or verified, so nothing here may be quoted as a finding of those
papers until the PDFs are read. What follows is what the snippets say and what
each item changes for us.

> Update: the Flamm et al. 2022 item below has since been read in full. Use
> `docs/literature.md` for its verified claims (Figure 3, abstract quote) and
> treat this file's version as superseded for that paper.

## 1. H2.1 is anticipated — Flamm et al. 2022

**Flamm, Wielach, Wolfinger, Badelt, Lorenz, Hofacker.** "Caveats to deep
learning approaches to RNA secondary structure prediction." *Frontiers in
Bioinformatics* 2:835422, 2022. DOI 10.3389/fbinf.2022.835422.
Code/resources: https://github.com/ViennaRNA/rnadeep

Verified from the abstract/author list (PubMed 36304289, PMC9580944): the paper
argues that ML models must be shown to learn the *biophysical model*, not just
fit phenomenologically, and tests several models on synthetic data from a
controlled model so that dataset biases can be removed. **To check in the full
text, not yet read:** the specific claim I expect to find — that models
overpredict base pairs as sequence length grows and depend heavily on the
structural families seen in training.

**Consequence.** H2.1 ("long-range false positives rise with length") is not a
new observation; it is a 2022 result obtained under stronger control than we
have (a known generative model). We must cite it and must not present H2.1 as
novel. What can still be ours:

- the **causal intervention** (P2-block: rewrite only the distal content and
  watch the *change* in long-range FPs), which is a different question from
  "long inputs produce more pairs";
- the **decode-time mechanism** (H2.4: does an exactly consistent decoder remove
  them), because a structural prior is a different lever from an architectural
  one;
- any **performance** consequence, which Flamm et al. do not report for a
  pretrained-backbone setting.

If the PI wants an H2-level novelty claim, the paper to beat is this one, not the
DL-prediction leaderboard.

## 2. UFold-X — a new strong baseline we are not currently comparing against

"UFold-X: an enhanced Dual & Dynamic U-Mamba model for long-range RNA secondary
structure prediction." *Nucleic Acids Research* 54(17):gkag887, 23 Sep 2026.
https://doi.org/10.1093/nar/gkag887 (author list not yet captured — check the
PDF before citing)

Abstract-level content: two branches — convolutional for local, Mamba for
long-range — with a length-adaptive gate; evaluated on RNAs up to 5,000 nt;
reports a **cross-clan benchmark for long RNAs** and explicit long-range base-pair
evaluation; claims to be comparable to state-of-the-art methods overall and to
outperform other DL methods on cross-clan generalisation. The exact benchmark
names in my first search snippet (RNAstralign-1800, BGSU1800, Rfam1800) are
**unverified** and must be read from the paper.

**Consequence.** Two concrete actions:

1. add UFold-X to the baseline list (weight download and licence check on the
   server), because a reviewer will ask;
2. read their cross-clan protocol — if it is usable, it is an external benchmark
   for our "cross-family generalisation" claim, which currently rests on our own
   S1′ split only.

## 3. structRFM is now a Nature paper

Zhu Heqin (first author), Zhou Shaohua and Xiong Peng (corresponding), USTC
MIRACLE lab. "An open structure-guided RNA foundation model for generalizable
prediction." *Nature*, online 4–5 Nov 2026. Verified via the USTC
announcement and Xinhua coverage: ~21 million computationally annotated
sequence–secondary-structure pairs, pairing information injected into
pretraining, open release, feature levels at sequence/nucleotide/pair.

Not yet verified: the checkpoint licence, the maximum input length (the
existing repo note says 512 nt), and checkpoint size. Check the repository
rather than repeating the note.

**Consequence.** The repo already flags structRFM as "MIT, contamination-flagged";
that flag now matters more, not less, and if the 512 nt input limit is real it
interacts directly with our length-band findings. A backbone comparison that
includes structRFM must state the truncation policy explicitly and check
whether the claimed 21M training pairs overlap our test clans.

## 4. Others noted, not yet acted on

- **DRfold2** (*PLoS Biology*, 2026 — exact volume/DOI not captured): RNA
  *tertiary* structure; relevant only as context for why secondary structure
  still matters.
- **arXiv 2511.02622**, "Machine Learning for RNA Secondary Structure Prediction:
  a review of current methods and challenges": useful for the related-work
  section and for its framing of the data bottleneck; check its coverage of
  decoder/consistency methods before writing ours.

## What is still unverified

Every claim above is from a search snippet. Before any of it enters a report or
a paper draft: obtain the PDF, read the methods, and record the specific claim
with a page/section. In particular the Flamm et al. numbers (how many extra
pairs, at what lengths) must be read from the paper, not inferred.
