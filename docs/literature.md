# Literature records — claims verified against full text

Each entry states what the paper actually says, where it says it, and what it
changes for us. Snippet-level notes live in `docs/literature_watch_*.md` and are
not citable.

## Flamm et al. 2022 — the length-scaling claim, and what they left open

Flamm C., Wielach J., Wolfinger M. T., Badelt S., Lorenz R., Hofacker I. L.
"Caveats to deep learning approaches to RNA secondary structure prediction."
*Frontiers in Bioinformatics* 2:835422 (2022). doi 10.3389/fbinf.2022.835422.
PMC9580944 (full text read 2026-11-12 via the Europe PMC XML).

### What they report

- **Abstract, verbatim:** "the number of predicted base pairs scales
  quadratically with sequence length, even though a secondary structure can only
  accommodate a linear number of pairs."
- **Figure 3** (results, "Predicting Pairedness" / output-matrix section): for
  2,000 sequences per length bin over 30–250 nt, the number of entries of the
  output matrix **above 0.5, i.e. before any postprocessing**, grows quadratically
  for all three network models, while the RNAfold ground truth grows linearly.
  Their explanation: the matrix has n² entries and the network predicts a
  constant fraction of all possible pairs.
- **Discussion, verbatim:** "While this behavior can easily be addressed during
  postprocessing, it is not clear whether that would correct or merely hide the
  underlying problem."

### What their postprocessing was

They deliberately avoided sophistication: keep the highest entry per row, round
to 0/1, remove pseudoknots. That is a one-partner, non-crossing rule with **no
energy term and no loop model** — the same family as our greedy nested decoder.
They did not test an energy-based decode.

### Consequences for our claims

1. **H2.1 is anticipated.** The phenomenon we measure (long-range FP rate per
   sequence rising with length) is the same one they report, with a better
   controlled generative process. H2.1 is replication; cite Figure 3 for the
   pair-count scaling, not for our statistic.
2. **H2.4 is a question they explicitly pose but do not answer.** They say
   postprocessing "can easily be addressed" but may "merely hide" the problem,
   and they did not run an energy-based decoder. Our restated H2.4
   (Amendment 3 §A3.1) is therefore the answer to their open question, not a
   rediscovery — *provided* we measure both the simple one-partner decode
   (which we do) and the energetic one at matched recall.
3. **Framing for the paper.** "Deep predictors over-produce long-range pairs"
   cannot be our finding. "Energetic inference over FM pair evidence changes
   which errors survive, measured at matched recall" can be, if it holds.

### Still to check

- Whether Flamm et al.'s follow-ups (or the `rnadeep` codebase) test a
  learned-score + exact-energy DP. The paper itself does not.
- Whether their Figure 3 statistic (matrix > 0.5, all cells) and ours
  (per-sequence long-range FP count) can be related quantitatively; if yes, we
  can reproduce their curve on our splits as a sanity check of the whole
  pipeline.
