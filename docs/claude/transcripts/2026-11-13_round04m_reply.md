# Round 4m reply — Claude Opus 5.5 (High), PI

Source: https://claude.ai/chat/cc2ad40a-406e-481e-8e7b-89a64b78f2f4
Captured: 2026-11-13, after E-C1 passed.

## Round 4m: E-C2 authorized, but one cheap specificity control has to run first

E-C1 passing on both proxies, with the analysis pre-registered and the
population frozen, is the first positive result this project has had. Before it
carries the name "pretraining coverage", though, it needs one more control.

### The open question: is the effect really about pretraining?

The outcome was ΔF1 = (FM head) − (RNAfold). There are two ways this correlation
with family abundance could appear without pretraining having anything to do
with it:

1. **The RNAfold side moves, not the FM side.** Abundant, well-studied families
   could be harder for thermodynamics, or be labelled in a more comparative
   style that any model trained on comparative labels would match.
2. **Any learned model would show it.** Abundant families may have cleaner,
   more consistent `SS_cons` labels, which benefits every model trained on
   bpRNA, pretrained or not.

Pretraining coverage predicts something specific: the relation is present for
the FM head and absent, or much weaker, for an identical head with no
pretraining. That is the control.

### E-C1b: the specificity test (pre-register in Amendment 13 before running)

**Decomposition** (descriptive, on the same 409 families):

- partial ρ of FM-head F1 with P1;
- partial ρ of RNAfold F1 with P1.

If most of the ρ = 0.34 comes from RNAfold doing worse on abundant families
rather than the FM head doing better, C is about thermodynamics, not pretraining.

**Matched no-pretraining control (primary):**

- Same decoder-B architecture and training recipe, with the FM embeddings
  replaced by learned one-hot nucleotide embeddings, trained on the same train
  split with 3 seeds. That's minutes per head.
- Compute `ΔF1_onehot = (one-hot head) − (RNAfold)` with the identical partial
  correlation.

**Criterion:** `Δρ = ρ_FM − ρ_onehot ≥ 0.15`, with a paired family-bootstrap CI
excluding 0.

- Passes: "pretraining coverage" is earned.
- Fails: the finding is reported as "learned-model advantage tracks family
  abundance", which is still true but weaker, and we don't attribute it to
  pretraining.

**Secondary control:** MXfold2 retrained on our split, run through the same
analysis. Run E-C1b on test_novel in the single unblinding as well.

### E-C2: authorized, in parallel, under Amendment 13

E-C2 only claims a performance improvement, so it can start before E-C1b
finishes. Its interpretation depends on E-C1b. Pin the following in Amendment
13 before any fit.

**Gate inputs:** only quantities computable from the sequence itself. Use
per-sequence pseudo-perplexity, not family medians and not Rfam counts. A
deployable decoder can't look up family abundance for a truly novel RNA.

- **G1 (primary, the only variant that can carry the confirmatory claim):**
  `s(x) = s_max · σ(α − β · z(pppl(x)))`.
- **G2 (secondary):** G1 with an extra `−γ · log L` term inside the sigmoid.
  This is motivated by the H2.4 length bins, which came from withheld, so it
  stays secondary.

**Fitting:** a grid over `(s_max, α, β)` using the clan-grouped 5-fold CV on
train from Amendment 7. Freeze and hash before scoring withheld. Compute train
pseudo-perplexity with the same 4 masks and seeds as E-C1.

**Success:** on withheld, G1 beats fixed-s decoder A and the best clean baseline
by ≥ 0.02 family-macro F1, with a CI excluding 0. Clean baselines are RNAfold,
EternaFold, RNAstructure and MXfold2 retrained on our split. BPfold is reported
but flagged, because TR0 contains our withheld clans, and it can't define
"best". The result then has to be confirmed on test_novel.

**Expectation stated before the run:** the mean FM advantage over RNAfold is
+0.005, and EternaFold usually beats RNAfold. So clearing +0.02 over the best
clean baseline is a long shot. Also report decoder A-energy-fixed's
family-macro F1 against RNAfold and EternaFold on withheld; it's needed for the
table either way.

**Also pre-register in Amendment 13:** a secondary confirmatory H2.4 subgroup:
for `L < 300`, decoder A reduces long-range FPs and increases TPs at the
pair-budget point. This hypothesis came from the withheld bins, so test_novel
is its only clean test. The ≥300 nt failure goes in as a known limitation.

### MXfold2 retraining: start now, as the top compute item

It's now both the last baseline needed for the unblinding and E-C1b's secondary
control. It's CPU-heavy (the DP is in C++), so submit it as a many-core Slurm
job today, with the train-split hash, MXfold2 commit and hyperparameters
(README defaults) recorded in the manifest.

### Paper A: draft in parallel, structured so C can go in or stay out

Write now:

- **methods:** splits, the v3 projection, the scorer convention including the
  corrected internal bug, pre-registration history;
- **the audit section (F1–F3);**
- **the negative mechanism results.**

Add P3 as stated: for this head, full FM context gives no gain over ±32 nt, in
distribution (`Δ −0.0095 [−0.025, +0.006]`) or OOD. Report the H2.4 length split
descriptively.

**Hold the title, abstract and discussion** until E-C1b, E-C2 and test_novel
are in. The headline is one of three:

1. pretraining coverage plus a gated decoder;
2. pretraining coverage as a mechanism only;
3. audit plus well-controlled negatives.

Leave a placeholder section for C with both wordings written in advance.

### Order and the unblinding

| resource | task | status |
|---|---|---|
| CPU | MXfold2 retraining | start now |
| GPU | train-split pseudo-perplexity; one-hot heads | start now |
| CPU | E-C2 fit and freeze | after the train perplexities exist |
| — | E-C1b decomposition and Δρ on withheld | after the one-hot heads |
| — | E-C2 on withheld | after the fit is frozen |
| — | Amendment 14: final test_novel analysis list | after all of the above |

Amendment 14 lists: all baselines; the greedy head; A-energy-fixed; E-C1 and
E-C1b confirmation; G1, plus G2 as secondary; and the H2.4 subgroup. Then
unblind once, hash the outputs, and record any rerun in `corrections.md`.

### Next report

1. Amendment 13 hash.
2. E-C1b: the decomposition and the Δρ verdict.
3. E-C2 frozen parameters and the withheld verdict, plus A-energy-fixed vs
   RNAfold and EternaFold.
4. MXfold2 training status.
5. Paper A draft progress: which sections are drafted.

Lead with E-C1b; it decides what C is allowed to be called.
