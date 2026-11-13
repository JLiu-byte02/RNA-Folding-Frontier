# Paper A — working draft (2026-11-13)

**Status.** Sections 1–4 and Tables T1–T6 are drafted from committed
artifacts. The title, abstract, discussion, and both conclusion wordings are
deliberately **held** pending (i) the frozen four-tool common cohort, (ii) the
MXfold2 retrain, and (iii) the Amendment 14 one-shot `test_novel` unblinding.
Nothing in this file is a submitted claim.

**Standing rule.** Every number below is traceable to a committed artifact and
was produced after the `pair_prf` zero-TP correction (Amendment 9). Numbers that
predate the correction appear only in `docs/corrections.md` and are never cited
here.

---

## 1. Introduction (draft)

RNA secondary-structure prediction is routinely benchmarked on ArchiveII and on
held-out Rfam families, and recent work reports large gains from RNA foundation
models (FMs) fine-tuned with a shallow pair head. Two properties of that
evaluation are rarely quantified. First, the benchmark sequences are drawn from
the same corpora that train both the FMs and the classical baselines, so a
fraction of the reported accuracy can be sequence memorisation rather than
folding. Second, the *labels* for Rfam-derived evaluation sets are produced by
projecting a consensus secondary structure from alignment coordinates onto an
ungapped sequence; when that projection is done at the string level, pairs are
silently lost or invented.

We report an audit of both layers, followed by three pre-registered mechanism
experiments that ask whether the residual FM advantage is explained by a
head-level shortcut, by evolutionarily consistent distal context, or by
energy-consistent decoding. All three are negative or unsupported under their
pre-registered bars. We report them as such, together with the one mechanism
result that did survive: pretraining corpus coverage predicts the FM's
per-family advantage over a thermodynamic baseline.

The contribution is therefore an audit plus well-powered negative mechanism
results, not a state-of-the-art claim.

## 2. Results

### 2.1 Benchmark overlap is large and clade-structured (F1, T1)

ArchiveII (3,370 sequences in our mirror) overlaps the union of two standard
trainers, RNAStrAlign and bpRNA TR0, at **61.20% exact** sequence duplicates and
**88.40% near-duplicates** (k = 15 containment >= 0.5). Against bpRNA-1m TR0
alone the exact rate is 20.47% and the near-duplicate rate 38.25%. ArchiveII
also contains 313 internal duplicate sequences. The apparent cleanliness of the
long families (16S, 23S, group II) is an artefact of the trainers capping length
at 600 nt, not of genuine separation.

Overlap alone would not inflate F1; overlap *plus label identity* does. Labels
agree at mean pair F1 **0.9937** between ArchiveII and RNAStrAlign and **0.9195**
between ArchiveII and bpRNA TR0 (3,370-record mirror; telomerase 0.6606 is the
outlier family). A duplicated sequence therefore arrives with an almost
identical target.

Clan structure compounds this. Rfam 15.0 contains 4,178 families in 147 clans
with 463 clan memberships; ArchiveII spans only **8 clans**. tRNA and tmRNA share
CL00001, so the standard leave-one-family-out protocol is not clan-disjoint on
that fold and cannot separate "generalises to a new family" from "generalises to
a sibling of a training family".

### 2.2 String-level label projection corrupts the target (F2, T2)

Rfam seed alignments store the consensus structure `SS_cons` in alignment
columns. Projecting it character-by-character onto the ungapped sequence and
then parsing the resulting string is the obvious implementation. It is wrong in
both directions: a gap that removes one helix endpoint lets the parser pair the
surviving endpoint with a different closing column (**invented** pair), and a
disturbed stack can lose a genuine pair.

On the 92,781 seed rows in our splits:

| quantity | value |
|---|---:|
| rows whose string-level (v1) projection does not parse | 16,540 (17.8%) |
| of those, repaired by v2 (per-character, gap-tolerant) | 7,670 |
| of those, still failing under v2 | 9,322 |
| rows parsed by pair-level v3 | 92,781 (100.00%) |
| pairs, v1 effective / v2 effective / v3 | 1,669,314 / 2,029,886 / 2,517,623 |
| rows with <4 v3 pairs | 5,372 (5.79%) |

On the 83,459 rows where both v2 and v3 parse, the pair sets still differ on
**6,757 rows**: 27,367 v2-only pairs and 58,637 v3-only pairs. "Both parse" is
therefore not "both correct"; v3 — parse `SS_cons` in column coordinates, map
both endpoints, drop the pair if either column is a gap, then apply the frozen
canonical/GU filter — is the specification.

The consequence for evaluation is not cosmetic. In the withheld split, 11,527 of
17,371 rows (66%) carry an Rfam seed label. Of those, **2,392 (20.75%) had an
empty v1 target**, which a scorer can read as "this RNA has no structure" rather
than "the projection failed". Under v3 only 112 are genuinely empty and 341
(2.96%) fall below the four-pair threshold. Any benchmark that scores a
string-projected Rfam evaluation set therefore contains a non-trivial block of
rows whose target is a bug.

The same issue exists on the bpRNA side: 15,794 of 66,715 bpRNA-1m sequences
(23.67%) contain at least one crossing pair, 6.18% of all pairs. A decoder or
loss that silently reduces to a nested set is optimising a different target, and
must say so.

### 2.3 Head-level context manipulation is null (F3, T3)

If the FM advantage were a family shortcut implemented in the head, then
restricting the head's input to a local window during *training* should not hurt
out-of-distribution accuracy. We tested this at matched everything: same
sequences, same backbone, same fit/validation split (16,716/544), same decoder,
same v3 nested truth, three seeds, on a common population of 16,163 withheld
sequences in 409 families.

| arm | family-macro F1 | 95% CI | micro F1 |
|---|---:|---|---:|
| a native | .5959 | [.5726, .6191] | .5968 |
| e e-block (approx. k=32 crop) | .5965 | [.5733, .6196] | .5852 |
| b b-block (dinucleotide-shuffled distal) | .5860 | [.5630, .6088] | .5766 |

Paired family-level deltas versus a: e **+.0007 [−.0079, +.0093]**, b
**−.0099 [−.0159, −.0036]**. The pre-registered falsification rule required
either block arm to beat a by >= .02 with a CI excluding zero; neither does.
**D1-core is falsified.** On id_test the ordering is weaker but does not
reverse (e −.0095 [−.0249, +.0058], b −.0208 [−.0359, −.0068]).

This rules out one specific account — that distal context is exploited *by the
head* in a way that transfers to unseen families. It does not rule out
pretraining contamination, and it does not show that full context is useless:
the b arm is worse, so distal context carries real information. It is a null on
the shortcut hypothesis, not a positive claim about locality.

### 2.4 Homologous test-time context transfers only partially (F4, T4)

Test-time substitution of the distal region gives a descriptive ladder on the
same 16,163-sequence withheld population: native a .5911, dinucleotide-shuffled
b .5629, same-family homolog c .6066, other-clan d .5543, N-mask e .5523.

The pre-registered H3 ratio rule was `R_c = (c-e)/(a-e) >= .90` and
`R_d <= .30` with CIs excluding the bar. Observed on the common a/c/e population
(8,476 sequences, 93 families): `R_c = .895 [.726, 1.052]`, `R_d = .225
[.019, .342]`. **H3 is not supported** — both intervals cross their bars. The
point estimates are directionally consistent with homologous context
transferring better than unrelated context, and we report them as descriptive.

An earlier P2-central reading of `R_c ~ 102%` was a population mismatch (the
conditions were scored on different sequence sets); on the common population it
is `.747 [.564, .960]` on kept–kept pairs. The corrected value, not the earlier
one, is the citable number.

### 2.5 Energy-consistent decoding does not remove long-range false positives (F5, T5)

Long-range false positives are the failure mode that motivates structure-aware
decoding. On the withheld split (p2fix seed 17) there are **4,391** long-range
FP helices; their mean isolated energy is **−1.43 kcal/mol** versus −1.23 for
true helices, and **25.6%** lie entirely inside the RNAfold MFE (35.7% share at
least one pair). They are therefore not free-energy artefacts that a physical
prior should erase.

We fitted a single global scale `s` and offset `b` on the train split under a
pre-registered refinement rule, froze `(s, b) = (.1875, 0)` with hash
`be54ddab...`, and decoded with an exact non-crossing DP at a pair budget
matched to greedy.

| split | greedy TP | DP TP | greedy LR-FP | DP LR-FP | verdict |
|---|---:|---:|---:|---:|---|
| id_test | 28,965 | 29,835 | 3,284 | 4,205 | not confirmed |
| withheld | 284,906 | 306,046 | 17,630 | 20,218 | not confirmed |

DP gains true positives at the matched budget on both splits, but it gains
long-range false positives at the same time. The pre-registered criterion
required fewer long-range FPs and at least as many TPs; **H2.4 is not
confirmed**. Per length band the trade is favourable below 300 nt and dominated
by the >=300-nt bin, where a 2,508-TP gain costs 2,928 extra long-range FPs.

A stronger prior was also tested and failed: a coverage-gated context term
collapsed to its grid boundary (`s_max = .1875`, `alpha = beta = gamma = 0`), so
the fitted gate is a pure constant scale that is **bit-identical** to its own
control and reverses sign on id_test (Δ `−.0064 [−.0122, −.0002]`). We report
that as a failed gate, not as a mechanism. The coarse grid restricted
`beta >= 0`, so negative dependence remains untested rather than excluded.

### 2.6 Pretraining coverage predicts the FM advantage (F6)

One mechanism result survived. On the 409-family common withheld population the
FM's per-family advantage over RNAfold correlates with two independent proxies
for how much of the family sits in the pretraining corpus:

| proxy | partial rho | 95% CI |
|---|---:|---|
| Rfam full-region membership | **.338** | [.244, .429] |
| RiNALMo pseudo-perplexity | **−.242** | [−.332, −.147] |

The pre-registered rule (P1 rho >= .30 with CI excluding 0; P2 rho <= −.20) is
met. The **specificity control did not pass**: the same slope for a one-hot
head is .253, so `d_rho = .085 [.009, .160]`, below the pre-registered .15. We
therefore report coverage as a *shared* predictor of accuracy on this benchmark,
not as something specific to the pretrained representation. The mean advantage
is also small (+.0051 family-macro across families, range −.742 to +.857), so
this is a mechanism result and not a performance result.

## 3. Methods

### 3.1 Data and splits

- **Corpora.** bpRNA-1m (TR0/VL0/TS0 as published), ArchiveII, RNAStrAlign, Rfam
  15.0 seed alignments and covariance models, and the HuggingFace
  `rouskinlab/bpRNA-1m` mirror (66,715 sequences) used as the label superset.
- **Family assignment.** Rfam families by accession; bpRNA sequences by cmscan
  against the pressed Rfam 15.0 CMs with `--cut_ga`.
- **Splits (`s1prime`, frozen).** 3,861 clans: 1,983 novel / 376 withheld / 1,502
  retained. `train` 17,746; `id_test` 1,402; `test_withheld` 17,371;
  `test_novel` 42,873. `test_novel` is embargoed and has never been scored.
- **Common populations only.** Every comparative number is computed on the
  intersection of the arms being compared, on identical ids, identical truth,
  and identical decoding. Per-tool populations are never compared side by side.

### 3.2 Label construction (v1 / v2 / v3)

`SS_cons` is parsed in alignment-column coordinates with all bracket types plus
letter pairs. Each column pair is mapped to sequence coordinates; if either
column is a gap in that sequence the pair is dropped; the surviving set is
filtered to canonical and GU pairs. This is `labels_v3` (92,781 rows, 2,517,623
pairs, sha256 `b8fe71ab...`), frozen by Amendment 8 as the primary *and*
secondary truth. v1 (string projection, parse-or-empty) and v2 (per-character
gap-tolerant projection) are provenance only.

Eligibility is decided at scoring time, not at split-construction time:
sequences with fewer than four true pairs under v3 are excluded from primary F1
and their predicted-pair counts are reported separately as a false-positive
diagnostic. This replaces an earlier dataset-level ">= 4 pairs in `SS_cons`"
condition that was applied when the split was built.

### 3.3 Scoring convention, including a disclosed bug

Pairs are scored as a set. `tp + fp + fn == 0` (both empty) yields F1 = 1.0;
any other case with `tp == 0` yields F1 = 0.0. An earlier implementation
returned 1.0 whenever `precision + recall == 0`, i.e. it scored *complete
disagreement* as perfect. The bug was found while auditing a template control,
is fixed in `src/rff/analysis/labels.py`, and has a regression test
(`test_pair_prf_zero_overlap_is_not_vacuous`). Every number in this paper was
produced after the fix; the pre-fix numbers and their corrected replacements are
listed in `docs/corrections.md`. The visible consequence was large — a
dinucleotide-shuffled control "scored" .90 before the fix and .017 after it.

Family-macro F1 is the unweighted mean of per-family mean F1; micro F1 pools
pairs across sequences. CIs are 95% family-bootstrap intervals (10,000
resamples over families, seed-fixed). Paired deltas are computed on the same
resampled families.

### 3.4 Models and training

Backbone: RiNALMo-giga (frozen). Head: a shallow pair-prediction head over the
outer-product of residue embeddings with a symmetric outer-sum, trained with
BCE on v1 pairs for the Round 4 experiments and on v3 pairs thereafter. Decoder:
threshold + one-partner + non-crossing greedy; a bare-threshold arm and a nested
DP arm are reported where they isolate a constraint. No backbone is fine-tuned;
no model is pretrained from scratch in this work.

Training runs are recorded with config, log, checkpoint, and evaluation JSON
under `experiments/` and `runs/`. Checkpoints produced under the buggy selection
metric (`runs/p2_seed*`) are quarantined and not used.

### 3.5 Context interventions

- **P2-block (test time).** Replace a 16-nt block plus 32-nt flanks in the distal
  region with: (a) native, (b) dinucleotide-shuffled block, (c) same-family
  homolog block, (d) other-clan block, (e) N-run.
- **P3 (train *and* test).** Train the head with only the local window
  (a vs e-block vs b-block) and evaluate on the full input. This separates "the
  head uses distal context" from "the distal tokens change the input
  distribution".

### 3.6 Decoders and energy calibration

Greedy: sort cells by probability, accept if both endpoints are free and the
pair does not cross an accepted pair. DP: exact maximum-weight non-crossing
matching, linear in memory and quadratic in time per sequence, cross-checked
against exhaustive enumeration on all components with <= 16 arcs (33,714
components, 0 mismatches; the old greedy reduction was suboptimal on 87
components and is kept as a regression test).

Energy decoding adds `s * (stacking energy) + b` to each cell's log-odds, with
`s, b` fitted on train by a pre-registered coarse-then-refined grid and frozen
before any evaluation. `(s, b) = (.1875, 0)` was selected at the refined grid's
interior and hashed before H2.4 was run.

### 3.7 Statistics

Family bootstrap for CIs and paired deltas; family-level partial correlation
(controlling for length) for E-C1; pre-registered support rules with explicit
bars. No test-set-driven model selection: every choice that could move a
reported number was fixed by an amendment before the corresponding number was
computed.

### 3.8 Pre-registration and amendment history

Fourteen amendments (`experiments/001_e1prime/AMENDMENT_*.md`) record the
protocol, the hypothesis, the success bar, and the analysis list *before* the
corresponding evaluation. Amendments are append-only; corrections are recorded
in `docs/corrections.md` with the affected numbers and their replacements.
Frozen artifact hashes (labels, splits, eval subsets, `(s, b)`) are recorded in
`state/amendments.lock.json` as hashes of git blob bytes, and a guard refuses to
score a dirty worktree.

`test_novel` has a one-shot unblinding list fixed by Amendment 14 and has not
been invoked.

### 3.9 Reproducibility

The repository is the record. A CI job clones it into a fresh directory,
installs the package, runs the test suite, and verifies every hash in
`state/amendments.lock.json` against the corresponding git blob. Data
provenance, download commands, and per-artifact hashes are in
`docs/data_sources.md` and `docs/remote_provenance_2026-11-12.md`.

## 4. Limitations

1. The audit covers ArchiveII, bpRNA, RNAStrAlign, and Rfam 15.0. We do not
   claim the same overlap rates for other benchmarks, and we have not verified
   which published papers used string-level projection — only that the failure
   mode exists and is large in this corpus.
2. D1-core is a null on one intervention (local-window training under
   RiNALMo-giga and this head). It does not exclude contamination, and it does
   not exclude shortcut exploitation by a different architecture or a
   fine-tuned backbone.
3. The H3 ratios are descriptive; the pre-registered bar was not met.
4. The energy decoder uses a single global scale. Family-specific or
   length-specific calibration was not tested.
5. `beta >= 0` was imposed on the gate grid, so negative dependence is untested.
6. Our head does not beat RNAfold or EternaFold on family-macro F1 on the
   current per-tool populations. The frozen cohort has not yet been scored, so
   no cross-tool statement is made.
7. RiNALMo-giga is AGPL-3.0. The repository licence is undecided, which blocks
   any public release of weights.

## C. Conclusion — HELD

Two candidate wordings are drafted; neither is final.

**C1 (audit-led).** RNA secondary-structure benchmarks contain two layers of
apparent performance that are not folding: heavy sequence overlap with the
training corpora, and labels produced by a projection that fails on one row in
five. After correcting both, the advantage of a frozen RNA foundation model over
a thermodynamic baseline is small, correlates with pretraining coverage, and is
not attributable to a head-level family shortcut, to evolutionarily consistent
distal context, or to energy-consistent decoding.

**C2 (negative-led).** Three plausible mechanisms for the residual advantage of
RNA foundation models on secondary-structure prediction — head-level shortcuts,
homologous distal context, and energy-consistent decoding — fail their
pre-registered tests under matched conditions, while corpus coverage predicts
per-family accuracy on two independent proxies but not specifically for the
pretrained representation. We argue the field should report label-projection
and overlap audits alongside F1.

## Claim table (T6)

| # | claim | status | strength |
|---|---|---|---|
| 1 | ArchiveII is 61% exact / 88% near-duplicate against common trainers | established | strong |
| 2 | ArchiveII spans 8 of 147 Rfam clans; leave-one-family-out is not clan-disjoint | established | strong |
| 3 | String-level seed projection fails on 17.8% of rows and corrupts 8.1% of the rest | established | strong |
| 4 | Head-level distal-context manipulation does not improve OOD F1 (D1-core) | falsified hypothesis | strong (pre-registered) |
| 5 | Homolog test-time context transfers more than other-clan context, below the pre-registered bar (H3) | descriptive | weak |
| 6 | Energy-consistent decoding adds TPs but also long-range FPs (H2.4) | not confirmed | strong (pre-registered) |
| 7 | Corpus coverage predicts FM advantage; specificity not shown (E-C1 / E-C1b) | supported / failed control | moderate |
| 8 | Any SOTA or "first" claim | **not available** | — |
| 9 | Any cross-tool comparison | **blocked** until the frozen cohort | — |
