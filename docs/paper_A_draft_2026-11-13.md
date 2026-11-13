# Paper A — working draft (2026-11-13)

**Status.** Sections 1–4 and Tables T1–T7 are drafted from committed
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

### 2.7 No head advantage on a frozen common-cohort comparison (F7, T7)

Every method was scored on one frozen cohort of identical sequence ids under the
primary v3 labels, nested truth, and `--min-true-pairs 4`: 16,059 sequences in
409 families on withheld, 1,287 in 155 families on id_test. CIs are 10,000
family-block bootstrap resamples (seed 20261113); paired deltas resample
families, so the pairing respects family structure on both sides.

**Table T7.** Family-macro and micro pair F1 on the frozen common cohort.

| method | withheld macro F1 [95% CI] | withheld micro | id_test macro F1 [95% CI] | id_test micro |
|---|---:|---:|---:|---:|
| RNAfold | .6004 [.5744, .6263] | .5560 | .5955 [.5520, .6387] | .5295 |
| EternaFold | .6050 [.5792, .6301] | .5908 | .5916 [.5479, .6338] | .5465 |
| RNAstructure | .5901 [.5644, .6159] | .5482 | .5868 [.5430, .6302] | .5169 |
| head, seed 17 | .6054 [.5808, .6290] | .6040 | .6009 [.5594, .6407] | .6057 |
| head, seed 23 | .5785 [.5539, .6020] | .5897 | .5784 [.5359, .6194] | .5839 |
| head, seed 41 | .5854 [.5616, .6088] | .5933 | .5772 [.5366, .6172] | .5872 |
| head, 3-seed mean | .5898 [.5660, .6127] | .5955 | .5855 [.5453, .6243] | .5920 |

The head has the highest micro F1 on both splits (it is more aggressive and
recovers more base pairs overall), but its family-macro F1 is at or below the
thermodynamic tools. No paired family-block interval excludes zero:

| comparison | withheld delta [95% CI] | p | id_test delta [95% CI] | p |
|---|---:|---:|---:|---:|
| head_mean - RNAfold | -.0106 [-.0312, +.0095] | .3092 | -.0101 [-.0391, +.0180] | .4884 |
| head_mean - EternaFold | -.0152 [-.0350, +.0046] | .1256 | -.0061 [-.0311, +.0185] | .6218 |
| head_mean - RNAstructure | -.0003 [-.0204, +.0196] | .9840 | -.0013 [-.0296, +.0264] | .9332 |
| seed17 - RNAfold | +.0050 [-.0171, +.0272] | .6580 | +.0053 [-.0273, +.0375] | .7518 |
| seed17 - EternaFold | +.0004 [-.0214, +.0224] | .9726 | +.0092 [-.0229, +.0409] | .5708 |

The single seed that nominally exceeds RNAfold and EternaFold does so by
between .0004 and .0092 depending on the comparison, with every interval
spanning zero and 43-51% of families better in the four raw comparisons. We therefore
report **no advantage over the thermodynamic baselines**, not a tie and not a
loss: the cohort is too small in family count to resolve a .005 difference. The
cohort is also capped at 1,014 nt (withheld) and 990 nt (id_test), two families
are absent entirely, and MXfold2 had not finished training when this table was
frozen, so the four-method cohort will be rebuilt once it does.

## 3. Methods

### 3.1 Data and splits

- **Corpora.** bpRNA-1m (TR0/VL0/TS0 as published), ArchiveII, RNAStrAlign, Rfam
  15.0 seed alignments and covariance models, and the HuggingFace
  `rouskinlab/bpRNA-1m` mirror (66,715 sequences) used as the label superset.
- **Family assignment.** Rfam families by accession; bpRNA sequences by cmscan
  against the pressed Rfam 15.0 CMs with `--cut_ga --rfam --nohmmonly`.
- **Splits (`s1prime_filtered`, frozen).** 3,861 clans: 1,983 novel / 376
  withheld / 1,502 retained. After the reverse-homology and 15-mer screens:
  `train` 17,746; `id_test` 1,402; `test_withheld` 17,371; `test_novel` 42,873.
  The pre-screen `s1prime` counts are 17,748 and 1,409. `test_novel` is
  embargoed and has never been scored.
- **Common populations only.** Every comparative number is computed on the
  intersection of the arms being compared, on identical ids and identical
  truth. Decoding is held fixed within model-arm comparisons (greedy vs DP vs
  energy arms share one probability map); each external tool uses its native
  decoder. Per-tool populations are never compared side by side.

### 3.2 Label construction (v1 / v2 / v3)

`SS_cons` is parsed in alignment-column coordinates with all bracket types plus
letter pairs. Each column pair is mapped to sequence coordinates; if either
column is a gap in that sequence the pair is dropped; the surviving set is
filtered to canonical and GU pairs. This is `labels_v3` (92,781 rows, 2,517,623
pairs, sha256 `b8fe71ab...`), frozen by Amendment 8 as the primary *and*
secondary truth. v1 (string projection, parse-or-empty) and v2 (per-character
gap-tolerant projection) are provenance only.

Eligibility is decided at scoring time, not at split-construction time:
sequences with fewer than four true pairs under v3 are excluded from primary F1.
The D1 scorer and the external-baseline scorer also report their predicted-pair
counts separately as a false-positive diagnostic. The P2/P3 block scorer
(`src/rff/analysis/block_eval.py`) excludes the same rows without emitting that
diagnostic, so P2/P3 tables cannot be used to argue about empty-truth rows.
This replaces an earlier dataset-level ">= 4 pairs in `SS_cons`" condition that
was applied when the split was built.

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

Backbone: RiNALMo-giga (frozen). Head: each residue embedding is projected to
32 dimensions, and cell (i, j) is represented by the concatenation
`[z_i, z_j, z_i * z_j]` plus a learned relative-distance embedding; eight
dilated 3x3 residual blocks then produce one logit per cell, symmetrised as
`0.5 * (L + L^T)`. Training is BCE on the frozen bpRNA-era `pairs` field of the
split for every head; v3 changes evaluation labels only (Amendment 10).
Decoder: threshold + one-partner + non-crossing greedy; a bare-threshold arm
and a nested DP arm are reported where they isolate a constraint. No backbone
is fine-tuned; no model is pretrained from scratch in this work.

Training runs are recorded with config, log, checkpoint, and evaluation JSON
under `experiments/` and `runs/`. Checkpoints produced under the buggy selection
metric (`runs/p2_seed*`) are excluded from every headline number; they were
retained only to regenerate diagnostic dumps that were then re-scored with the
fixed metric.

### 3.5 Context interventions

- **P2-block (test time).** Replace a 16-nt block plus 32-nt flanks in the distal
  region with: (a) native, (b) dinucleotide-shuffled block, (c) same-family
  homolog block, (d) other-clan block, (e) N-run.
- **P3 (train *and* test).** Three arms train under different distal
  conditions -- native (a), N-run (e-block) and dinucleotide-shuffled
  (b-block) -- and are evaluated on the matching condition cache as well as on
  the native input. This separates "the head uses distal context" from "the
  distal tokens change the input distribution".

### 3.6 Decoders and energy calibration

Greedy: sort cells by probability, accept if both endpoints are free and the
pair does not cross an accepted pair. Decoder A: exact maximum-weight
non-crossing matching under the nearest-neighbour decomposition, O(n^3) time
and O(n^2) memory per sequence. The 33,714-component exhaustive check belongs to
a different routine -- the polynomial maximum-nested-subset solver used to
project reference structures (`max_nested_subset_of_arcs`), which matched
exhaustive enumeration on every component with <= 16 arcs (0 mismatches) and
showed the old greedy reduction was suboptimal on 87 of them; that check is a
regression test for the label solver, not for decoder A.

Energy decoding passes each pair's log-odds to ViennaRNA as a per-pair
pseudo-energy `-(s * score_ij + b)`; Turner stacking remains inside ViennaRNA.
`s, b` were fitted on train by a pre-registered coarse-then-refined grid and
frozen before any evaluation. `(s, b) = (.1875, 0)` sits at the upper boundary
of the refined grid (`s_max`) and was hashed before H2.4 was run. The fit used
5-fold cross-validation with folds assigned by `crc32(record_id) mod 5`, not by
clan as Amendment 7 specified; this deviation is recorded in
`docs/corrections.md`. E-C2 failed independently of it (the fitted gate is a
constant scale, bit-identical to its control), so no reported conclusion rests
on the leakage record-level folds could introduce; a clan-grouped refit has not
been run.

### 3.7 Statistics

Family bootstrap for CIs and paired deltas; family-level partial correlation
for E-C1, adjusting for median length, GC content, pair density, multiloop
presence, and label source; pre-registered support rules with explicit bars. No test-set-driven model selection: every choice that could move a
reported number was fixed by an amendment before the corresponding number was
computed.

### 3.8 Pre-registration and amendment history

Fourteen numbered amendments record the protocol, the hypothesis, the success
bar, and the analysis list *before* the corresponding evaluation: Amendment 1
is a section of `experiments/001_e1prime/PREREGISTRATION.md`, Amendments 2-14
are standalone files. Amendments are append-only; corrections are recorded in
`docs/corrections.md` with the affected numbers and their replacements.
`state/amendments.lock.json` pins the amendment files, the `test_novel`
evaluation list, the energy fit, and the E-C2 fit as hashes of git blob bytes;
label and split hashes live in `experiments/011_labels_v3/README.md` and
`docs/remote_provenance_2026-11-12.md`. The embargo guard refuses a
confirmatory read when any locked file differs from `HEAD`; it does not reject
unrelated worktree changes.

`test_novel` has a one-shot unblinding list fixed by Amendment 14 and has not
been invoked.

### 3.9 Reproducibility

The repository is the record for code, protocols, and evaluation artifacts;
weights, raw logs, and the corpora themselves stay on the server (`.gitignore`).
A CI job clones the repository into a fresh directory, installs the test
dependencies, runs the test suite, and verifies every hash in
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
6. Our head does not beat RNAfold, EternaFold, or RNAstructure on family-macro
   F1 on a single frozen common cohort (withheld 16,059 sequences / 409
   families; id_test 1,287 / 155; sequences of at most 1,014 nt). No paired
   family-block bootstrap interval excludes zero. This is one cohort under one
   label policy and does not generalise to other benchmarks or to long RNAs.
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
| 9 | Cross-tool comparison on one frozen common cohort | available: no head advantage | moderate (single cohort, <=1,014 nt) |
