# Claude Round 0 — Direction Selection

- Date: 2026-11-11
- Model: Claude Opus 5.5, High effort
- Session: https://claude.ai/chat/cc2ad40a-406e-481e-8e7b-89a64b78f2f4
- Executor prompt: PI briefing with objective, constraints, and request for three
  directions plus first experiment batch.

## PI summary of the literature

- Cross-family failure is documented (Szikszai et al., Bioinformatics 2022).
  CHANRG (arXiv 2603.22330, March 2026, non-peer-reviewed per PI) reports FM mean
  F1 0.6731 in-distribution vs 0.1796 OOD, while structured decoders retained
  92.3%; failure analysis says FMs find local stems but fail global topology.
- bpRNA TS0 and RNAStrAlign random splits are intra-family and uninformative for
  generalization; ArchiveII is only meaningful leave-one-family-out.
- RiNALMo inter-family gains are uneven (telomerase F1 reported as 0.12).
- Three contamination channels: FM pretraining corpora, structRFM's 21M
  supervised pairs, and comparative labels that encode family consensus.
- Differentiable folding already exists (Matthies et al., NAR 2024; Krueger &
  Aviran 2025; EternaFold).

## Direction D1 (selected): Local evidence, global physics

**Hypothesis.** The transferable part of FM knowledge about pairing is local
(stacking and loop context within +/-k nt). Long-range context mostly encodes
family identity, a shortcut. Global wiring should come from exact combinatorial
inference, not a 2D CNN.

**Mechanism.**
1. Embed each position with a frozen FM on a local crop (cannot see global
   context).
2. Map crop embeddings to additive corrections on a nearest-neighbor energy
   model.
3. Decode with Zuker/McCaskill dynamic programming.

**Closest prior work.** MXfold2 (learned scores + Turner + global features),
EternaFold, Krueger & Aviran 2025 (global parameters, not context-dependent),
CHANRG (diagnosis without a fix). Claimed gap: nobody has measured context length
as a causal variable for the ID/OOD gap, nor fed crop-embedded FM features into a
nearest-neighbor DP. _[PI marked both as VERIFY.]_

**Minimal decisive experiment.** Sweep k on family-held-out splits. Prediction:
OOD F1 peaks at finite k* while ID F1 keeps rising with k.

**Cost.** ~200 GPU-h plus substantial CPU for DP.

**Failure modes.** Local features add nothing beyond Turner; DP too slow for
L > 500; pseudoknots out of scope.

## Other directions

- **D2: structure-free supervision from chemical probing** (latent-ensemble
  learning). Train energy corrections only through ensemble likelihood of
  SHAPE/DMS reactivity (Ribonanza). Novelty medium; ~500 GPU-h; weak
  identifiability because reactivity constrains pairing status, not partners.
  Kept as pivot.
- **D3: evidence-stratified ground truth** (R-scape/CaCoFold-supported pairs).
  Cheap (<30 GPU-h) but primarily a benchmark/audit contribution.

## Comparison

| | Originality | Feasibility | Evidence | Incremental risk |
|---|---|---|---|---|
| D1 | High if k-curve holds | Good | Mechanistic, falsifiable | Medium ("MXfold2 + FM") |
| D2 | Medium | Heavy | Strong, label-independent | Low leak, high identifiability risk |
| D3 | Medium | Excellent | Strong but descriptive | High |

**Pick: D1.** D3's stratified evaluation becomes standard protocol. D2 is the
planned pivot (same energy-correction model, different loss).

## First batch

Splits after leakage audit E0:

- S1: CHANRG GenC (clan) and GenF (family), if released.
- S2: ArchiveII leave-one-family-out (9 folds, Szikszai protocol).
- S3: bpRNA-new.
- S4: temporal PDB set (post-cutoff, <80% identity, no shared Rfam clan).

Training: bpRNA-1m filtered to remove every test clan (cmscan vs Rfam 15.0 +
MMseqs2).

Statistics: family is the unit of analysis; family-macro F1, family-level
bootstrap 95% CIs, paired Wilcoxon across families, 3 seeds.

**E0 leakage audit (prerequisite).** Max train-test identity, clan overlap,
test hits in FM pretraining corpora. Rebuild any split with >1% clan overlap.

**E1 context x decoder factorial (decisive).**
- Features: frozen RiNALMo + RNA-FM; crop half-width
  k in {8,16,32,64,128,full}.
- Decoders: (A) nearest-neighbor DP with learned corrections, built on MXfold2's
  DP; (B) matched-capacity 2D ResNet head.
- Baselines: RNAstructure, EternaFold, MXfold2, RiNALMo full-context fine-tune,
  CHANRG reported numbers.
- Metrics: base-pair F1, P, R, Stem-F1, OOD/ID retention, F1 by length bin.
- Budget ~12 cells x 3 seeds x ~4 GPU-h = 150 GPU-h.

**E2 shortcut probe.** Linear Rfam-clan classifier on per-position embeddings at
each k; prediction: clan decodability rises with k and tracks OOD drop.
Counterfactual transplant of hairpins/stems between families. ~20 GPU-h.

**E3 evidence/contamination-stratified evaluation.** Score E1 models on
R-scape-supported vs unsupported pairs, S4, and clean synthetic designs.

**Go/no-go after E1** (all three):
1. OOD F1 peaks at k* < full and beats full context by >= 0.05, CI excluding 0.
2. DP decoder beats ResNet head at k* by >= 0.05 OOD.
3. Best model beats EternaFold by >= 0.03 OOD on S1 and S2, retention >= 80%.

Pivot to D2 if (1)-(2) hold but (3) fails.

**Falsification.** Monotone OOD F1 in k falsifies the locality hypothesis;
DP <= ResNet at matched k falsifies the global-inference claim; flat clan
decodability in k falsifies the mechanism.

## PI requests before E0

1. Data: Rfam 15.0 CMs/seeds; bpRNA-1m(90), bpRNA-new, ArchiveII, RNAStrAlign;
   SPOT-RNA TS1-3 and PDB sets; CHANRG splits and license; Ribonanza (Kaggle
   terms).
2. Tools/weights: RiNALMo micro/mega/giga, RNA-FM; MXfold2 source and DP
   interface/license; Matthies/Krueger code; Infernal, MMseqs2, R-scape,
   RNAstructure, EternaFold.
3. Contamination: public pretraining sequence lists for RiNALMo and RNA-FM.
4. Literature verification: context-length vs OOD prior work; FM embeddings into
   MXfold2/nearest-neighbor DP; BPfold's exact mechanism; Krueger & Aviran scope;
   CHANRG peer-review status and reproducibility; CaCoFold venue.
5. Cluster: >=300 GB scratch, CPU cores per node, preemption policy.

## Standing evidence format

One directory per round `reports/roundNN/`:

- `SUMMARY.md` (<=300 words): what ran, deviations, pass/fail per pre-registered
  threshold, anomalies.
- `manifest.json`: git SHA, sha256 of split files, environment lock hash, Slurm
  job IDs, seeds, GPU-hours.
- Tables: T1 main results, T2 per-family F1, T3 leakage audit, T4 baseline
  reproduction (tolerance +/-0.02).
- Plots: k-curve ID vs OOD with CIs; per-family paired scatter; F1 by length bin;
  clan-probe accuracy vs k.
- Failure cases: five worst OOD examples with sequence and dot-bracket.
- Controls: label-shuffle run and Turner-only DP.

