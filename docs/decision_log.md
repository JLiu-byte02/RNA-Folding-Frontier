# Decision Log

Append-only. Each entry: date, decision, evidence, alternatives, consequences.

## 2026-11-11 — Project bootstrap

- **Decision:** Start a clean, independent project (`RNA-Folding-Frontier`) rather
  than extending prior RNA_ThermoFusion work.
- **Rationale:** User requirement; avoids inheriting prior design choices and
  conclusions.
- **Consequences:** All baselines, splits, and claims must be rebuilt from scratch.

## 2026-11-11 — Compute envelope

- **Decision:** Use Slurm partition `g-3090-2`, max 6 GPUs concurrently, remote
  root `/home/junxuanliu/projects/RNA-Folding-Frontier`.
- **Rationale:** User-specified resource limit and server policy.
- **Consequences:** Favor frozen backbones, small adapters, and parallel short runs.

## 2026-11-11 — Round 0 direction: D1 selected

- **Decision:** Pursue **D1 — local evidence, global physics**: frozen FM
  embeddings on local crops -> additive corrections to a nearest-neighbor energy
  model -> Zuker/McCaskill DP. D3's evidence-stratified evaluation is adopted as
  standard protocol; D2 is the planned pivot.
- **Rationale:** Only D1 makes a causal, falsifiable claim about why FMs fail
  out of distribution (context length as a shortcut variable) and fixes it with a
  structural change rather than a new head.
- **Evidence status:** PI's literature claims are provisional and must be
  independently verified before E1 (see `docs/literature.md`).
- **Consequences:** E0 leakage audit is a hard prerequisite; no E1 number is
  trusted before T3 exists.

