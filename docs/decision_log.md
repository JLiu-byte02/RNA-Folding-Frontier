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

