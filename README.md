# RNA-Folding-Frontier

An independent research project toward a **novel, generalizable RNA secondary-structure
prediction method** with publishable scientific contribution.

This repository is the single source of truth for the project: code, experiment
configs, evaluation scripts, negative results, literature notes, decision records,
and paper drafts. Large datasets and model checkpoints live on the compute server.

> Status (2026-11-11): **bootstrapping**. No results yet. No performance claims are
> made in this repository until they are backed by reproducible artifacts.

## Research loop

The project runs as a long-horizon loop:

1. **PI (Claude Opus 5.5, High effort)** proposes directions, hypotheses, and
   experiment plans grounded in current literature.
2. **Agent (Codex)** validates feasibility, implements, runs experiments on the
   Slurm cluster, and reports real evidence back.
3. **PI** interprets evidence, kills or redirects weak lines, and defines the next
   round.

Every exchange and decision is persisted under `docs/` and `state/` so the project
can resume after context loss, browser failure, or rate limits.

## Repository layout

```
docs/            charter, living research state, decision log, protocols
state/           machine-readable state + next actions (resume point)
src/rff/         library code (data, models, training, evaluation)
configs/         experiment configurations
scripts/         local tooling (remote execution, sync, cycle management)
experiments/     one directory per experiment: config, logs, results, report
reports/         cross-experiment analyses and paper drafts
```

## Compute

- Slurm cluster, partition `g-3090-2` (NVIDIA RTX 3090), hard cap of **6 GPUs**
  across all of this project's jobs.
- Remote root: `/home/junxuanliu/projects/RNA-Folding-Frontier`.
- All jobs must record config, logs, checkpoints, and evaluation outputs.

## Reproduction

Reproduction instructions will be added as the first baseline is established.
Each experiment directory must contain a `README.md` with the exact command,
environment, data version, and expected outputs.

## Integrity policy

- No fabricated results, citations, or SOTA claims.
- Train/validation/test splits are fixed before modeling and never tuned on test.
- Sequence-level and structural homology leakage are audited and reported.
- Negative and null results are kept and published in this repository.

