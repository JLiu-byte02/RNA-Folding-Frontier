# Research State

Updated: 2026-11-11

## Phase

`bootstrap` — infrastructure, literature grounding, and first direction selection.

## Established facts

- GitHub repository created: https://github.com/JLiu-byte02/RNA-Folding-Frontier
- Slurm partition `g-3090-2`: 4 nodes x 2 RTX 3090 GPUs (8 total, shared with
  other users). Project cap: 6 GPUs.
- Remote project root: `/home/junxuanliu/projects/RNA-Folding-Frontier` (empty).
- Login-node Python is 3.6.8; Conda modules are available. A dedicated environment
  must be created before any experiment.
- E0 first pass: ArchiveII overlaps bpRNA-1m TR0 by 20.47% exact sequence
  duplicates and 38.25% near-duplicates (k=15 containment >= 0.5). TS0, VL0 and
  bpRNA-new are sequence-clean at this threshold.

## Open questions

- Which of the PI's three candidate directions is most novel and feasible?
- Which datasets and splits give a defensible, leakage-audited benchmark?
- What is the strongest reproducible baseline achievable under the compute cap?

## Current blockers

- None.

## Do not claim

- Any performance number, SOTA statement, or novelty claim that is not backed by a
  committed artifact in this repository.

