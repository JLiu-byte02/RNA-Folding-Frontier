# Project Charter

## Objective

Develop an original RNA secondary-structure prediction method whose scientific
contribution and empirical evidence could support submission to a top venue
(Nature Communications, Nature Methods, Nucleic Acids Research, or equivalent
ML/comp-bio venue).

The target is **scientific novelty plus reliable, generalizable improvement**, not
module stacking or a single benchmark win.

## Scope

- Task: RNA secondary-structure prediction (base-pair probability / structure
  decoding), including generalization across families and datasets.
- Backbone strategy: prefer frozen or parameter-efficient adaptation of existing
  RNA foundation models (e.g., RiNALMo, RNA-FM). From-scratch large-scale
  pretraining is out of scope unless a compelling case is approved.
- Evaluation: bpRNA, ArchiveII, RNAStrAlign and successors as appropriate, with a
  protocol justified by current literature.

## Constraints

- At most 6 RTX 3090 GPUs in total for this project, scheduled via Slurm.
- Never disrupt other users' jobs.
- Prefer parallelizable, low-cost experiments; optimize the plan before asking for
  more compute.
- All training runs must persist config, logs, checkpoints, and evaluation output.

## Scientific standards

1. **Independent baselines.** Build and verify a credible baseline before claiming
   improvements.
2. **Leakage audit.** Check exact and near-duplicate sequences, family homology,
   and pretraining-corpus contamination risk. Report what cannot be ruled out.
3. **Fair comparison.** Same splits, preprocessing, decoding, and compute budget
   where possible.
4. **Ablations.** Every claimed component must have a removal experiment.
5. **Statistics.** Report uncertainty (bootstrap CIs, paired tests, multiple seeds).
6. **Generalization.** Cross-family and cross-dataset evaluation, including
   distribution shift.
7. **Mechanism.** Explain why the method works and when it fails.
8. **Test discipline.** The final test set is touched once per major claim.
9. **Honesty.** Failed hypotheses and negative results are recorded, not hidden.

## Roles

- **PI:** Claude Opus 5.5 (High effort) via claude.ai. Owns direction, hypotheses,
  success criteria, and interpretation.
- **Executor:** Codex agent. Owns feasibility review, implementation, experiment
  execution, verification, and evidence reporting.

The executor must challenge infeasible, unsafe, or methodologically invalid
instructions rather than execute them blindly.

