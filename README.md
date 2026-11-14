# RNA-Folding-Frontier

An independent research project toward a **novel, generalizable RNA
secondary-structure prediction method** with a publishable scientific
contribution. The repository is the single source of truth for code, protocols,
experiments, negative results, decision records, and the paper draft; datasets
and checkpoints stay on the compute server.

> **Status 2026-11-14 — negative result, frozen evidence.** A frozen
> four-method comparison shows that our frozen-RiNALMo pair head does **not**
> beat the thermodynamic baselines. We are not claiming SOTA. The paper in
> progress is an audit-plus-mechanism study, and the next direction is an open
> PI decision.

## What is established

- **Benchmark overlap.** ArchiveII overlaps RNAStrAlign + bpRNA TR0 at 61.2%
  exact and 88.4% near-duplicate identity (3,966-record leakage corpus);
  Rfam-15.0 clan structure makes leave-one-family-out non-clan-disjoint.
- **Label projection bug.** String-level `SS_cons` projection fails to parse on
  17.8% of Rfam seed rows and disagrees with pair-level projection on 8.1% of
  rows where both parse. `labels_v3` (pair-level, gap-aware) is the frozen
  primary and secondary truth.
- **No head advantage on one frozen cohort.** 16,040 sequences / 408 families
  (withheld) and 1,272 / 153 (id_test); the MXfold2 600 nt cap binds, so the
  cohort stops at 597/572 nt. Family-macro F1:

  | method | withheld | id_test |
  |---|---:|---:|
  | RNAfold | .6003 | .5997 |
  | EternaFold | .6049 | .5952 |
  | RNAstructure | .5900 | .5908 |
  | MXfold2 (retrained on our split) | .5518 | .5491 |
  | head, 3-seed mean | .5905 | .5901 |

  Paired family bootstrap: head−RNAfold −.0098 [−.0301,+.0107] and −.0097
  [−.0384,+.0187]; no interval against a thermodynamic tool excludes zero.
  The retrained MXfold2 is a 10-epoch control, not the published system, so the
  head's margin over it is not a competitiveness claim.
- **Mechanism tests were mostly negative.** Head-level distal-context
  manipulation (D1-core) is falsified; homologous test-time context (H3) misses
  its pre-registered bar; an energy-consistent decoder adds long-range false
  positives as fast as true positives; a coverage-gated scale collapses to a
  constant. Pretraining coverage correlates with per-family advantage, but a
  no-pretraining control fails the specificity bar, so it is reported as a
  shared predictor, not a pretraining mechanism.
- **A scorer bug was found and corrected.** `pair_prf` scored complete
  disagreement as perfect; every affected number is re-derived in
  `docs/corrections.md`.

## Repository layout

```
docs/            charter, research state, decisions, corrections, paper draft
state/           machine-readable state and next actions (resume point)
src/rff/         library code (data, models, training, evaluation, decoders)
configs/         experiment configurations
scripts/         local tooling and cluster job scripts
experiments/     one directory per experiment: config, logs, results, report
tests/           unit, regression, paper-number, and repository-hygiene checks
```

## Reproduction

```bash
python -m pytest -q          # full suite, including the gates below
```

- `tests/test_paper_numbers.py` re-parses every table in the paper draft and
  fails if a value does not match its committed JSON.
- `tests/test_public_repo_hygiene.py` fails if a credential or private key ever
  enters the repository.
- The CI job clones the repository into a fresh directory, verifies every lock
  hash against git blob bytes, and runs the suite.

Cluster-side regeneration of the headline table:

```bash
sbatch scripts/jobs/mxfold2_fold.sbatch          # retrained MXfold2 inference
sbatch scripts/jobs/rebuild_cohort_4tools.sbatch # cohort, scores, paired bootstrap
python scripts/render_cohort_tables.py ...       # markdown from the paired JSON
```

## Open decisions (PI)

1. Whether to spend the one-shot `test_novel` unblinding (Amendment 14).
2. Paper A framing: audit-led, negative-led, or coverage-mechanism-plus-decoder.
3. Whether to refit the energy scale with the clan-grouped folds Amendment 7
   specified.
4. Whether any method contribution is still defensible from this evidence, or
   the next cycle should test a new hypothesis with a pre-declared kill rule.

## Integrity policy

- No fabricated results, citations, or SOTA claims.
- Splits are frozen before modeling; no test-set-driven selection.
- Sequence and structural homology leakage are audited and reported.
- Negative and null results stay in the repository.
- Credentials never enter the repository (`docs/public_repo_hygiene.md`).
