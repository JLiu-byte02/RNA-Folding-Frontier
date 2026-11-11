# Experiment Protocol

This protocol is binding for every experiment in this repository.

## 1. Splits and leakage

- Define train/validation/test splits **before** any model selection.
- Report exact-sequence duplicates across splits; remove or explicitly account for
  them.
- Run a near-duplicate / homology audit (e.g., CD-HIT-EST or k-mer Jaccard) and
  report thresholds.
- Document known contamination risk from backbone pretraining corpora; state what
  cannot be excluded.

## 2. Baselines

- Reproduce at least one strong published baseline under our splits.
- Keep preprocessing, maximum length, and decoding identical when comparing.
- Record the exact command and environment for every baseline number.

## 3. Metrics

- Primary: micro/macro F1 over base pairs, and per-family breakdown.
- Secondary: precision, recall, MCC; base-pair probability calibration if
  applicable.
- Report mean ± std across seeds where stochasticity exists.

## 4. Statistics

- Bootstrap confidence intervals over sequences (not base pairs).
- Paired tests (e.g., Wilcoxon or paired bootstrap) for method comparisons.
- Pre-register the primary metric and success threshold for each experiment.

## 5. Ablations

- One factor at a time; retrain where the component affects optimization.
- Report compute-matched controls.

## 6. Test discipline

- Validation drives all iteration.
- Test is evaluated only for pre-registered final claims; each evaluation is logged
  with a timestamp and reason.

## 7. Artifacts

Every experiment directory contains:

```
README.md        question, hypothesis, command, environment, data version
config.yaml      exact configuration
results.json     machine-readable metrics
report.md        interpretation, including negative findings
```

