# Amendment 11 — retrain the P3 condition-a arm on the common split

Status: written 2026-11-12 after the P3 b/e heads finished, but **before any
D1-core score was computed**. It is a procedural correction to Amendment 10
§A10.5, not a new hypothesis.

## A11.1 — The mismatch

Amendment 10 said the condition-a arm could reuse `runs/p2fix_seed{17,23,41}`
because it used the same cache and hyperparameters. The cache ID sets are in
fact identical (`giga_a_train`, `giga_block_b_train`, `giga_block_e_train` all
cover 17,260 records; pairwise differences are empty), but the recorded
training/validation split counts are not:

| arm | fit | val | total |
|---|---:|---:|---:|
| p2fix a (seed 17) | 16,774 | 486 | 17,260 |
| P3 b (seed 17) | 16,716 | 544 | 17,260 |
| P3 e (seed 17) | 16,716 | 544 | 17,260 |

The p2fix result was produced before the current `train_head_s1.py`/split
snapshot was fixed; its validation family selection is not reproducible from
the current pair of files. Reusing it would make the a arm's early-stopping
data differ from both block arms.

## A11.2 — Corrected P3 a arm

Retrain condition a with the current `train_head_s1.py`, the same
`splits/s1prime_filtered/train.jsonl`, and `data/cache/giga_a_train`, using the
same hyperparameters and seeds as the block arms:

```
runs/p3_a_seed{17,23,41}/best.pt
```

The retraining command is `scripts/jobs/p3_train_a.sbatch <seed>`. The only
changed input relative to p2fix is the validation-family split, which makes it
identical to the b/e arms; no evaluation split is used.

All three P3 arms are then evaluated on the intersection of records covered by
the three caches, and the D1-core rule remains exactly Amendment 10 §A10.5:
at least one block arm must have a seed-averaged withheld all-pairs
family-macro ΔF1 ≥ .02 with a 95% family-bootstrap CI lower bound above zero.

## A11.3 — What does not change

- P2-block test-time H3 still uses the p2fix heads; that analysis is a
  test-time intervention and was separately pre-registered.
- No threshold, seed, cache geometry, population rule, or D1-core threshold is
  changed.
- The P3 block arms already trained are not retrained; only the a arm moves to
  the common split.
- D1-core is not scored until this amendment is committed and its hash is in
  `state/amendments.lock.json`.
