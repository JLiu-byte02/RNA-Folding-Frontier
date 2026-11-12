# P2-block and P3 train-and-test

Status: pre-registered in `experiments/001_e1prime/AMENDMENT_10.md` before any
block cache was read. This directory records the implementation and the
eventual tables; it contains no result yet.

## What is being tested

- **P2-block test-time arm.** One condition-a-trained head per seed, with the
  distal context replaced per 16-nt block. Conditions a–e on `id_test` and
  `test_withheld` at fixed threshold 0.5. H3 asks whether a same-family homolog
  (c) retains at least 90% of the a−e difference on withheld, while an
  other-clan donor (d) retains at most 30%.
- **P3 train-and-test.** Separate heads trained on the train split with
  full context (a), an e-block cache (approximately a k=32 crop), and a
  b-block cache. Each is evaluated on the same condition in withheld. D1-core
  is falsified unless at least one block-trained head beats the a-trained head
  by ≥0.02 family-macro all-pairs F1 with a 95% CI excluding 0.

## Commands

Build one test or train cache:

```bash
sbatch scripts/jobs/p2_block_build_one.sbatch test_withheld b
sbatch scripts/jobs/p2_block_build_one.sbatch train e
```

Dump a condition-a head on all five test conditions:

```bash
sbatch scripts/jobs/p2_block_dump.sbatch 17 test_withheld
```

Train and dump a P3 block arm:

```bash
sbatch scripts/jobs/p3_train_block.sbatch e 17
sbatch scripts/jobs/p3_dump_block.sbatch e 17
```

Score on CPU:

```bash
sbatch scripts/jobs/p2_block_score.sbatch test_withheld
sbatch scripts/jobs/p3_score.sbatch
```

## Artifacts

- `data/cache/giga_block_{condition}_{split}/manifest.json` — geometry,
  coverage, passes, and timing for each cache.
- `experiments/002_e1prime_runs/dump_p2block_seed*/` — P2-block probability
  dumps.
- `experiments/002_e1prime_runs/p2block_h3_*.json` and
  `docs/p2block_h3_*_2026-11-12.md` — H3 table.
- `runs/p3_{b,e}_seed*/best.pt` — P3 heads.
- `experiments/002_e1prime_runs/d1core_test_withheld.json` and
  `docs/d1core_2026-11-12.md` — D1-core table and verdict.

No `test_novel` path is read by any of these commands.
