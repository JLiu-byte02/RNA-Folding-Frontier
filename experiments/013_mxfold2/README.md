# MXfold2 retrain (A13.1 secondary control)

The published MXfold2 weights were trained on a corpus that overlaps our test
families, so using them would mix training data. This directory holds a MixC
model retrained on **our** frozen train split under the MXfold2 authors'
README defaults, together with the manifests needed to reproduce the run
without shipping the 200 KB checkpoint.

## Training data

| quantity | value |
|---|---:|
| train rows | 17,746 |
| exported | 17,155 |
| dropped, longer than 600 nt | 577 |
| dropped, fewer than one pair | 14 |
| rows with crossing pairs (nested projection) | 1,846 |
| pairs before / after projection | 454,222 / 439,299 |

- split `splits/s1prime_filtered/train.jsonl`, sha256 `b8af2335...f94ecd`
- exported list sha256 `d2ba5c96...b5291c`
- labels: the frozen split `pairs` field (the same labels every head was
  trained on), projected to the maximum nested subset because MXfold2's folding
  layer cannot represent crossing pairs. Note that the export manifest's
  `label_policy` string says "vim3"; that string was wrong and is corrected in
  `docs/corrections.md`, the exported structures are the intended ones.

## Run

| field | value |
|---|---|
| Slurm job | 15631 `rff_mx_train` |
| partition / GPU | `g-3090-2`, one RTX 3090 |
| model | MixC (Turner-integrated), MXfold2 commit `51b2136` |
| epochs / batch | 10 / 1 |
| max length | 600 nt |
| wall time | 29,606 s (8 h 13 m) |
| state | COMPLETED, exit 0 |
| param sha256 | `54b93dcb...2881dd` |
| config sha256 | `446d79c0...f0bd7e` |

`model.pth` stays on the server; `train_manifest.json` and `model.conf` are
committed so the checkpoint can be re-derived and audited.

## Folding evaluation splits

`scripts/predict_mxfold2.py` wraps the upstream CLI: chunked, resumable via a
`.partial.jsonl` cache, explicit 600 nt cap with the excluded ids written next
to the output, and a strict parse that fails rather than silently dropping a
row. `scripts/jobs/mxfold2_fold.sbatch` folds `test_withheld` and `id_test`
with 16 CPU workers (the DP dominates, so CPU beats one GPU here).

Outputs land in `mxfold2_<split>.jsonl` here and are then fed to
`scripts/jobs/rebuild_cohort_4tools.sbatch`, which builds the four-method
common cohort and every downstream table.
