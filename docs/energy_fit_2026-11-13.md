# (s, b) energy calibration — 2026-11-13

Fitted under Amendment 7 §A7.4 on the frozen train split, with v3 labels,
maximum-nested-subset truth, at least four true pairs, and 5-fold
clan-grouped CV. Objective: family-macro F1 as the mean of per-family mean
per-sequence F1.

## Coarse grid

The original 5×5 grid selected `s=0.25`, `b=0`, the lower edge of the s axis.
That triggered the pre-registered refinement rule.

## Refinement

Grid: `s ∈ {0.125, 0.1875, 0.25, 0.3125, 0.5}`,
`b ∈ {−3, −1.5, 0, 1.5, 3}`.

| s | b | mean fold F1 | sd |
|---:|---:|---:|---:|
| **0.1875** | **0** | **.6059** | .0065 |
| 0.25 | 0 | .6053 | .0050 |
| 0.3125 | 0 | .6051 | .0056 |
| 0.125 | 0 | .6047 | .0071 |
| 0.5 | 0 | .6015 | .0054 |

Selected and frozen: **s = 0.1875, b = 0**, fit sha256
`be54ddab6552102b645a42dff5ac60eaa92072680db062ae437c8c4889ad5ad5`,
16,558 training sequences. The selected point is interior on both axes.

Artifact: `experiments/002_e1prime_runs/energy_fit_refined_p2fix.json`;
hash recorded in `state/amendments.lock.json`.
