# P2-central retention ratios with family-bootstrap CIs

Computed by `scripts/bootstrap_p2central.py` from the corrected p2fix
dumps; fixed threshold 0.5. The ratio is `(X - e) / (a - e)` on
family-macro F1; CIs are 10,000 family-bootstrap resamples.

| split | readout | ratio | n | families | point | 95% CI | undefined |
|---|---|---|---:|---:|---:|---|---:|
| id_test | kept_kept | R_c | 449 | 25 | 0.6122 | [0.4128, 0.8530] | 0.0000 |
| id_test | kept_kept | R_d | 973 | 128 | -0.1729 | [-0.6170, 0.0798] | 0.0000 |
| id_test | all_pairs | R_c | 526 | 30 | 0.4950 | [0.3762, 0.6108] | 0.0000 |
| id_test | all_pairs | R_d | 1287 | 155 | -0.1124 | [-0.1748, -0.0642] | 0.0000 |
| test_withheld | kept_kept | R_c | 7572 | 80 | 0.7471 | [0.5638, 0.9602] | 0.0000 |
| test_withheld | kept_kept | R_d | 13834 | 358 | 0.0820 | [-0.0399, 0.2030] | 0.0000 |
| test_withheld | all_pairs | R_c | 8476 | 93 | 0.5039 | [0.3900, 0.6057] | 0.0000 |
| test_withheld | all_pairs | R_d | 16163 | 409 | -0.0689 | [-0.1028, -0.0360] | 0.0000 |

## Condition means on the common population

These are the seed-averaged family-macro F1 values that enter each
ratio, so a reader can see the population restriction directly.

| split | readout | numerator | n | a | numerator | e |
|---|---|---|---:|---:|---:|---:|
| id_test | kept_kept | c | 449 | 0.5970 | 0.5616 | 0.5058 |
| id_test | kept_kept | d | 973 | 0.5759 | 0.5035 | 0.5142 |
| id_test | all_pairs | c | 526 | 0.5649 | 0.4380 | 0.3135 |
| id_test | all_pairs | d | 1287 | 0.5837 | 0.3224 | 0.3488 |
| test_withheld | kept_kept | c | 7572 | 0.6116 | 0.5927 | 0.5368 |
| test_withheld | kept_kept | d | 13834 | 0.5910 | 0.5231 | 0.5171 |
| test_withheld | all_pairs | c | 8476 | 0.6125 | 0.4927 | 0.3712 |
| test_withheld | all_pairs | d | 16163 | 0.5911 | 0.3463 | 0.3620 |

The ratio is undefined when `a - e <= 0`; the undefined fraction counts
bootstrap draws in that case. These CIs are diagnostic and do not
replace the P3 D1-core verdict.
