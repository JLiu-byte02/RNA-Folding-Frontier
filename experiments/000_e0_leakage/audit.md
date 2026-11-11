# Leakage audit

Reference split: `train_tr0`

| split | n | len median | len max | families |
|---|---:|---:|---:|---:|
| train_tr0 | 10814 | 105 | 498 | 9595 |
| valid_vl0 | 1300 | 106 | 497 | 1147 |
| test_ts0 | 1305 | 109 | 499 | 1131 |
| archiveII | 3966 | 120 | 1800 | 10 |
| bprna_new | 5401 | 92 | 489 | 1 |

Cross-split exact duplicate sequences: **897**

## valid_vl0
- exact duplicates vs reference: 0 (0.00% of test sequences)
- near-duplicate query sequences (k=15): 4 (0.31%); 4 pairs
- shared families: 1 (0.1% of test families)
- shared clans: 0 (0.0% of test clans)

## test_ts0
- exact duplicates vs reference: 0 (0.00% of test sequences)
- near-duplicate query sequences (k=15): 6 (0.46%); 6 pairs
- shared families: 1 (0.1% of test families)
- shared clans: 0 (0.0% of test clans)

## archiveII
- exact duplicates vs reference: 812 (20.47% of test sequences)
- near-duplicate query sequences (k=15): 1517 (38.25%); 1622 pairs
- shared families: 0 (0.0% of test families)
- shared clans: 0 (0.0% of test clans)

## bprna_new
- exact duplicates vs reference: 0 (0.00% of test sequences)
- near-duplicate query sequences (k=15): 3 (0.06%); 3 pairs
- shared families: 0 (0.0% of test families)
- shared clans: 0 (0.0% of test clans)
