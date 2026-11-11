# Leakage audit

Reference split: `rnastralign`

| split | n | len median | len max | families |
|---|---:|---:|---:|---:|
| rnastralign | 37052 | 122 | 1850 | 1 |
| archiveII | 3966 | 120 | 1800 | 10 |
| test_ts0 | 1305 | 109 | 499 | 1131 |

Cross-split exact duplicate sequences: **1776**

## archiveII
- exact duplicates vs reference: 2005 (50.55% of test sequences)
- near-duplicate query sequences (k=15): 3246 (81.85%); 10401 pairs
- shared families: 0 (0.0% of test families)
- shared clans: 0 (0.0% of test clans)

## test_ts0
- exact duplicates vs reference: 66 (5.06% of test sequences)
- near-duplicate query sequences (k=15): 149 (11.42%); 342 pairs
- shared families: 0 (0.0% of test families)
- shared clans: 0 (0.0% of test clans)
