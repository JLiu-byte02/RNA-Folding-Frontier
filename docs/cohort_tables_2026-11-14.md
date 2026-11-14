### splits/s1prime_filtered/test_withheld.jsonl — 16,040 sequences / 408 families

| method | family-macro F1 [95% CI] | micro F1 |
|---|---:|---:|
| RNAfold | 0.6003 [0.5736, 0.6266] | 0.5570 |
| EternaFold | 0.6049 [0.5792, 0.6301] | 0.5915 |
| RNAstructure | 0.5900 [0.5634, 0.6165] | 0.5490 |
| MXfold2 (retrained) | 0.5518 [0.5258, 0.5772] | 0.5161 |
| head, seed 17 | 0.6061 [0.5818, 0.6303] | 0.6072 |
| head, seed 23 | 0.5790 [0.5548, 0.6028] | 0.5921 |
| head, seed 41 | 0.5863 [0.5624, 0.6100] | 0.5959 |
| head, 3-seed mean | 0.5905 [0.5671, 0.6139] | 0.5982 |

| comparison | delta [95% CI] | p | frac families better |
|---|---:|---:|---:|
| head, 3-seed mean - RNAfold | -0.0098 [-0.0301, +0.0107] | 0.3368 | 0.409 |
| head, 3-seed mean - EternaFold | -0.0144 [-0.0341, +0.0052] | 0.1478 | 0.387 |
| head, 3-seed mean - RNAstructure | +0.0005 [-0.0193, +0.0208] | 0.9680 | 0.431 |
| head, 3-seed mean - MXfold2 (retrained) | +0.0387 [+0.0225, +0.0549] | 0.0000 | 0.591 |
| head, seed 17 - RNAfold | +0.0058 [-0.0165, +0.0284] | 0.6120 | 0.444 |
| head, seed 17 - EternaFold | +0.0013 [-0.0202, +0.0232] | 0.9088 | 0.434 |

Bootstrap: 10,000 family resamples, seed 20261113; min true pairs 4; labels `experiments/011_labels_v3/labels_v3_pairs.json.gz`.

### splits/s1prime_filtered/id_test.jsonl — 1,272 sequences / 153 families

| method | family-macro F1 [95% CI] | micro F1 |
|---|---:|---:|
| RNAfold | 0.5997 [0.5559, 0.6427] | 0.5437 |
| EternaFold | 0.5952 [0.5522, 0.6374] | 0.5582 |
| RNAstructure | 0.5908 [0.5468, 0.6343] | 0.5321 |
| MXfold2 (retrained) | 0.5491 [0.5048, 0.5928] | 0.4979 |
| head, seed 17 | 0.6057 [0.5652, 0.6449] | 0.6198 |
| head, seed 23 | 0.5828 [0.5411, 0.6235] | 0.5982 |
| head, seed 41 | 0.5817 [0.5411, 0.6214] | 0.6007 |
| head, 3-seed mean | 0.5901 [0.5512, 0.6275] | 0.6059 |

| comparison | delta [95% CI] | p | frac families better |
|---|---:|---:|---:|
| head, 3-seed mean - RNAfold | -0.0097 [-0.0384, +0.0187] | 0.5152 | 0.444 |
| head, 3-seed mean - EternaFold | -0.0052 [-0.0304, +0.0211] | 0.6708 | 0.438 |
| head, 3-seed mean - RNAstructure | -0.0008 [-0.0289, +0.0272] | 0.9624 | 0.464 |
| head, 3-seed mean - MXfold2 (retrained) | +0.0409 [+0.0160, +0.0663] | 0.0008 | 0.556 |
| head, seed 17 - RNAfold | +0.0060 [-0.0271, +0.0387] | 0.7068 | 0.510 |
| head, seed 17 - EternaFold | +0.0104 [-0.0221, +0.0438] | 0.5230 | 0.516 |

Bootstrap: 10,000 family resamples, seed 20261113; min true pairs 4; labels `experiments/011_labels_v3/labels_v3_pairs.json.gz`.
