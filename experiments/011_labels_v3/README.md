# labels_v3 — pair-level seed projection

This directory holds the label correction frozen by `AMENDMENT_8.md`: Rfam seed
primary **and** secondary labels are `labels_v3`; v1 and v2 are provenance only.

## Why v3

v2 projects `SS_cons` character by character and parses the projected string. If
a gap removes one helix endpoint, the projected string can pair the surviving
endpoint with a different closing column and invent a base pair. v3 parses
`SS_cons` in alignment-column coordinates, maps both endpoints to the ungapped
sequence, drops the pair if either column is a gap, and then applies the frozen
canonical/GU filter.

## Artifacts

- `labels_v3_pairs.json.gz` — full replacement map, key
  `<RF accession>|<seed name>`, value a canonical pair list or `null`.
- sha256 `b8fe71abfab30c97298bc90140bd82fc7b787f2cc439b838d86a1b0b6b895ec1`
- 4,236,726 compressed bytes; 25,233,881 uncompressed bytes.
- `labels_v3_summary.json` — parse, agreement, sparse-label, eval-subset, and
  server-side per-split counts.

The artifact is byte-reproducible on Windows and Linux. `gzip.compress(...,
mtime=0)` still lets zlib choose the gzip header OS byte (`0x0a` on this machine,
`0xff` on the server), so the compressed stream was identical while the file
hash was not. `rff.data.io.gzip_deterministic` pins that metadata byte to
`0x0a`; the server reproduces the frozen hash exactly.

## Frozen numbers

| quantity | value |
|---|---:|
| seed rows | 92,781 |
| v3 pair-level projection parses | 92,781 (100.00%) |
| v3 parse failures | 0 |
| pairs, v1 effective | 1,669,314 |
| pairs, v2 effective | 2,029,886 |
| pairs, v3 | 2,517,623 |
| pairs added v3 vs v2 effective | +487,737 |
| rows with 0 v3 pairs | 4,245 |
| rows with 1–3 v3 pairs | 1,127 |
| rows with <4 v3 pairs | 5,372 (5.79%) |

On the 83,459 rows where both v2 and v3 parse, the pair sets differ on 6,757
rows: 27,367 v2-only pairs and 58,637 v3-only pairs. The expected “identical”
result is false: v2 can fabricate a pair after a gap, and can lose a genuine
pair when the projected stack is disturbed. v3 is the specification.

## Per-split counts

Counts below come from the server membership files
(`--split-dir splits/s1prime`). `n_seed_rows` is the number of records that map
to an Rfam seed label; the rest are bpRNA rows, whose labels are unchanged.

| split | records | seed rows | v1 empty truth (seed) | v1 <4 (seed) | v3 parse ok | v3 empty | v3 <4 | v3 excluded |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `train` | 17,748 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| `id_test` | 1,409 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| `test_withheld` | 17,371 | 11,527 | 2,392 | 2,601 | 11,527 | 112 | 341 | 1.96% |
| `test_novel` | 42,873 | 42,873 | 10,470 | 10,743 | 42,873 | 3,675 | 3,675 | 8.57% |

The withheld seed subset is the population the PI flagged: **20.75% of its v1
labels were empty truth**. Under v3, 112 rows are genuinely empty and 341
(2.96% of seed rows) fall below the Amendment 8 `<4` threshold. The full
`test_novel` split loses 8.57% under the same rule, below the 10% escalation
gate; the frozen 2,372-sequence eval subset loses 0/2,372.

## Reproduction

```bash
# local (no membership lists required)
python scripts/build_labels_v3.py

# server, with per-split counts
/home/junxuanliu/.conda/envs/rff/bin/python scripts/build_labels_v3.py \
    --split-dir splits/s1prime
```

`tests/test_labels_v3.py` checks the artifact hash, the total row count, the
frozen eval-subset counts, and the server-side per-split counters.
