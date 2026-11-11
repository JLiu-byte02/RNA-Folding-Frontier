# bpRNA-1m provenance (HuggingFace mirror)

The authors' `bpRNA_1m_90.zip` is 4.02 GB and served at ~8 KB/s from
`bprna.cgrb.oregonstate.edu`, which is not viable. We therefore use the
`rouskinlab/bpRNA-1m` mirror and document exactly what it is.

## What the mirror is

- Source: <https://huggingface.co/datasets/rouskinlab/bpRNA-1m>
- Pinned revision: `8b973adeba0218113a696dd1c22ada682f26324a` (2024-03-17)
- File: `data.json`, 103,163,072 bytes
- SHA-256: `d36c62ee92eda8b6db90ac4cd6e9d9804aedbaa8a860dd1fcd9b116d5475dff4`
- Records: 66,715 (the mirror's own conversion report says 102,318 input
  datapoints, 5,064 invalid, 1,482 duplicate sequences with differing
  structure/DMS/SHAPE dropped)

## Is it faithful? Measured, not assumed

Containment of the BPfold bpRNA release inside the mirror
(`scripts/check_bprna_mirror.py`):

| split | n | sequence found | bpseq pairs ⊆ mirror pairs | identical after canonical filter | missing bpseq pairs |
|---|---:|---:|---:|---:|---:|
| TR0 | 10,814 | 98.78% | 100% | 91.2% | 0 |
| VL0 | 1,300 | 98.92% | 100% | 92.4% | 0 |
| TS0 | 1,305 | 98.70% | 100% | 90.8% | 0 |

Interpretation: the mirror is a strict **superset** of the bpseq labels. It never
drops a pair that bpseq has (0 missing across 13,419 sequences). Where it differs
it *adds* pairs: 42,106 extra pairs across TR0, of which 80.6% are non-canonical
(G-A, C-U, G-G, ...). That is consistent with the mirror keeping bpRNA's full
3D-derived annotation while `.bpseq` can only encode canonical nested pairs.
Filtering the mirror to canonical pairs plus GU reproduces the bpseq labels
exactly for 91.2% of TR0.

## Pair policy adopted

Canonical pairs (AU, UA, GC, CG) plus GU/UG; non-canonical pairs and pseudoknots
are recorded but excluded from primary F1. With that policy the mirror and the
bpseq release agree, and the mirror is usable as the training pool.

## What we could not determine

- The exact filter the mirror's converter applied beyond the counts above.
- Whether any source-database field was dropped (the mirror keeps only id,
  sequence and pairs; bpRNA's source/species/family annotations are not present,
  which is why family assignment must come from cmscan).

## Re-derivation of a 90% subset

`cd-hit-est -c 0.9 -n 5` (CD-HIT 4.8.1) on the mirror's sequences, because
MMseqs2 binaries segfault against this host's glibc/libstdc++. This is a
documented re-derivation, not the authors' `bpRNA_1m_90` file; the leakage
control that matters is clan holdout plus containment filtering, not the
clustering step.
