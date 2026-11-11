# Data sources and mirrors

Every external dataset used by this project, with the exact locator and a
checksum, so that a local copy can be verified rather than trusted.

Large files stay out of git (`data/` is ignored); this file is the record.

## ArchiveII (evaluation benchmark)

| field | value |
|---|---|
| source | HuggingFace `rouskinlab/archiveII` (`data.json`) |
| mirror used locally | `https://hf-mirror.com/datasets/rouskinlab/archiveII/resolve/main/data.json` |
| size | 3,047,984 bytes |
| sha256 | `d43276a64beb2d992aa40b060df48ac6dca288c67a1c3ed32a16fc3fce96e28a` |
| records | 3,370 |

Format: JSON object keyed by record id (`<family>_<name>-<n>`); each value has
`sequence` and `structure` (a list of 0-based `[i, j]` pairs).

Family counts: 5s 1107, srp 697, tRNA 493, RNaseP 429, tmRNA 404, grp1 93,
16s 69, telomerase 37, 23s 30, grp2 11.

**Caveat, recorded before use.** This is a third-party mirror of ArchiveII, not
the Doty lab distribution. It contains 3,370 records while the cmscan input in
E0 held 3,966; the difference has not been explained. Do not quote counts from
this file as if they were the canonical distribution without checking against
the cluster copy.

## bpRNA-1m (training pool, 90% identity clustering)

| field | value |
|---|---|
| source | HuggingFace `rouskinlab/bpRNA-1m` |
| pinned revision | `8b973adeba0218113a696dd1c22ada682f26324a` (2024-03-17) |
| mirror used locally | `https://hf-mirror.com/datasets/rouskinlab/bpRNA-1m/resolve/<rev>/data.json` |
| size | 103,163,072 bytes |
| sha256 | `d36c62ee92eda8b6db90ac4cd6e9d9804aedbaa8a860dd1fcd9b116d5475dff4` |
| records | 66,715 |

Verified 2026-11-12: the local download is byte-identical to the pinned hash
recorded in `experiments/000_e0_leakage/bprna_mirror/README.md`. The mirror is a
strict superset of the BPfold bpseq labels (0 missing pairs across TR0/VL0/TS0);
its provenance limits are documented in that README.

## Rfam 15.0

| file | locator | size | sha256 |
|---|---|---|---|
| `Rfam.seed.gz` | `https://ftp.ebi.ac.uk/pub/databases/Rfam/15.0/Rfam.seed.gz` | 5,083,331 | `b06ac74d772609e863c5e1906eb491d699e9d325c5f22f4fc78126a808ad0bed` |
| `Rfam.clanin` | `https://ftp.ebi.ac.uk/pub/databases/Rfam/15.0/Rfam.clanin` | 5,418 | (not yet hashed) |

Sanity check: `Rfam.seed.gz` contains **4,178** family records, which matches the
family count recorded in `docs/research_state.md` (E0c).

Not yet mirrored: `Rfam.cm` (the covariance models, needed by `cmscan`). The
split rebuild in Amendment 4 §A4.5 requires it plus an Infernal binary.

## Network constraint observed 2026-11-12

`huggingface.co` is unreachable from this machine (TLS reset on :443), while
`hf-mirror.com`, `ftp.ebi.ac.uk`, `rna.urmc.rochester.edu` and `zenodo.org`
answer 200. All mirror URLs above are the ones that actually worked, not the
canonical ones. Any hash comparison is against the pinned value, never against
the mirror's own claim.
