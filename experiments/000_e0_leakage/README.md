# E0 — Leakage audit of standard RNA secondary-structure splits

## Question

How much sequence-level overlap exists between the standard bpRNA-1m training
split (TR0) and commonly used test sets (bpRNA VL0/TS0, ArchiveII, bpRNA-new)?

## Data

- Source: BPfold data release v0.2 (`BPfold_data.tar.gz`),
  `https://github.com/heqin-zhu/BPfold/releases/download/v0.2/BPfold_data.tar.gz`
- Size: 166,902,534 bytes
- SHA-256: `819acdea53c6b093df14cab35bf079c4678325d2d5ea4296b904068d35959403`
- Local path (server): `data/raw/bpfold/BPfold_data/`

This release is used for pipeline development and baseline reproduction (T4).
Claims will be re-run on raw-source splits (S1′-S3) as specified by the PI.

## Method

- Sequences normalized (upper-case, T -> U).
- Exact duplicates: SHA-256 of the full normalized sequence.
- Near-duplicates: shared 15-mers, candidates from an inverted index, verified
  with containment = shared / min(|kmers_query|, |kmers_reference|) >= 0.5.
- Family labels parsed from bpRNA identifiers (`bpRNA_RFAM_<id>`) or from the
  ArchiveII filename prefix.

## Commands

```bash
PYTHONPATH=src python scripts/audit_splits.py \
  configs/e0_bpfold_audit.json \
  --out experiments/000_e0_leakage/report
```

## Caveats

- Family labels are not comparable across bpRNA (`RFAM_*`) and ArchiveII
  (functional prefixes), so the family-overlap columns are not yet meaningful.
  A common Rfam 15.0 mapping via cmscan is the next step.
- The near-duplicate test is a conservative k-mer heuristic, not an alignment.
  It can miss remote homology and can flag conserved short motifs; sequence-level
  hits are inspected before any claim.
- Lengths above 500 nt are outside the main E1′ scope and are reported separately.

