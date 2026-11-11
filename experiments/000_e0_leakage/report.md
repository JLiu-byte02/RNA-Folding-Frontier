# E0 report — Standard splits are not clean

Date: 2026-11-11

## Result

Reference training split: bpRNA-1m TR0 (10,814 sequences).

| test set | n | exact duplicates vs TR0 | near-duplicate sequences (k=15, containment >= 0.5) |
|---|---:|---:|---:|
| bpRNA VL0 | 1,300 | 0 (0.00%) | 4 (0.31%) |
| bpRNA TS0 | 1,305 | 0 (0.00%) | 6 (0.46%) |
| ArchiveII | 3,966 | **812 (20.47%)** | **1,517 (38.25%)** |
| bpRNA-new | 5,401 | 0 (0.00%) | 3 (0.06%) |

## Interpretation

1. **ArchiveII is not an out-of-distribution benchmark when trained on bpRNA-1m.**
   One in five ArchiveII sequences is an exact sequence duplicate of a TR0
   training sequence, and more than one in three is a near-duplicate by a
   conservative 15-mer containment criterion. Any method trained on bpRNA-1m and
   evaluated on full ArchiveII is reporting a partially in-distribution number.
2. **bpRNA TS0/VL0 are sequence-clean relative to TR0**, with <0.5% near-duplicate
   hits. They remain useful for in-distribution and length-controlled analyses,
   but not for cross-family generalization claims.
3. **bpRNA-new is sequence-clean relative to TR0** at this threshold (0.06%),
   consistent with its design as a non-redundant test set. Family/clan overlap
   still needs an Rfam mapping before it can be called OOD.
4. The same audit applied to BPfold's protocol implies that its reported
   ArchiveII F1 (0.820) is measured on a test set with ~20% exact train overlap.
   This does not invalidate BPfold's method, but it weakens the OOD claim and
   motivates the rebuilt S1′/S2 splits.

## Limitations

- Cross-dataset family/clan overlap is unresolved because the naming schemes are
  disjoint. Next: cmscan every split against Rfam 15.0 CMs and map families to
  clans.
- The 15-mer criterion is not a substitute for structural homology detection.
- Only bpRNA TR0 is used as reference here; S1′ will use a clan-aware training
  pool, not TR0.

## Artifacts

- `audit.md`, `summary.json` — machine-readable summary (full hit lists remain on
  the server under `experiments/000_e0_leakage/report/`).

