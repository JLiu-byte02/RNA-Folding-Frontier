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
- **The family/clan columns in `audit.md` are currently vacuous.** bpRNA records
  carry no family metadata, so "shared families: 0" there means "unmapped", not
  "no overlap". Clan numbers must be treated as missing until cmscan is run.

## Artifacts

- `audit.md`, `summary.json` — machine-readable summary (full hit lists remain on
  the server under `experiments/000_e0_leakage/report/`).

## E0b — label agreement on the exact duplicates

For the 812 ArchiveII sequences that appear verbatim in bpRNA TR0 we can compare
the two annotations directly. This bounds how well *any* predictor can score
against either label set on that subset, and it is free of model error.

| family | n | mean pair F1 | median pair F1 | identical annotations |
|---|---:|---:|---:|---:|
| tRNA | 175 | 0.9971 | 1.000 | 174 |
| srp | 201 | 0.9703 | 1.000 | 174 |
| grp1 | 34 | 0.9649 | 0.967 | 0 |
| RNaseP | 134 | 0.9340 | 0.960 | 9 |
| 5s | 85 | 0.8768 | 0.957 | 18 |
| tmRNA | 173 | 0.7960 | 0.863 | 3 |
| telomerase | 10 | 0.6606 | 0.689 | 0 |
| **all** | **812** | **0.9191** | **0.982** | **378 (46.6%)** |

Direction of the disagreement: bpRNA has *fewer* base pairs than ArchiveII in 394
of 812 cases, more in 30, and the same count in 388. Mean precision (bpRNA scored
against ArchiveII) is 0.966, mean recall 0.886. Median shifted-register pairs is
0, so the disagreement is not a small alignment offset; bpRNA is systematically
missing pairs, most likely pseudoknots and family-consensus projections that the
comparative model did not resolve.

Consequence: on the overlapping subset, F1 is capped near 0.92 on average and
near 0.66 for telomerase. A model reporting F1 above these numbers on those
sequences is partly memorising one annotation convention. This is the seed for
the D3 evidence-stratified evaluation.

## E0c — ArchiveII's nine "families" are eight clans

Rfam 15.0 ships 4,178 families, 147 clans and 463 clan memberships; 3,715
families are singletons. Mapping ArchiveII's families through the release's
`clan_membership.txt`:

| ArchiveII family | Rfam family | clan |
|---|---|---|
| 5s | RF00001 | CL00113 |
| 16s | RF00177 | CL00111 |
| 23s | RF02541 | CL00112 |
| tRNA | RF00005 | **CL00001** |
| tmRNA | RF00023 | **CL00001** |
| RNaseP | RF00010 | CL00002 |
| srp | RF00017 / RF00169 | CL00003 |
| telomerase | RF00025 | CL00004 |
| grp1 | RF00028 | singleton |
| grp2 | RF00029 | singleton |

tRNA and tmRNA share a clan, so ArchiveII leave-one-family-out has a
same-clan sibling in the training pool for those folds. Every other fold is
cross-clan. Any S2 fold must therefore drop the whole clan from training, not
just the named family, and the tRNA/tmRNA pair is the specific case where a naive
implementation leaks.

## Reproduce

```bash
python scripts/label_discordance.py \
  --archiveii data/raw/bpfold/BPfold_data/archiveII/archiveII \
  --out experiments/000_e0_leakage/label_discordance_bpfold_archiveII
python scripts/audit_splits.py --config configs/e0_bpfold_audit.json
python scripts/build_rfam_index.py
```

