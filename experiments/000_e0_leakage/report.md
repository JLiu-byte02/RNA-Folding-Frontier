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

## E0d — The other trainer: ArchiveII vs RNAStrAlign

Verified from the BPfold paper text: trainable models there are trained on
**RNAStrAlign (19,313 sequences, <= 600 nt) *and* bpRNA-1m** (TR0, 12,114;
released as TR0 10,814 + VL0 1,300). Auditing only bpRNA therefore audits half
the training pool.

| test set | n | exact duplicates vs RNAStrAlign | near-duplicates (k=15, containment >= 0.5) |
|---|---:|---:|---:|
| ArchiveII | 3,966 | **2,005 (50.55%)** | **3,246 (81.85%)** |
| bpRNA TS0 | 1,305 | 66 (5.06%) | 149 (11.42%) |

Read together with E0's bpRNA numbers, ArchiveII is more than 50% exact-overlap
in-distribution under BPfold's own protocol, and 81.85% near-duplicate to the
RNAStrAlign trainer alone.

That does **not** contradict BPfold: the paper explicitly calls ArchiveII and
TS0 "sequence-wise" datasets, reserves "cross-family" for Rfam 12.3-14.10
(10,791 RNAs), and never calls ArchiveII out-of-distribution. The correct claim
is about the *field*: ArchiveII is routinely reported as a generalisation
benchmark while being in-distribution for any model trained on the standard
RNAStrAlign + bpRNA pool.

Caveat: this is measured against the released RNAStrAlign directory (37,052
files), not against BPfold's filtered 19,313-sequence training subset, so the
training-subset overlap is <= these numbers. Recomputation against the exact
filtered subset belongs to T4.

## E0e — Headline: ArchiveII against the union of both trainers

| quantity | value |
|---|---:|
| ArchiveII sequences | 3,966 |
| exact duplicates of bpRNA TR0 | 812 (20.47%) |
| exact duplicates of RNAStrAlign | 2,005 (50.55%) |
| **exact duplicates of the union** | **2,427 (61.20%)** |
| **near-duplicates of the union** (k=15, containment >= 0.5) | **3,506 (88.40%)** |
| duplicate sequences *inside* ArchiveII | 313 sequences, 515 extra copies |

Per-family, the contamination is concentrated in the short families and absent
in the long ones, because the standard trainers cap length at 600 nt:

| family | n | exact TR0 | exact RNAStrAlign | exact union | near union |
|---|---:|---:|---:|---:|---:|
| tmRNA | 462 | 173 | 432 | 444 (96.1%) | 455 |
| 5s | 1,283 | 85 | 852 | 867 (67.6%) | 1,197 |
| srp | 928 | 201 | 601 | 673 (72.5%) | 804 |
| tRNA | 557 | 175 | 97 | 250 (44.9%) | 461 |
| RNaseP | 454 | 134 | 22 | 148 (32.6%) | 449 |
| grp1 | 98 | 34 | 1 | 35 (35.7%) | 91 |
| telomerase | 37 | 10 | 0 | 10 (27.0%) | 37 |
| 16s | 106 | 0 | 0 | 0 | 12 |
| 23s | 30 | 0 | 0 | 0 | 0 |
| grp2 | 11 | 0 | 0 | 0 | 0 |

The 313 within-ArchiveII duplicates matter for reporting: per-sequence means
count those sequences up to eight times, so sequence-micro F1 over-weights a few
structures while family-macro F1 does not.

## E0f — Why overlap translates into F1: the labels are nearly identical

Same procedure as E0b, but comparing ArchiveII against **RNAStrAlign**:

| comparison | duplicates | identical annotations | mean pair F1 | median | p10 |
|---|---:|---:|---:|---:|---:|
| ArchiveII vs bpRNA TR0 | 812 | 378 (46.6%) | 0.9191 | 0.982 | 0.710 |
| ArchiveII vs RNAStrAlign | 1,828 | 1,667 (91.2%) | 0.9946 | 1.000 | 1.000 |

Per family (ArchiveII vs RNAStrAlign): srp 0.999, tmRNA 0.9999, 5s 0.991, tRNA
0.985, RNaseP 0.959, grp1 0.875; only 5 of 1,828 sequences fall below 0.8.

This is the mechanism that turns overlap into accuracy. RNAStrAlign and
ArchiveII share a curation source, so their labels agree at 0.995, whereas bpRNA
— a separately derived consensus annotation — agrees with ArchiveII at only
0.919. A model trained on RNAStrAlign can therefore reproduce ArchiveII's own
labels almost verbatim on the 50.55% of the test set it has memorised, and no
amount of held-out *sequence* identity prevents that. Held-out *clan* is the only
fix, which is now the S1′ design.

## Reproduce

```bash
python scripts/label_discordance.py \
  --archiveii data/raw/bpfold/BPfold_data/archiveII/archiveII \
  --out experiments/000_e0_leakage/label_discordance_bpfold_archiveII
python scripts/audit_splits.py --config configs/e0_bpfold_audit.json
python scripts/build_rfam_index.py
```

