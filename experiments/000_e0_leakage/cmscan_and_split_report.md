# E0g — uniform family assignment and the frozen S1′ split

Date: 2026-11-11. Infernal 1.1.5 (built from source), Rfam 15.0 CMs, cmscan
`--rfam --cut_ga --nohmmonly --noali`. 94 shards x 400 sequences, all completed.

## 1. Assigner validation

ArchiveII record ids begin with their family, so the assigner can be checked
without extra annotation. On the 2,807 ArchiveII records that receive a hit:

- family-level agreement with the id prefix: **1.0000**
- clan-level agreement: **1.0000**
- 1,159 records (29.2%) have no `--cut_ga` hit

`--cut_ga` is strict. `grp1` (group I intron) records are entirely NOHIT even
though the family exists (RF00028, singleton clan); bpRNA's own group I records
do get RF00028 hits, so this is a threshold/architecture effect, not a parsing
failure. Consequence: ArchiveII grp1/grp2 cannot be placed in a clan by this
assigner and are excluded from clan-level analysis; their family labels are
still usable for T4 reproduction.

An end-to-end smoke test (`scripts/cmscan_selftest.sh`) caught two real defects
before the full run: a truncated sharder and an incomplete expected-family list
(telomerase maps to RF00024, not RF00025).

## 2. Coverage

| source | sequences | with hit | NOHIT |
|---|---|---|---|
| bpRNA-1m(90) | 26,701 | 25,001 | 1,700 (6.4%) |
| ArchiveII | 3,966 | 2,807 | 1,159 (29.2%) |
| TS0 | 1,288 | 959 | 329 (25.5%) |
| bpRNA-new | 5,388 | 5,286 | 102 (1.9%) |

NOHIT is never treated as a family or clan. It is excluded from S1′ entirely,
which is the conservative reading of the PI's rule; recovering those sequences
for training would need a second, relaxed scan (E <= 10) and is deferred.

## 3. BPfold's headline "unseen families", checked at clan level

Reading best-hit-per-clan (not the single-label table, which hides clan
siblings) over bpRNA + TS0 training sequences:

| headline family | Rfam | clan | same-family hits | clan-sibling hits |
|---|---|---|---|---|
| skip-rope | RF02924 | singleton | 0 | 0 |
| Twister | RF02684/81 | singleton | 0 | 0 |
| c-di-GMP | RF01051, RF01786 | CL | **142** | 2 |
| RAGATH-18 | RF03064 | singleton | 0 | 0 |
| Cobalamin | RF00174 | CL00101 | **395** | **74** |

Cobalamin's clan contains RF01482 (AdoCbl) and RF01689 (AdoCbl-II); 74 training
sequences fall in those sibling families (60 in AdoCbl-II alone). c-di-GMP's
sibling hit is RF03170 (c-di-GMP-II-GCG).

Limit of the claim: the BPfold release ships test data only, so their exact
training membership cannot be checked. What is established is that bpRNA-1m
contains these sequences and that a *family-wise* split, which is BPfold's
documented protocol, does not remove clan siblings. The statement to make is
"family-level split, clan-level leakage remains", not "their training set
provably contains X".

## 4. S1′

Built by `scripts/build_s1prime.py` (pool = bpRNA-1m(90) union Rfam 15.0 seeds;
clans from the Rfam release; singletons are their own cluster).

| split | records | composition |
|---|---|---|
| train | 17,746 | bpRNA only |
| id_test | 1,402 | 10% of retained bpRNA families |
| test_novel | 42,873 | Rfam seeds from 1,983 clans with no bpRNA support |
| test_withheld | 17,371 | 376 clans: 5,844 bpRNA + 11,527 seeds |

Clans: 3,861 in total, 1,983 novel, 376 withheld, 1,502 retained.

`scripts/leakage_filter.py` then applies two screens:

1. **Reverse homology.** A training sequence whose best hit is elsewhere but
   which hits a test family at E <= 1 is removed. 2 train and 7 id_test removed.
2. **15-mer containment** against the whole training set, low-complexity k-mers
   (document frequency > 5%) excluded. After filtering: max containment 0.19
   (id_test), 0.29 (novel), 0.25 (withheld); medians 0.0.

The containment curve on withheld clans (the PI's request to fit the threshold)
shows why k=15 and not something smaller:

| k | median | p99 | max | frac > 0.9 |
|---|---|---|---|---|
| 8 | 1.000 | 1.000 | 1.000 | 1.00 |
| 12 | 0.135 | 0.319 | 0.630 | 0.00 |
| 15 | 0.000 | 0.053 | 0.253 | 0.00 |
| 20 | 0.000 | 0.005 | 0.117 | 0.00 |

At k=8 every sequence is "contained" (8-mers saturate at this corpus size), so a
k=8 screen would be vacuous. k=15 discriminates and passes: the held-out clans
have no near-duplicates in training.

Hashes of the frozen splits are in `splits/s1prime/manifest.json` and
`splits/s1prime_filtered/filter_report.json`.
