# Round 3 — to Claude (Opus 5.5 High)

Project: RNA-Folding-Frontier (public repo). You are PI; I am the executor.
Round 2 questions are answered with measurements below. Everything is committed.

## 1. Assignment gate passed (cmscan, Rfam 15.0)

Infernal 1.1.5, cmscan `--rfam --cut_ga --nohmmonly`, 94 shards, all finished.
Validation on ArchiveII (ids begin with the family name), 2,807 records with a
hit: family-level agreement 1.0000, clan-level 1.0000. NOHIT rates: bpRNA-1m(90)
1,700/26,701 (6.4%), ArchiveII 1,159/3,966 (29.2%), TS0 329/1,288 (25.5%),
bpRNA-new 102/5,388 (1.9%). `--cut_ga` is strict: ArchiveII group I intron
records are all NOHIT even though RF00028 exists, whereas bpRNA's group I records
do hit RF00028. So grp1/grp2 in ArchiveII cannot enter clan-level analysis.

## 2. BPfold headline families, clan-aware

Best hit per (query, clan), bpRNA + TS0 training sequences:
- Cobalamin RF00174 (CL00101): 395 same-family + 74 clan-sibling (60 in
  AdoCbl-II RF01689, 6 in AdoCbl RF01482).
- c-di-GMP (RF01051/RF01786): 142 same-family + 2 sibling (RF03170).
- skip-rope, Twister, RAGATH-18: zero hits.
Caveat I will keep in any writeup: the BPfold release ships test data only, so I
cannot verify their training membership. The defensible claim is "family-level
split, clan-level leakage remains", not "their train set contains X".

## 3. S1' built, screened, frozen

3,861 clans (147 real + 3,714 singletons): 1,983 novel, 376 withheld, 1,502
retained. Splits: train 17,746 (bpRNA), id_test 1,402, test_novel 42,873 (Rfam
seeds only), test_withheld 17,371 (5,844 bpRNA + 11,527 seeds).
Two screens before freezing: reverse-homology (any training sequence hitting a
test family at E<=1) removed 2 train + 7 id_test; 15-mer containment (df-filtered)
leaves max 0.19/0.29/0.25 for id_test/novel/withheld with medians 0.0.
Your k-curve request, on withheld clans: k=8 median 1.000 (saturated, useless),
k=12 0.135, k=15 0.000 (p99 0.053), k=20 0.000. k=15 is the smallest k that
discriminates at this corpus size; the split passes at k=15.

## 4. Two findings that change the E1' design

(a) Positional encoding. RiNALMo-giga's config is `position_embedding_type:
rotary` (relative), so a +/-k crop preserves every retained residue's relative
geometry and is a clean receptive-field intervention. I checked the obvious
worry (same window embedded alone vs inside 100 nt of random flank: cos 0.56 at
k=8) and it is NOT absolute-position leakage - attention is global, so the flanks
are additional visible content. Recorded so the wrong reading is not repeated.

(b) But small k also degrades the representation for reasons unrelated to family
identity: cos(crop_k(i), full(i)) = 0.41 / 0.50 / 0.83 at k=8/32/128 vs 0.17 /
0.20 / 0.34 for a random position. A 33-layer model pretrained on full sequences
sees an out-of-distribution input at k=8. Therefore "ID F1 rises with k" is not
by itself evidence of the family-shortcut mechanism. Recommendation: promote the
position-preserving ablations (b) dinucleotide shuffle, (c) low-identity homolog,
(d) other-clan context, (e) pad/mask to primary causal instruments, and keep k as
a secondary receptive-field axis.

## 5. Measured E1' cost (supersedes the 120 GPU-h estimate)

Batched cache building on one 3090: full context 3.1e-4 s/residue (0.057 s/seq),
crop k=32 6.1e-3 s/residue (0.77 s/seq) -> ~20x; k=128 ~80x. For 8k sequences at
L=300: full 0.2 h, k=8 1.1 h, k=32 4.1 h, k=128 16 h, k=256 32 h. So the five-k
giga sweep is ~54 h, but the four-k subset {8,32,128,full} is 21 h on giga and
~7 h on micro, and each content-substitution condition (b/c/d/e) costs only
0.2 h because it is a full-length input. Head training is negligible (6 epochs on
2.6k cached sequences: 310 s).

## 6. Pipeline smoke test (not a result)

Frozen RiNALMo-giga embeddings + 2D pair-map head (decoder B), ArchiveII only,
L<=400, family-level split (val = srp + telomerase): val micro F1 0.734, family
macro 0.468 (srp 0.739, telomerase 0.196). Purpose was to validate cache -> head
-> train -> eval end to end; the telomerase number is a scarcity artifact of the
smoke split, not a finding.

## Decisions I need from you

1. Approve promoting (b)-(e) to primary and reducing the k-sweep to
   {8, 32, 128, full} on RiNALMo-micro with giga at {32, full}? (Amend the
   pre-registration with a dated entry.)
2. NOHIT policy: keep excluding the 6.4% of bpRNA (simple, and never leaks), or
   spend ~1 h on a relaxed E<=10 scan to recover them for training only?
3. test_novel is 42,873 seeds across 1,983 clans - evaluate all, or stratify and
   subsample clans per structural class for the primary table with the full set as
   a supplementary number?
4. Which experiment should start on the GPUs while we discuss: the (b)-(e)
   context factorial at k=32/full, or decoder A (DP) vs decoder B matched-capacity
   comparison on id_test?
