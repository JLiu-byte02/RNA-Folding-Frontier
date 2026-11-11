# P2 — test-time distal-context intervention (RiNALMo-giga, decoder B)

Status: **first complete run.** Three seeds (17/23/41), decoder-B pair head trained
on native inputs only (condition a, 8 epochs, family-disjoint internal validation
slice), then evaluated with the distal region rewritten at test time. Nothing
about the weights changes between rows.

Artifacts: `runs/p2_seed{17,23,41}/results.json` (server),
`p2_analysis.json`, `diagnostics_{id_test,test_withheld}.json`,
`diag_*_perseq.jsonl`.

## The table

Family-macro F1 is the pre-registered primary metric; Δ is the paired family
bootstrap (10,000 resamples) against **the same split's condition-a baseline**.
`c` is evaluated on the subset that has a verified <60%-identity homolog, so its
baseline is `a_csub` on that subset, not `a` on the full split.

### id_test (1,392 sequences, 168 families)

| condition | n | micro F1 | family-macro F1 | Δ macro vs native [95% CI] |
|---|---|---|---|---|
| a native | 1392 | 0.6562 | 0.6662 | — |
| b dinucleotide-shuffled | 1392 | 0.6251 | 0.6936 | +0.0274 [−0.0144, +0.0706] |
| c homolog <60% identity | 532 | 0.4685 | 0.5460 | **−0.0915 [−0.1412, −0.0402]** |
| d other clan | 1392 | 0.6553 | 0.7190 | **+0.0528 [+0.0086, +0.0977]** |
| e N-masked | 1392 | 0.7103 | 0.7549 | **+0.0888 [+0.0448, +0.1342]** |

`a_csub` (the 532-sequence subset that has a homolog) is 0.6518 micro / 0.6352
macro.

### test_withheld (16,793 sequences, 442 families)

| condition | n | micro F1 | family-macro F1 | Δ macro vs native [95% CI] |
|---|---|---|---|---|
| a native | 16793 | 0.5920 | 0.6028 | — |
| b dinucleotide-shuffled | 16793 | 0.5305 | 0.6108 | +0.0080 [−0.0119, +0.0290] |
| c homolog <60% identity | 8558 | 0.4730 | 0.5335 | **−0.0757 [−0.1014, −0.0502]** |
| d other clan | 16793 | 0.5338 | 0.6097 | +0.0069 [−0.0135, +0.0279] |
| e N-masked | 16793 | 0.5695 | 0.6433 | **+0.0405 [+0.0190, +0.0633]** |

## What it says

1. **Replacing the distal region with plain "unknown" (e) improves family-macro
   F1 on both splits**, +0.089 (ID) and +0.041 (OOD), with intervals that exclude
   zero. The native distal context is therefore a net *negative* for most
   families, not a helpful prior.
2. **Injecting plausible-but-wrong distal content is worse than injecting
   nothing.** Other-clan (d) and shuffled (b) sequences land below (e) almost
   everywhere in micro terms; the model does not ignore the swapped sequence.
3. **(c) is the single worst condition, significantly so on both splits.** A
   low-identity same-family homolog, spliced into the distal region, costs 0.09
   (ID) and 0.08 (OOD) family-macro F1 relative to its own baseline. This is the
   opposite of the pre-registered P2 prediction that homolog context should
   *retain* the native advantage.
4. **Length decides the sign.** By band, e − a is +0.08 at 80–149 nt and +0.11 at
   150–299 nt, but −0.03 at ≥300 nt (micro, id_test means over seeds):

| band | a | b | d | e |
|---|---|---|---|---|
| L<80 | 0.748 | 0.731 | 0.738 | 0.750 |
| 80–149 | 0.671 | 0.633 | 0.667 | 0.744 |
| 150–299 | 0.543 | 0.584 | 0.618 | 0.657 |
| L≥300 | 0.550 | 0.430 | 0.479 | 0.518 |

   Below 300 nt the distal context mostly contributes a misleading prior; above
   300 nt it also carries genuine long-range structure, and removing it hurts.

## Against the pre-registration

- **P2 as written is not supported and cannot be evaluated as phrased.** Its
  criterion is "(c) retains ≥70% of the (a)−(e) ID gain; (d) retains ≤30%".
  The (a)−(e) gain is *negative* on both splits (−0.089, −0.041), so a retention
  ratio is undefined, and (c) goes the wrong way by a wide margin.
- **P3 (OOD: b ≥ a) is not supported**: b is slightly *below* a on withheld
  (+0.008 macro, −0.062 micro).
- **The Amendment 1.3 decodability control has almost no power here.** The clan
  probe reaches 0.0526 balanced accuracy on native id_test embeddings against a
  0.042 chance level over the 24 clans shared with training, and it is
  *undefined* on withheld because no clan is shared by construction. Per
  sequence, distortion explains 5–11% of F1 variance within a condition, and
  adding the probe's p(true clan) term adds only 0.2–0.3 points of R². So the
  clause "clan decodability predicts F1 beyond distortion" fails, but mostly
  because the quantity carries no signal on these splits, not because distortion
  explains everything.

## Why the probe is at chance, and why that matters

`probe_sanity_clan_train.json`: the same probe, trained and tested on a random
80/20 split of *training* sequences, reaches **0.9216** accuracy over 590 clans
(chance 0.0017; balanced 0.8571; nearest-centroid 0.7315). The representation
therefore carries identity information, but it is **family-specific**: the split
holds out whole families inside retained clans, and on those held-out families
the probe collapses to chance.

That is a statement about the split, not about the model: with family-disjoint
training and evaluation, a family-identity shortcut has nothing to exploit. If
we want to test "distal context is a family shortcut", it has to be tested on a
split where near-duplicates survive (a sequence-level split), which is exactly
what E0 showed the field's benchmarks do.

## Caveats

- The head is small (8 dilated conv blocks, 32 channels) and is *not* the
  final decoder. Absolute F1 is not a benchmark number.
- `test_withheld` native F1 (0.592 micro) is lower than `id_test` (0.656) as
  expected for a clan-held-out split; N differs slightly (16,793 of 17,371)
  because 578 sequences exceed the 1,024-nt cache limit.
- 10.3% of training and 8.5% of withheld sequences are shorter than the ±32
  window, so their distal region is empty and conditions a/b/d/e are identical
  for them. They dilute every effect toward zero; the by-band table is the
  honest read.
- Condition (c) covers 39% of id_test and 49% of withheld — only records with a
  verified 0.40–0.60 identity homolog. Those subsets are not representative.
- Condition (e) uses the `N` token and is a distortion anchor, not a natural
  input; its F1 should not be quoted as a method result.
