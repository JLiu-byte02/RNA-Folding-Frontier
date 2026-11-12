# Round 4i reply — Claude Opus 5.5 (High), PI

Source: https://claude.ai/chat/cc2ad40a-406e-481e-8e7b-89a64b78f2f4
Captured: 2026-11-12, after the v3/D1/K1/D2 report.

## Round 4i: continue H2, but two cheap controls come before the P2-block run

v3 labels parsing at 100% with the embargo intact, the server under git with the
frozen hashes recorded, and K1 computed on regenerated probability dumps: that’s
a solid round. My rulings follow, but first a problem with the D1 result itself.

### The e effect on id_test can’t be what it appears to be until one control is run

Under condition e, for sequences of 80–299 nt, half or more of the sequence is
replaced with N, and the head sees nothing but FM embeddings. Yet id_test nested
F1 rises from 0.675 to 0.757, and that’s at the oracle threshold, so calibration
doesn’t explain it. There is a plausible way for this to happen: a template
prior.

- id_test consists of new sequences from training clans, which have stereotyped
  structures at stereotyped positions (tRNA, 5S, SRP).
- The head has a learned relative-position embedding.
- With rotary attention, the N-masked embeddings still encode distance to the
  kept window.

So under e the head can output the structure typical of a given length and
position without reading any sequence. That works in distribution and poorly on
unseen clans, which matches the pattern of +0.08 on id_test against +0.01 on
withheld. If that’s right, the id_test e effect is a shortcut of its own, based
on geometry rather than sequence, and it says nothing about distal context
misleading the model. Two controls on the existing checkpoints, under 1 GPU-h
together:

1. **Region-stratified D1, condition a vs e:** F1 separately for pairs with both
   ends in the kept window, one end, and neither. If e’s “neither” cell is well
   above zero on id_test, the head is predicting pairs it has no sequence for,
   which means a template.
2. **Sequence-free inputs:** an all-N input of the true length, plus a
   whole-sequence dinucleotide shuffle, both evaluated with the native head on
   id_test and withheld. Each is one full-length cache, about 0.2 h. The all-N F1
   measures how much structure length and position alone recover.

How to read them:

- **High all-N F1 on id_test** (say ≥ 0.3) and a nonzero “neither” cell:
  P2-central’s e condition measures template fallback. H2 is then judged only on
  kept–kept cells and in P2-block.
- **All-N near zero:** the e effect is real, and the spurious-partner reading
  gains support.

### H2’s status under the pre-registered rules

| Item | Status | Note |
|---|---|---|
| K1 | Not fired | Fires only if Δ < 0.01 and the CI includes 0. Withheld is 0.0119, so it narrowly misses. Record it as a near miss, not support. |
| H2.1 | Partial replication | Holds on withheld (3.89× per candidate cell), fails the 2× threshold on id_test (1.62×) once normalised. Normalising per candidate cell was the right call; note it was added after the amendment. |
| H2.2 | Supported | Long-range false positives sit in stacked helices about 4× more often than the null, on both splits. |

So H2 continues. H2.2 raises an issue for H2.4: if those false-positive helices
are energetically competitive alternative structures, a Turner-consistent
decoder will keep them, and H2.4 could fail for a principled reason. Add this
exploratory check on id_test and withheld (CPU only):

- Compute the stacking free energy of each false-positive helix with RNAeval on
  the helix alone, and compare its distribution with that of true helices.
- Count the fraction of false-positive helices that also appear in the RNAfold
  MFE structure.

If most of them appear in the MFE, expect A-energy-fixed to fail on long-range
false positives, and record that expectation before H2.4 runs.

### Your spending question: neither option as posed

Retraining heads on v3 labels isn’t needed. Training is bpRNA only, so v3 doesn’t
change any training label. Keeping crossing pairs as training positives is a
method choice, and evaluation already uses the nested subset. Retraining happens
only where already planned: the P3 train-and-test heads and A-energy-learned.
Order, using both resource types in parallel:

1. **GPU, about 1 h:** the two controls above. They decide how P2-block results
   get read.
2. **CPU, starting now:** dump the native head’s logits on the train split and
   fit `(s, b)` under Amendment 7. Freeze and hash them, then run H2.4 on id_test
   and withheld (exploratory there). This is the claim the paper depends on, and
   it doesn’t depend on P2-block.
3. **GPU, about 20 h, after step 1:** the P2-block cache for conditions a–e, then
   H2.3. If step 1 shows a strong template prior, add the all-N cache to P2-block
   as condition f, so the template baseline is measured inside the confirmatory
   geometry.
4. **Then:** train-and-test heads, the positive control, baseline retraining,
   A-energy-learned.

### One check on v3

v3 being correct by construction doesn’t settle the v2/v3 disagreement on 6,757
rows. There is a free cross-check: sequences that appear both in bpRNA-1m’s
Rfam-sourced entries and among the Rfam 15.0 seeds, as exact sequence matches.
Their bpRNA labels were projected by another group from an older Rfam release, so
agreement won’t be perfect, but v3 should agree with them clearly better than v2
does.

- Report pair-level F1 of v2 and v3 against the bpRNA label on that overlap.
- Inspect 10 of the rows where v2 and v3 differ by hand.

### Next report

1. The template-prior controls: region-stratified a/e and all-N/shuffle F1 on
   both splits.
2. The frozen `(s, b)` and the H2.4 verdict on id_test and withheld, with its
   basis (recall-matched, pair-budget or null).
3. The false-positive helix energetics.
4. The v3 cross-check.

Lead with item 1; it determines how much of Round 4 still stands.
