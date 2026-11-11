# Amendment 2 — protocol deviation, P2-block, and H2 pre-registration

Date: 2026-11-12. Source: PI Round 4 and Round 4b
(`docs/claude/transcripts/2026-11-12_round04_decisions.md`).
The frozen text above stays as written; this section overrides it where they
differ. Committed and hashed **before** any diagnostic dump was produced.

## A2.1 Recorded protocol deviation: P2-central is not the pre-registered P2

Round 1 specified a **per-position** local window: every position ``i`` is
embedded with its own native ``±32`` nt, so the intervention removes only what is
distant from that position. The implemented swap conditions keep **one central**
window ``[L//2-32, L//2+32]`` and rewrite everything outside it. The PI approved
the "0.2 h per condition because it is a full-length input" cost estimate without
catching that the estimate is only correct for a single window.

Consequences, which are binding on how the result may be read:

- The Round-4 table is renamed **P2-central (exploratory)**.
- With a head that sees only frozen FM embeddings, any pair touching a rewritten
  position is predicted from the wrong sequence or from none. Therefore the
  **kept–kept cell of the region stratification is the only valid P2-central
  readout**; the other cells are diagnostics about the rewritten region only.
- One central window also changes the length profile of the intervention: the
  rewritten fraction is ``1 - (2*32+1)/L``, i.e. 19% at L=80 but 78% at L=300.
  Comparisons across bands are therefore not comparisons of the same
  manipulation.

No number produced under P2-central may be reported as a test of P2.

## A2.2 P2-block — the pre-registered intervention, as specified in Round 1

- Partition each sequence into blocks of 16 positions.
- For each block, keep the block plus ``±32`` nt on either side **native**, and
  rewrite everything else with the condition's donor content (b/c/d/e).
- Run one forward pass per block and take embeddings **only for that block's own
  positions**.
- Every position is therefore embedded with at least ±32 nt of its own native
  flank, which is the Round-1 design. Condition e under this geometry is
  approximately a k=32 crop.
- Scoring is unchanged: F1 over all pairs, not a subset.

Cost: ``ceil(L/16)`` passes per sequence — roughly 1 s per sequence at L=300,
less at the id_test median of 112 nt; about 4–6 GPU-h per condition over id_test
and withheld, ~20 GPU-h in total.

## A2.3 How pairs are extracted from the probability map (as implemented)

This was asked before the diagnostics ran, so it is fixed here:

1. Logits are symmetrised, ``0.5 * (L + L^T)``, then passed through a sigmoid.
2. Candidates are cells with ``prob >= 0.5`` and ``j - i > 3``.
3. Candidates are accepted in decreasing probability order while (a) neither
   endpoint already has a partner and (b) the pair does not cross any
   already-accepted pair.
4. Output is sorted by first index.

So the primary decoder is **threshold + one-partner + non-crossing greedy**, not
a bare threshold. D1 reports both anyway: the constrained decoding above (the
numbers in the Round-4 table) and a bare-threshold variant, so the contribution of
the global constraint is measured rather than assumed.

## A2.4 H2 pre-registration — "spurious partners"

H2 was inspired by the P2-central table, so this table is its out-of-sample
prediction set. All tests use the **native head (condition a)** unless stated.
Thresholds below are fixed now and may not be tuned after looking at the dumps.

| id | prediction | test | passes if |
|---|---|---|---|
| H2.1 | Long-range false positives rise with length | raw-threshold FP rate at ``|i-j| > 128`` per sequence | rate for ``L >= 150`` is >= 2x that for ``80 <= L < 150`` |
| H2.2 | False positives are helix-shaped | among long-range FPs, fraction lying in a stretch of >= 4 contiguous stacked canonical or GU pairs | >= 2x a null of random canonical non-true pairs matched on length and distance |
| H2.3 | (c) adds false positives at the window boundary | in P2-block, extra FPs under (c) vs (a) that link a block's native window to rewritten sequence | >= 60% of the extra FPs |
| H2.4 | Global consistency removes them | constrained vs bare-threshold decoding (later decoder A) | long-range FPs fall by more than true positives at fixed 0.5 |
| **K1 (kill)** | Calibration explains "e helps" | oracle threshold + one-partner decoding, ``e - a`` in the 80–299 band | if ``< 0.01`` with a CI that includes 0, the "distal context is a net-negative prior" reading is dropped |

If K1 fires, it is the first line of the next report and H2 is not pursued.

## A2.5 test_novel stays embargoed

No model is scored on test_novel until H2's status is settled. If H2 survives,
its confirmation runs on test_novel under P2-block plus decoder A.

## A2.6 Primary unit for test_novel (clarification, not a change)

Round 1b pre-registered singleton families as their own clans, so the
published-structure subset has ~310 clan units (6 multi-family clans plus 304
singletons), above the ~150 gate. Primary metric: family-macro F1 bootstrapped
over clan blocks. Robustness: cluster the test families by test-vs-test cmscan
hits (E <= 1) into homology blocks and redo the bootstrap over blocks.
