# Round 4k reply — Claude Opus 5.5 (High), PI

Source: https://claude.ai/chat/cc2ad40a-406e-481e-8e7b-89a64b78f2f4
Captured: 2026-11-12, after the corrected P2-central table.

## Round 4k: authorized, with P2-block’s purpose reframed and the positive control deprioritised

The corrected table is clean, the invariants pass, and it’s good to know that the
E0 audit numbers and the independent scorer held up. Here is what the table says
before deciding what to run.

### Reading the corrected P2-central table (exploratory: central-window geometry)

Retention on kept–kept pairs at the fixed threshold, as defined in P2:
`(X − e)/(a − e)`.

| Condition | id_test | withheld | P2 criterion |
|---|---:|---:|---|
| c (homolog distal context) | (.5616−.5142)/(.5759−.5142) ≈ **77%** | (.5927−.5171)/(.5910−.5171) ≈ **102%** | ≥70% |
| d (other-clan distal context) | ≈ **−17%** | ≈ **8%** | ≤30% |

So the exploratory version meets the P2 numbers on both splits. Add bootstrap
CIs for both ratios to the P2-central report; ratios of differences can be
unstable.

This doesn’t support the shortcut hypothesis. The pattern points to something
different. D1’s core claim was that distal context encodes family identity,
which helps in distribution and hurts out of distribution. The table shows:

- **Distal context helps OOD too.** On withheld kept–kept pairs, a (.591) beats
  e (.517) and b (.528). P3’s prediction (b ≥ a on OOD) fails.
- **Homologous context substitutes fully for native context on clans the head
  never saw** (102% retention on withheld). Unrelated context does not.

The simplest reading: the FM uses distal context that is evolutionarily
consistent with the sequence (partner-like, structurally compatible), and that
transfers to clans the head wasn’t trained on. That’s closer to a legitimate
signal than to a shortcut.

One caveat keeps the shortcut alive in a different place. Withheld clans are
unseen by the head, but almost certainly not by the FM, which was pretrained on
RNAcentral. Any family-identity effect could live in the pretraining, which no
split of ours can remove. Note it now. It’s also a plausible account of
CHANRG’s failure (OOD on families the FM saw less of in pretraining), and could
become an exploratory analysis later.

**What would settle D1 is train-and-test, not more test-time swaps.** If heads
trained without distal context, or with scrambled distal context, don’t beat the
native-trained head on withheld, D1’s core claim is falsified as originally
stated.

### Authorizations

1. **P2-block cache: authorized**, extended to the train split for the decisive
   train-and-test arms.
   - Test-time arm: conditions a–e on id_test and withheld, scored on all pairs
     under the pre-registered P2 criteria.
   - Train-and-test arm (P3): heads trained and tested on condition a, on
     e-block (≈ a k=32 crop) and on b-block, 3 seeds each.
   - Budget: about 40 GPU-h across 6 GPUs, roughly a day.
2. **Amendment 10, pinned before the P2-block cache is read:**
   - **D1-core falsification rule:** D1’s core claim is falsified if neither
     e-block-trained nor b-block-trained heads beat the a-trained head on
     withheld all-pairs family-macro F1 by ≥0.02 with a CI excluding 0. This
     makes P3 decisive.
   - **Alternative hypothesis H3 (“evolutionarily consistent context
     transfers”)**, pre-registered now that we’ve seen the central-window
     pattern but before P2-block exists: in P2-block on withheld, c retains ≥90%
     of (a−e) and d retains ≤30%, both with CIs.
   - **Exploratory dose-response:** split condition c into homolog-identity bins
     (0.40–0.50 vs 0.50–0.60) from the existing c data. H3 predicts retention
     rises with identity.
3. **Positive control: deferred, not cancelled.** It was meant to show the
   instrument can detect dependence on family context, and c ≫ d on S1′ already
   shows that. Its only remaining value is to test whether the head uses clan
   identity beyond structurally compatible context in a near-duplicate setting.
   Revisit after P3.
4. **(s,b) refit with the retrained heads, then H2.4 on id_test and withheld:
   authorized now on CPU**, in parallel with the GPU work. Include the
   false-positive helix energetics check from Round 4i, and record its
   expectation for H2.4 before H2.4 runs. Order: `(s,b)` freeze and hash,
   energetics expectation recorded, H2.4.

### Where the paper stands

- D1-core falsified and H3 holds: the mechanism half becomes “FM failure OOD is
  not caused by head-level family shortcuts; distal context is evolutionarily
  informative and transfers.” A useful negative result plus a positive one. The
  method half depends on H2.4 and P5, with D2 as fallback.
- D1-core supported (cropped or scrambled training beats native OOD): the
  original story is back, with a corrected instrument.
- Both fail: the paper becomes the audit plus H2.4, and the D2 decision gets
  made then.

### Next report

1. Amendment 10 hash.
2. The P2-block test-time table (a–e, all pairs, ID and withheld) with P2 and H3
   ratios and their CIs.
3. The P3 train-and-test table, with the D1-core verdict.
4. The frozen `(s,b)`, the energetics expectation, and the H2.4 verdict on
   id_test and withheld, with its basis (recall-matched, pair-budget or null).
5. CIs for the P2-central ratios, plus the dose-response split.

Lead with the D1-core verdict. `test_novel` stays embargoed until Amendment 10’s
outcomes are recorded.
