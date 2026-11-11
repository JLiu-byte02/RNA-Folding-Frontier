# Next Actions

Updated: 2026-11-11 (after PI Round 3)

The pre-registration is now `PREREGISTRATION.md` plus Amendment 1. The amendment,
not the frozen text, governs: swaps are primary evidence, the k-sweep is
descriptive, P2 is a test-time intervention while P3/P5 are train-and-test per
condition.

## Running now

- [ ] `cache_emb.sbatch` x3: RiNALMo-giga, full context, on train (17,746),
      test_withheld (17,371) and id_test (1,402) -> `data/cache/giga_full_*`.
- [ ] NOHIT relaxed-scan audit (E <= 10, no GA cutoff): report both counts and
      escalate if training sequences with secondary hits to a test clan exceed
      0.5%.

## Immediate

- [ ] Check RNA-FM `position_embedding_type`; if absolute, pass the original
      `position_ids` or restrict RNA-FM to the swap conditions.
- [ ] `scripts/build_conditions.py`: conditions b (dinucleotide-shuffled distal),
      d (other-clan distal), e (pad/mask distal), local window fixed at +/-32 and
      length preserved. Condition c needs `build_homolog_pool.py`.
- [ ] Cache giga conditions b/c/d/e (full-length inputs, ~0.2 h each).
- [ ] Cache giga k=32 crops for the k-arm (train ~4 h), then {8, 128, full}.
- [ ] P2 head: train decoder B on condition a, 3 seeds; evaluate under a-e on
      id_test and withheld; record per-position distortion (1 - cos) and clan
      decodability; regress F1 on both.
- [ ] Build the canonical test_novel catalogue under Amendment 1.6: 310
      published-structure families, 2,720 sequences after the 10-per-family cap
      (above the 150-clan gate); commit its hash before scoring it.

## CPU / development

- [x] Decoder A core: nearest-neighbour DP with traceback, verified against
      exhaustive enumeration for n <= 12 (`docs/decoder_a.md`).
- [ ] Turner parity for decoder A: load `rna_turner2004.par` and match
      `RNAfold -d0` on 100 sequences. Gate for calling it a physical model.
- [ ] Withheld label calibration: F1 on the bpRNA vs seed subsets of the same
      withheld clans.
- [ ] Fold the bpRNA mirror provenance table into the Round 4 report.

## Standing

- [ ] Memorization-gap runs need UFold/MXfold2/BPfold checkpoints (the MXfold2
      repo URL 404'd earlier; find the correct one).
- [ ] structRFM pretraining manifest and overlap with every test set.
- [ ] Repository license decision (RiNALMo-giga AGPL vs MIT rest).
