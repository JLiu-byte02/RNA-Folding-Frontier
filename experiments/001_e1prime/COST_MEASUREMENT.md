# E1' compute cost, measured

Status: measured 2026-11-11 on one RTX 3090 through `scripts/build_cache.py`
(batched, which is how the caches are actually produced). This supersedes the
single-sequence estimate in `bench_rinalmo_giga.json`: the benchmark embeds one
sequence per call and is **9x pessimistic** for full context, so it must not be
used for grid planning.

| condition | s / sequence | s / residue | 8k seqs @ L=300 |
|---|---|---|---|
| full context (L avg 184) | 0.057 | 3.1e-4 | 0.2 h |
| crop k=32 (L avg 125) | 0.766 | 6.1e-3 | 4.1 h |

Crop cost scales as `O(L * (2k+1))` residues embedded per sequence, so relative
to k=32: k=8 ~ 0.26x, k=128 ~ 4.0x, k=256 ~ 8.0x. Extrapolated grid for one
backbone over the pre-registered 8k cap:

| k | est. GPU-h (8k seqs, L<=300) |
|---|---|
| 8 | 1.1 |
| 32 | 4.1 |
| 128 | 16 |
| 256 | 32 |
| full | 0.2 |

## What this implies

1. The five credited k values cost ~54 GPU-h per backbone. On RiNALMo-giga that
   is under half the 120 GPU-h estimate but leaves no slack for the conditions,
   the decoders and retraining. The pre-registered ladder exists for this case;
   the cheap subset is {8, 32, 128, full} = 21 GPU-h on giga, ~7 GPU-h on micro.
2. **The content-substitution conditions (b/c/d/e) are full-length inputs** and
   therefore cost 0.2 GPU-h per condition per 8k sequences, i.e. 250x cheaper
   than the k=128 cache. Combined with the fact that small k also degrades the
   representation for reasons unrelated to family identity
   (`DIAGNOSTIC_context_probe.md`), the cheap conditions are both the stronger
   causal instrument and the cheaper experiment.
3. Decoder training is negligible by comparison: 6 epochs over 2.6k cached
   sequences at L<=400 took 310 s on one GPU, and the head has ~1M parameters.
