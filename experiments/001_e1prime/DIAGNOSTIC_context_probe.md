# Diagnostic: is a +/-k crop a clean distal-context intervention?

Status: measurement complete for RiNALMo-giga, 2026-11-11. Artifact:
`crop_sanity_giga.json`, produced by `scripts/crop_sanity.py` on one RTX 3090.

## Why this check exists

E1' treats "context half-width k" as a knob that removes *distal sequence
context*. That interpretation requires the frozen model's positional encoding to
be translation invariant. With absolute positional embeddings a crop also moves
the residue's absolute index, so crop-vs-full mixes "less context" with "wrong
position" and P1/P3 would be uninterpretable.

## Positional encoding (source of truth, not inference)

`data/models/rinalmo-giga/config.json` sets
`"position_embedding_type": "rotary"` and `multimolecule`'s attention applies
`RotaryEmbedding` to q/k (`modeling_rinalmo.py:856,894`). Rotary attention is a
function of *relative* offsets, so a window centred on residue i keeps every
retained residue at the same relative distance. **Crops are clean for RiNALMo.**

## Measured agreement (6 sequences, 180-400 nt, 144 positions per k)

| k | cos(crop_k(i), full(i)) | cos(crop_k(i), full(j != i)) |
|---|---|---|
| 8 | 0.405 | 0.173 |
| 32 | 0.503 | 0.200 |
| 128 | 0.827 | 0.345 |

The same position is recovered far better than a random position at every k, so
the crop embedding remains position-specific. Agreement is low at small k
because a 33-layer transformer pools over the whole visible window: at k=8 the
representation is genuinely a different object, not a mildly perturbed one.

## A test that looks informative but is not

Embedding a window alone versus the same window inside 100 nt of random flank
gives cos 0.556 (k=8), 0.749 (k=32), 0.868 (k=128). It is tempting to read this
as absolute-position leakage. It is not: attention is global, so the flanks are
*additional visible content*. This number measures content sensitivity, which is
the treatment, not a confound. It is recorded here so the naive interpretation is
not repeated.

## Consequences for the design

1. The k factor is a **receptive-field / distal-content** axis for RiNALMo and
   can be reported as such. It must be re-checked per backbone; RNA-FM ships
   learned absolute embeddings in `multimolecule` and would be confounded.
2. Because small k also degrades the representation for reasons unrelated to
   family shortcuts (the model was pretrained on full sequences), an ID-F1
   increase with k is *not* by itself evidence of the shortcut mechanism. The
   position-preserving ablations (b) dinucleotide shuffle, (c) low-identity
   homolog, (d) other-clan context, (e) pad/mask keep input length and the
   positional distribution fixed and are the stronger causal instruments.
   Recommendation to the PI: promote (b)-(e) to primary, keep k secondary.
