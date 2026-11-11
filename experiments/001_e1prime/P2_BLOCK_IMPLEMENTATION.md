# P2-block implementation notes

Implements Amendment 2 §A2.2. Written before any block cache exists; no number
below is a result.

## What is computed

`scripts/build_block_donors.py` writes a full-length donor string per record
(`{id: donor}`) for b/c/d/e. `scripts/build_block_cache.py` then, for each
16-nt block `[s, e)`, builds the composite

    donor everywhere except [s - 32, e + 32), native there

embeds it once, and keeps only positions `s..e`. Concatenating the blocks gives
a `[L, d]` embedding in which every position was produced by a pass that had at
least `min(32, distance to the sequence end)` nt of its own native flank.
Flanks are clipped at the ends, the same convention as `embed_crops`.

## Equivalences that avoid redundant compute

- **Condition (a):** donor == native, so every composite is the native string.
  The existing `giga_a_*` caches are therefore already the P2-block (a) caches;
  no re-embedding is needed.
- **Condition (e):** the composite is native within the block's flank and `N`
  elsewhere, which is a crop of half-width `block/2 + flank = 40` around each
  block rather than a fixed ±32 crop. The PREREGISTRATION calls it
  "approximately a k=32 crop"; this is the exact statement.

## Cost

`ceil(L/16)` passes per sequence: 7 at the id_test median (112 nt), 19 at
L=300. Passes are batched per sequence, so one sequence is one batch of
identical-length strings. The measured seconds/pass is written to each cache's
`manifest.json`; the 20 GPU-h figure in Amendment 2 is an estimate to be checked
against it, not a budget the run is allowed to exceed silently.

## Scope decision left open

H2.3 only needs the (a) vs (c) contrast, and (a) is free. So the first GPU work
is (c) on id_test alone (~1.4k sequences). b/d/e and the withheld split are
queued behind it and should be scoped by what H2.3 actually needs rather than
run to completion by default.

## Verification

`tests/test_block_context.py` checks that the block spans partition the
sequence, that the flank is the only native region in a composite, end
clipping, and, with a stub backbone, that each position in the assembled
embedding comes from its own native pass and not from a donor-contaminated one.
