# D1 and K1 — how the diagnostic is computed

Written before any dump exists, so the definitions cannot drift towards whatever
the data turns out to say. Amendments: `experiments/001_e1prime/AMENDMENT_2.md`
§A2.3 (decoder), §A2.4 (K1).

## Why three decoders

The production decoder (`pair_map_to_pairs`) is threshold + one partner +
non-crossing. Two cheaper decoders isolate what each constraint contributes:

| arm | constraints | what it measures |
|---|---|---|
| `bare` | threshold only | raw probability mass |
| `partner` | + one partner | removal of degree violations |
| `nested` | + non-crossing | the production decoder |

K1 is written against `partner`, because Amendment 2 names one-partner decoding.
`nested` is reported alongside so the K1 verdict can be checked against the
decoder that actually produced the P2 table.

## Metric

Family-macro F1, computed as: per-sequence base-pair F1 → mean within family →
mean over families. Pooling pairs instead would let long sequences dominate and
would not be comparable with the P2 table.

## Oracle threshold

One cut per (condition, decoder, length band), chosen on a fixed grid
`0.20, 0.25, ..., 0.90`. Fixed grid rather than score quantiles so two
conditions are not given different grids by their own score distributions.
`boundary: true` in the output means the argmax sat on a grid edge and the
oracle is therefore grid-limited, which is reported rather than hidden.

The oracle is deliberately generous. The head was trained with `pos_weight=20`,
so 0.5 is an arbitrary cut; letting each condition pick its best cut asks
whether the e−a gap is anything more than a scale shift.

## Commands

```
python scripts/e1prime_dump.py --checkpoint runs/p2_seed17/best.pt \
    --out-dir experiments/002_e1prime_runs/dump_seed17 \
    --eval a_id_test=... --eval e_id_test=...
python scripts/d1_scoring.py --dump-dir experiments/002_e1prime_runs/dump_seed17 \
    --out experiments/002_e1prime_runs/d1_seed17.json
python scripts/k1_check.py \
    --a experiments/002_e1prime_runs/d1_seed{17,23,41}.json.per_sequence/a_id_test.jsonl \
    --e experiments/002_e1prime_runs/d1_seed{17,23,41}.json.per_sequence/e_id_test.jsonl \
    --out experiments/002_e1prime_runs/k1_seed17_23_41.json
```

Seeds are averaged within a sequence before family means are taken, so a
sequence contributes one number. The bootstrap is paired over families with
10,000 draws, the same construction as the P2 table.

## Verdict rule (from Amendment 2, not from this document)

K1 fires — and the "distal context is a net-negative prior" reading is dropped —
iff `partner`-arm e−a in the 80–299 band is `< 0.01` with a 95% interval that
includes 0. `scripts/k1_check.py` prints exactly that predicate and nothing
softer.

## What would make D1 uninformative

- The dump floor (0.02) truncating the candidate set at an oracle cut below the
  floor. The oracle grid bottoms out at 0.20, well above the floor, so this
  cannot bite unless the grid edge is hit; boundary flags are reported.
- Fewer than ~20 families in the 80–299 band, which would make the family-macro
  number and its bootstrap unstable. `n` is reported per band.
