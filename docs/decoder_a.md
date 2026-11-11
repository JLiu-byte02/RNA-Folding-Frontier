# Decoder A — nearest-neighbour DP over learned scores

Status: **implemented and unit-verified 2026-11-12; Turner parity not yet done.**

Code: `src/rff/decoders/nn_dp.py`. Tests: `tests/test_nn_dp.py`. Use:
`scripts/h24_decoder_compare.py` (H2.4) and, later, BPfold-energy runs.

## What it does

Exact DP over the standard nearest-neighbour decomposition — hairpin,
interior/bulge, multiloop, exterior — with the score of each pair supplied as a
matrix. For H2.4 the matrix is the log-odds of the head's probability map, so
the decoder answers "what is the best structurally consistent structure under
these scores", not "which cells pass 0.5".

Recursion (all arrays `n x n`, `d = j - i`):

| array | meaning |
|---|---|
| `V(i,j)` | best score of a structure closed by pair `(i,j)` |
| `A0(i,j)` / `A1(i,j)` | multiloop content on `[i,j]` with `>= 0` / `>= 1` branches |
| `W(i,j)` | exterior structure on `[i,j]` |

`A0/A1` exist separately because a closing pair must see at least one branch,
while the tail after the last branch may be empty. Two bugs found during
testing are worth remembering: interior-loop offsets must have `dk, dl >= 1`
(the inner pair is strictly inside), and cells whose pair score is disallowed
must still run the `A`/`W` recursions — only `V` is skipped.

## Verification

`tests/test_nn_dp.py` enumerates **every** secondary structure for `n <= 12`
and checks that the DP optimum equals the brute-force optimum and that the
traceback reproduces it, for random pair scores with and without loop
penalties. 96 seed/penalty combinations were checked manually during
development; the committed suite covers the same contract. This validates the
recursion and traceback independently of any energy parameter — which is
precisely what a ViennaRNA comparison cannot validate.

Cost: 0.17 s at 120 nt, 1.1 s at 300 nt, 3.4 s at 500 nt (single CPU core,
sparse score matrix). id_test (1,402 sequences, median 112 nt) is minutes.

## What it is not

- **Not Turner-faithful.** No special hairpins, terminal AU/GU, Ninio/logistic
  multiloop, or dangle models. Loop terms are explicit arrays; the default used
  for H2.4 is zero penalty, i.e. a maximum-weight non-crossing matching over
  positive log-odds cells.
- **Not validated against `RNAfold`.** That check (Turner-2004 parameters,
  `-d0`, matched options, 100 sequences) is still outstanding and is the gate
  for calling decoder A a physical model. Until then it is a structural prior.

## Reading H2.4

The pre-registered test is "long-range FPs fall by more than true positives at
fixed 0.5". Note a structural fact that limits how large that effect can be: a
long-range pair only conflicts with pairs that **interleave** its endpoints.
Pairs nested inside it, or lying outside it, coexist with it, so DP often keeps
a long-range FP and adds the local TPs as well. The test is therefore an
empirical question about how much of the greedy error is genuine interleaving
conflict, not a foregone conclusion — and a null result is informative.
