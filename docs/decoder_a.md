# Decoder A — nearest-neighbour inference over learned pair scores

Status: **A-struct implemented and enumeration-verified 2026-11-12. A-energy
verified against real ViennaRNA 2.8.0 locally (60 sequences; zero-score parity
exact, bonus attraction 60/60). Turnover calibration of `--energy-scale` still
outstanding.**

Code: `src/rff/decoders/nn_dp.py` (A-struct), `src/rff/decoders/vienna_soft.py`
(A-energy), `src/rff/training/maxmargin.py` (training). Tests:
`tests/test_nn_dp.py`, `tests/test_vienna_soft.py`, `tests/test_maxmargin.py`,
`tests/test_h24.py`. Use: `scripts/h24_decoder_compare.py --energy {struct,vienna}`.

Amendment 3 §A3.1 splits the two variants and only A-energy can confirm H2.4:
A-struct adds structural validity, which the greedy nested decoder already has;
A-energy adds *energetic* consistency over Turner terms.

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

## A-energy: Turner terms from ViennaRNA, not from us

PI Round 4c rejected re-implementing Turner 2004 (weeks of parity work, and no
claim depends on owning it). `ViennaSoftDecoder` therefore builds
`RNA.fold_compound(seq)`, adds `sc_add_bp(i, j, -scale * score)` for every
allowed cell the head scores, and takes `fc.mfe()`. With all scores zero the
decode is exactly `RNAfold`, so parity is structural.

Verified locally on ViennaRNA 2.8.0 (`scripts/check_vienna_parity.py`, record in
`experiments/003_decoder_a/parity_vienna_2.8.0.json`):

- **index base is 1-based.** The docstring's matrix overload says
  `1 <= i < j <= n`, and a 0-based call silently constrained a different pair:
  with 0-based indices 59 of 60 bonuses failed to attract the pair they were
  given. The wrapper now passes `i+1, j+1` and counts rejected constraints.
- **sign: `e = -scale * score`** attracts a pair and the decoded structure moves
  as intended (56 of 60 probes changed the MFE structure, 60/60 included the
  asked pair).
- **zero scores reproduce `RNAfold` exactly** on 60 sequences, structure and
  energy (mismatches 0).

Still open:

- a defensible `--energy-scale` (log-odds and kcal/mol are different units; the
  default 1.0 is a placeholder).

`tests/test_vienna_soft.py` pins the wiring with a stub; `tests/test_vienna_real.py`
re-runs the parity and attraction checks against the installed package (skipped
when ViennaRNA is absent).

## Training: structured max-margin

`rff.training.maxmargin` implements the MXfold2-style step: loss-augmented
decode (`+margin` on every allowed non-true cell), subgradient
`1[predicted] - 1[true]`, descent on the head's scores. Tests check the sign —
one step must raise the true-minus-predicted score and a toy must reach zero
loss — because a sign error here would silently train the complement of the
truth.

## What it is not

- **A-struct is not a physical model.** Loop terms are explicit arrays; with the
  defaults it is a maximum-weight non-crossing matching over log-odds cells, and
  it must be reported as an exploratory ablation.
- **A-energy is not yet validated.** The stub tests pin the wiring, not the
  physics; the zero-score parity check and the sign check above are the gate.

## Reading H2.4

The pre-registered test is "long-range FPs fall by more than true positives at
fixed 0.5". Note a structural fact that limits how large that effect can be: a
long-range pair only conflicts with pairs that **interleave** its endpoints.
Pairs nested inside it, or lying outside it, coexist with it, so DP often keeps
a long-range FP and adds the local TPs as well. The test is therefore an
empirical question about how much of the greedy error is genuine interleaving
conflict, not a foregone conclusion — and a null result is informative.
