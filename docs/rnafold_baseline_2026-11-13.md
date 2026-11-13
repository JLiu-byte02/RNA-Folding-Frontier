# RNAfold baseline on the v3 scoring protocol — 2026-11-13

Purpose: a clean, same-protocol thermodynamic baseline for the P2-fix heads and
for the E-C2 gate target ("beat the best clean baseline"). Nothing here is a new
method; it is a reference point.

## Protocol

- Tool: ViennaRNA 2.8.0 `RNAfold -p` (single-sequence MFE), run on the server
  (`scripts/jobs/rnafold_baseline.sbatch`, job 15568).
- Labels: `experiments/011_labels_v3/labels_v3_pairs.json.gz`.
- Scoring: `scripts/score_external_baseline.py` — `--min-true-pairs 4` applied to
  the **raw** v3 pair set, truth then reduced with `max_nested_subset`
  (Amendment 8/9 conventions). Same script and flags as every other baseline.
- Restrict set: the P2-fix `a` dump ids, so the population is identical to the
  head's "all" population (id_test 1,287; withheld 16,163).

## Results

| split | n seqs | families | family-macro F1 | 95% CI | micro F1 | TP | FP | FN | mean pred pairs |
|---|---:|---:|---:|---|---:|---:|---:|---:|---:|
| test_withheld | 16,163 | 409 | **0.6003** | [0.5733, 0.6265] | 0.5548 | 284,021 | 317,917 | 137,883 | 37.24 |
| id_test | 1,287 | 155 | **0.5955** | [0.5530, 0.6375] | 0.5295 | 27,153 | 33,291 | 14,972 | 46.97 |

Exclusions (counted, not hidden): withheld 118 empty-truth + 630 `<4` pairs;
id_test 5 + 105. Parse failures: 0. Predictions for 578 / 10 excluded
sequences fall outside the restrict set and are dropped (`excluded_predicted_pairs`
14,291 / 2,039) — they are not scored in either direction.

## Comparison with the P2-fix `a` head (same population, v3 nested)

| split | head macro | RNAfold macro | d(macro) | head micro | RNAfold micro | d(micro) |
|---|---:|---:|---:|---:|---:|---:|
| test_withheld | 0.5911 | 0.6003 | **-0.0092** | 0.5942 | 0.5548 | **+0.0394** |
| id_test | 0.5837 | 0.5955 | **-0.0118** | 0.5937 | 0.5295 | **+0.0642** |

Head numbers are the 3-seed mean of the fixed-0.5 "all" readout in
`experiments/002_e1prime_runs/p2central_fixed_summary.json`.

The head is **below plain RNAfold in family-macro F1 on both splits** and above
it in micro F1 on both. The pattern is informative rather than contradictory:
the head trades balance for bulk — it wins on large/easy families (weight in the
micro average) and loses on the many small families that the macro average
weights equally. E-C1 already showed the per-sequence FM advantage over RNAfold
is only +.0051 on average; these numbers show it does not survive family
re-weighting.

So any claim that the head "beats thermodynamics" is not supported. What is
supported is: (i) it is competitive per sequence, (ii) its advantage is
concentrated in families where the FM has coverage, which is exactly the
heterogeneity E-C1 measures.

## Consequence for E-C2

The gate must clear a *harder* bar than the head itself: withheld family-macro
0.6003, id_test 0.5955. A gate that only recovers the head's own level (~0.59)
does not pass, whatever its internal ablation says. This should be stated in the
Round 4n report so the PI can price that bar correctly.

Recorded: `experiments/002_e1prime_runs/rnafold_test_withheld_score.json`,
`experiments/002_e1prime_runs/rnafold_id_test_score.json`.
