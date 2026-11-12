# Amendment 9 — scorer convention, invalidated evidence, and the corrected rerun

Status: written 2026-11-12 after the `pair_prf` defect was found in Round 4j.
It supersedes every earlier instruction that relied on the defective convention.
Amendment 8 remains in force for labels; this amendment governs scoring,
checkpoint selection, and the rerun.

## A9.1 — Frozen zero-true-positive convention

For a sequence-level pair comparison:

1. truth empty and prediction empty: `F1 = 1`;
2. zero true positives with either set non-empty: `F1 = 0`;
3. otherwise `F1 = 2PR/(P+R)`.

The defect returned `1.0` in case 2 whenever `precision + recall == 0`. That
made a completely wrong prediction on a non-empty truth look perfect.

Primary F1 also excludes sequences with fewer than four true pairs, per
Amendment 8. Consequently the empty-truth branch cannot enter the primary table.
The branch remains defined for diagnostics that explicitly include empty rows.

`test_old_zero_overlap_convention_is_rejected` in `tests/test_scoring.py` fails
under the old convention and passes under A9.1.

## A9.2 — What the defect invalidates

Every published or internally quoted family-macro number produced with
`pair_prf` is withdrawn. This includes, at minimum:

- the Round 4 P2-central table for conditions a–e, including the “P2 inverts
  the prediction” verdict, condition c being worst, and any claim that d ≥ a;
- D1 family-macro F1 and the K1 conclusions derived from it;
- the all-N and dinucleotide-shuffle template controls from the first run;
- v2/v3 cross-check mean pair-F1 from the first run;
- any checkpoint selection whose validation metric used `pair_prf`.

The Round 4 P2 verdict is withdrawn, not replaced by its opposite. H2 is
retired as a mechanism because its premise (“e helps”) was false. H2.4 remains
a standalone decoder question. H2.3 folds into P2-block.

## A9.3 — Scorer-scope audit is a gate

Before any new GPU work, `docs/corrections.md` must contain one table listing
every artifact or number that could have used a scoring function, with:

- artifact / report location;
- scoring function used;
- old value;
- corrected value, or `unaffected` with a reason.

The table must cover E0 label discordance, T4 reproduction, the Turner control,
the smoke test, P2-central, the K1 history, the template controls and the v3
cross-check. A number not in the table may not appear in a report.

## A9.4 — Native heads are retrained, not reused

The current `runs/p2_seed{17,23,41}/best.pt` checkpoints were selected with the
defective metric and are not reused for new claims.

The replacement run is condition a only, with the same frozen S1′ training
split, same cache, same architecture, same optimizer, same `pos_weight=20` and
same hyperparameters as Round 4. Seeds are 17, 23, 41. The only changed input is
the corrected validation/selection metric. Output directories are
`runs/p2fix_seed{seed}`. The retraining must start only after this amendment is
committed and its hash is recorded in `state/amendments.lock.json`.

## A9.5 — Corrected P2-central rerun

After A9.4, regenerate the native-head probability dumps for conditions a–e on
id_test and withheld. Conditions b, c and d use the existing embedding caches;
a uses the full-context cache; e is the masked control. No model is trained on
the evaluation split.

Score with:

- Rfam seed labels v3 and the maximum-nested-subset primary truth;
- fixed threshold 0.5 and oracle threshold per condition/decoder/band;
- three readouts: kept–kept pairs only (the valid causal readout), all pairs,
  and pooled micro F1; family-macro F1 is reported alongside micro everywhere.

The expected e result is much worse than a. That is a sanity check, not a
finding. The informative contrasts are b, c and d on kept–kept pairs.

H2.4 and P2-block are not launched until the PI has seen the corrected
P2-central table.

## A9.6 — Mandatory reporting invariants

Every report generated after this amendment must print family-macro and pooled
micro F1 side by side. A report fails if they differ by more than .15 without a
stated reason.

Every split must include a dinucleotide-shuffle floor. If the shuffle floor F1
is above .05, the report fails: that indicates a scoring or cache defect, not a
scientific result.

These are reporting checks. They do not add a new hypothesis or change the
pre-registration.

## A9.7 — Calibration and downstream order

The old `(s,b)` fit is cancelled because it used old checkpoints. Its output may
be kept as a record but may not be used for H2.4. After A9.4, the native-head
train logits are dumped again, `(s,b)` is fitted under Amendment 7, frozen and
hashed, and only then may H2.4 run on id_test/withheld.

P2-block follows the corrected P2-central table and the positive control. A
decision to pivot or continue D1 waits for those results.

## A9.8 — Embargo

`test_novel` remains unscored. No path under it may be read by a scoring script
without `--confirmatory` and an intact lockfile. This amendment does not change
the embargo or the membership of any split.
