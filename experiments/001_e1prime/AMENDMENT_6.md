# Amendment 6 — provenance procedure, label and pseudoknot policy, slip metrics, embargo guard

Status: written 2026-11-12 before any real evaluation number exists. It carries
the PI's Round 4f rulings. Amendment 4 is restored to the bytes that hash to
`08d6771a…`, and its former §A4.6 is reproduced here as A6.1; the episode is in
the git history, and A6.7 makes it a rule.

## A6.1 — Provenance resolution procedure (was A4.6)

Before any new computation when the cluster returns:

1. `git status` and `git log origin/main..HEAD` on the server, then push
   everything, including per-split sequence-ID membership lists and per-split
   sha256 hashes.
2. Record the modification times of the split files and of the test_novel freeze
   file, plus a statement that no scoring output exists for test_novel.
3. Determinism check: rebuild every split from committed code, seeds and the
   pinned inputs into a fresh directory and compare hashes with the on-disk
   originals. Match → the freeze is credible after the fact. Mismatch → stop and
   report before going further.

Verified locally on 2026-11-12: `experiments/001_e1prime/test_novel_eval.jsonl`
hashes to `9d1ab36c…` over 2,372 records, as declared. The 42,873-record
`test_novel` split's membership list and hash remain uncommitted.

## A6.2 — What the label pipeline actually does (measured, not assumed)

Verified by reading `scripts/build_s1prime.py`, `scripts/build_rfam_index.py`,
`src/rff/data/io.py` and `src/rff/analysis/labels.py`:

| question | bpRNA-1m labels | projected Rfam SEED labels |
|---|---|---|
| non-canonical pairs removed? | yes, `filter_pairs` keeps AU/UA/GC/CG/GU/UG only | yes, same function |
| are `<>`, `[]`, `{}` parsed? | n/a — bpRNA ships a single bracket level | yes: `pairs_from_dbn` has separate stacks per bracket type, and `project_ss` keeps all four |
| WUSS letter-pair notation (`Aa`, `Bb`, …) | n/a | **no** — `project_ss` maps every non-bracket character to `.`, so letter-encoded pairs are silently dropped |
| are crossing pairs kept? | **yes** | yes, when a bracket level encodes them |

Consequence, measured with `scripts/pseudoknot_prevalence.py` (record:
`experiments/005_label_policy/pseudoknot_prevalence.json`):

| source | sequences | with ≥1 crossing pair | pairs | crossing pairs |
|---|---:|---:|---:|---:|
| bpRNA-1m mirror (66,715 records) | 66,715 | **15,794 (23.67%)** | 5,990,889 | 369,951 (**6.18%**) |
| Rfam 15.0 seeds, ≤20 per family (30,370) | 30,370 | 5 (0.017%) | 675,574 | 46 (0.0068%) |

So crossing pairs are **not** rare in the bpRNA-derived labels. The sentence in
the split manifest, "non-canonical and pseudoknots dropped", is wrong about
pseudoknots and is recorded here as a documentation defect.

## A6.3 — Frozen scoring policy

1. **Labels.** Canonical plus GU pairs only, for every source. This is a
   scoring-time transform; split membership and hashes are unaffected.
2. **Primary F1.** Against the **maximum nested subset** of the true structure
   (`rff.analysis.labels.max_nested_subset`). Where a structure contains
   crossing pairs, take the largest non-crossing subset; ties are broken
   lexicographically on the sorted pair list, so the reduction is deterministic.
   Rationale: every decoder in our comparison (greedy nested, decoder A, and
   RNAfold) emits nested structures only.
3. **Secondary F1.** Against the full canonical set including crossing pairs, so
   a baseline that can emit them gets credit. Reported for every model on every
   split, never as the primary.
4. **Pseudoknot prevalence** is reported per split alongside both F1s, because
   the primary/secondary gap is proportional to it.
5. **Non-canonical prediction candidates.** The greedy head does not mask
   non-canonical cells today; it stays unmasked for H2.4, because that is the
   condition that was measured. The number of non-canonical predicted pairs is
   logged next to the result.
6. **Letter-pair WUSS notation** is a known gap in the builder, not a deliberate
   policy. It must be added to `project_ss` and the split rebuilt **before** any
   confirmatory number is quoted, since a rebuild changes split hashes. This is
   an item for the cluster queue, recorded rather than patched from memory.

## A6.4 — Two different "±1" diagnostics, neither primary

Round 4f caught an ambiguity in the earlier wording. They are now separate:

* **`slip_allowed`** — a truth pair counts as found when an unused prediction has
  both indices within 1. One-to-one, deterministic (nearest prediction, ties by
  lexicographic order). This metric can only be **≥** the exact F1; measured
  0.5773 vs exact 0.5407 on the ArchiveII/RNAfold control.
* **`shifted_prediction`** — the prediction is moved by ±1 and scored exactly.
  This measures how badly an index slip destroys the score (0.0055 on the same
  control) and is a sensitivity check, not a tolerance metric.

Both are reported for the scorer cross-check. Neither is used for any claim.

## A6.5 — Embargo guard for the confirmatory set

The confirmatory scripts (H2.4 and P5) refuse to read any path under
`test_novel` unless **both** hold:

1. `--confirmatory` is passed;
2. every amendment hash in `state/amendments.lock.json` matches the file on
   disk, and the `test_novel` evaluation-subset hash matches its recorded value.

Tested in both directions: refusal without the flag, and refusal when an
amendment's bytes have changed. The lockfile itself is committed, so the check
cannot be satisfied by editing an uncommitted file alone.

## A6.6 — What this does not change

No split membership, no model, no hyper-parameter. A6.2–A6.4 are scoring and
reporting rules; A6.5 is an access control. None of them has been applied to a
confirmatory result, because none exists yet.

## A6.7 — Amendments are append-only as files

Once an amendment's sha256 has been quoted to the PI, the file is frozen. Any
further change — including an addition — is a new amendment. This rule exists
because Amendment 4's hash moved after being quoted, which weakened the very
guarantee the hash was meant to provide.
