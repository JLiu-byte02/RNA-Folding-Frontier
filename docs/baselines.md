# External baselines — runner and what must be checked on the server

Code: `scripts/run_baselines.py`, `src/rff/analysis/dotbracket.py`.
Tests: `tests/test_baselines.py` (9 tests, including an end-to-end run with a
fake tool that must score 1.0).

The runner is deliberately tool-agnostic because the four baselines print
differently and one of them writes CT files:

- it takes a command template (`{fasta}`, `{outdir}` substituted);
- collects stdout plus any text file left in the output directory;
- matches predictions by `<id><TAB><db>`, `<id> <db>`, the RNAfold
  sequence-then-structure layout, or `.ct` files named after the sequence id;
- reports unmatched sequences as `n_missing` and never guesses an assignment;
  a sequence that appears twice in the split is never matched by sequence text,
  only by an explicit id.

Scoring: family-macro base-pair F1 with a family-block bootstrap
(10,000 resamples), plus sequence-micro F1 and F1 by length band.

## Turner-only reference on the local ArchiveII mirror (not a protocol number)

`scripts/archiveii_turner_control.py` runs `RNAfold` on the 3,370-sequence
mirror and scores it with the same `prf` function the model table uses. Record:
`experiments/004_turner_control/archiveii_rnafold.json`.

| quantity | value |
|---|---|
| family-macro F1 | **0.5198** |
| micro precision / recall / F1 | 0.5062 / 0.5802 / 0.5407 |
| mean pairs predicted / true per sequence | 66.31 / 57.85 |

Per family: tRNA 0.677, 23s 0.671, srp 0.626, 5s 0.607, RNaseP 0.531,
grp1 0.511, telomerase 0.460, 16s 0.414, tmRNA 0.418, grp2 0.283.

This is a **pipeline exercise on contaminated, unpacked data**: the frozen
split files are on the cluster (Amendment 4 §A4.5), canonical/GU filtering is
not applied, and ArchiveII overlaps the BPfold training union by design. It
exists to prove that the reader, the dot-bracket parser and the scorer agree,
and to have a Turner reference that does not depend on the cluster. It must not
be quoted as our RNAfold baseline.

## Commands to verify on the server (not yet run)

The exact flags below are guesses and must be confirmed with `--help` before
any number is quoted. Nothing in this section has been executed.

```bash
# RNAfold (ViennaRNA), dangle model 0 to match the decoder-A convention later
python scripts/run_baselines.py --split splits/s1prime_filtered/id_test.jsonl \
  --tool rnafold --out results/baselines/id_test/rnafold.json \
  --outdir work/rnafold --cmd "RNAfold --noPS -d0 < {fasta} > rnafold.txt"

# MXfold2
python scripts/run_baselines.py ... --tool mxfold2 --cmd "mxfold2 predict -f {fasta} -o {outdir}"

# RNAstructure: Fold writes .ct, then ct2dot; the .ct alone is enough for us
python scripts/run_baselines.py ... --tool rnastructure --cmd "Fold {fasta} {outdir}/out.ct"

# EternaFold (ships an RNAfold-compatible binary)
python scripts/run_baselines.py ... --tool eternafold --cmd "eternafold < {fasta} > eternal.txt"
```

BPfold is the important one for P5, but its interface is the least certain: it
takes a fasta and a model directory and emits its own format. Check `--help`
and, if it emits something the generic matcher cannot read, add a converter
rather than loosening the matcher.

## Fairness rules

- Same split files (`splits/s1prime_filtered/*.jsonl`), same metric function.
- Report the tool version and command line in the JSON (`run.command` is stored).
- Report `n_missing` next to every number; a baseline that silently scored 60%
  of a split is not comparable.
- Do not tune any baseline on id_test or withheld. If a tool needs parameters,
  choose them from its paper/defaults before looking at the split.

## Provenance table (Amendment 5 §A5.3)

Every baseline needs items 1–3 filled before it gets a number. Verified
2026-11-12 from package metadata, the authors' own README/LICENSE files, and the
paper text where noted. `training overlap with our clans` stays empty until
`cmscan` is reachable.

| tool | version / licence | released parameters trained on | training code / retrain on our split | clan overlap |
|---|---|---|---|---|
| **RNAfold** (ViennaRNA) | 2.8.0; ViennaRNA licence: free for research/education/commercial use, **no redistribution for a fee**, credit required | nothing sequence-specific: Turner 2004 nearest-neighbour parameters (Mathews et al., *PNAS* 101:7287, 2004) | n/a — fixed parameters, no retraining | none by construction |
| **RNAstructure** | "free of charge, GNU public license" (download page); version to be recorded on the server | Turner 2004 parameters, same as above; no sequence training | n/a | none by construction |
| **EternaFold** | v1.3.1; Stanford BSD-3-clause-style licence (LICENSE in repo) | Eterna crowdsourced high-throughput probing data (Wayment-Steele et al., *PNAS* 2022) — **not** bpRNA or RNAStrAlign | `contrafold` training code ships with the repo; retraining on our split would mean abandoning its reference parameters, so treat as fixed-parameter for P5 | unmeasured; expected small, must be checked not assumed |
| **MXfold2** | 0.1.1 (PyPI), MIT | default parameters are **TrainSetA + TrainSetB** as defined in Sato, Akiyama & Sakakibara, *Nat Mach Intell* 3:601 (2021) — **not** bpRNA-1m TR0. [CHECK: exact composition of the two sets from the paper PDF] | `mxfold2 train --model MixC --param ... data/<set>.lst` trains from BPSEQ, so retraining on our frozen train split is feasible | unmeasured; the Rivas/Rfam-era composition matters and is not the same as BPfold's |
| **BPfold** | 0.3.1 (PyPI), MIT | **RNAStrAlign (19,313 sequences, ≤600 nt) and bpRNA-1m TR0 (12,114)**, verified from the paper text in `experiments/000_e0_leakage/report.md` §E0d | training code in the repo; retraining on our split to be assessed when the cluster returns | **known to be large**: E0 measured ArchiveII at 61.20% exact / 88.40% near-duplicate against this union |

### A correction to the premise behind Amendment 5

Amendment 5 §A5.1 groups MXfold2 with BPfold as "trained on TR0 and/or
RNAStrAlign". That holds for **BPfold**, not for **MXfold2's default
parameters**, which are trained on TrainSetA/TrainSetB. The P5 consequence is
the same in direction (released weights are not trained on our frozen split), so
the rule stands, but the entry for MXfold2 must say what it actually saw rather
than inheriting BPfold's description.
