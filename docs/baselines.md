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
