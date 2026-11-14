# Artifact map

Which file backs which element of Paper A, how to regenerate it, and which CI
check keeps it honest. Everything here runs from the public repository; the only
inputs that stay on the compute server are raw corpora, model weights, and
per-sequence dumps.

| paper element | committed artifact | regenerate with | CI gate |
|---|---|---|---|
| Table T7 / Figure 7 (four-method cohort) | `experiments/002_e1prime_runs/cohort/paired_{test_withheld,id_test}.json` | `scripts/jobs/mxfold2_fold.sbatch` then `scripts/jobs/rebuild_cohort_4tools.sbatch` | `tests/test_paper_numbers.py`, `tests/test_paper_figures.py` |
| Length-band table (2.7) | `cohort/band_{split}.json` | `scripts/cohort_band_eval.py` | `tests/test_paper_numbers.py` |
| Within-family interaction (2.7) | `cohort/band_interaction_{split}.json` | `scripts/cohort_band_interaction.py` | prose only (no gate) |
| Family-size and concentration (2.7) | `cohort/family_delta_structure.json`, `family_delta_concentration.json` | `scripts/family_diagnostics.py` | `tests/test_family_diagnostics.py` |
| Seed variance (2.7, Limitation 8) | `cohort/seed_variation.json` | `scripts/family_diagnostics.py` | `tests/test_family_diagnostics.py` |
| Figure 1 (overlap, label agreement) | `experiments/000_e0_leakage/union_overlap/summary.json`, `label_discordance_*/summary.json` | `scripts/e0_union_overlap.py`, `scripts/label_discordance.py` | `tests/test_paper_figures.py` |
| Figure 2 (label projection) | `experiments/011_labels_v3/labels_v3_summary.json` | `scripts/build_labels_v3.py` | `tests/test_paper_figures.py` |
| v2/v3 cross-check (2.2) | `experiments/011_labels_v3/v2_v3_bprna_crosscheck.json` | `scripts/labels_v2_v3_crosscheck.py` | prose only |
| D1-core (2.3) | `experiments/002_e1prime_runs/d1core_{split}.json` | `scripts/jobs/p3_train_{a,block}.sbatch` + scoring scripts | prose only |
| H3 and dose-response (2.4) | `p2block_h3_{split}.json` | `scripts/jobs/p2_block_*.sbatch` | prose only |
| H2.4 / energy (2.5) | `h24_seed17_{split}_v3primary.json`, `ec2_*_score.json` | `scripts/jobs/h24_energy.sbatch`, `scripts/jobs/ec2_*.sbatch` | prose only |
| Coverage and one-hot control (2.6, Figure 6) | `ec1_p1_p2.json`, `ec1b_specificity.json` | `scripts/jobs/ec1_*.sbatch`, `ec1b_score.sbatch` | `tests/test_paper_figures.py` |

## Inputs that stay on the server

- raw corpora (`data/raw/...`), RiNALMo-giga weights, and per-sequence head dumps;
- MXfold2 checkpoint `experiments/013_mxfold2/run/model.pth`, pinned by
  `train_manifest.json` (`param_sha256`);
- Slurm logs under `logs/`.

Their hashes are recorded in the per-experiment manifests, so a reader can tell
which bytes produced each artifact without downloading them.

## What CI proves

`fresh-clone.yml` clones the repository into a new directory, verifies every
amendment hash against git blob bytes, runs the full test suite, re-renders the
figures, regenerates the family diagnostics byte-for-byte, re-parses every paper
table against its JSON, and refuses a repository that contains credentials. A
green run therefore means the committed record is self-consistent, not that the
experiments were re-run.
