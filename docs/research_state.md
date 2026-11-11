# Research State

Updated: 2026-11-11

## Phase

`e0-contamination-audit` — benchmark auditing and split construction, with the
E1′ apparatus being built in parallel. No model results exist yet.

## Established facts

- GitHub repository created: https://github.com/JLiu-byte02/RNA-Folding-Frontier
- Slurm partition `g-3090-2`: 4 nodes x 2 RTX 3090 GPUs (8 total, shared with
  other users). Project cap: 6 GPUs.
- Remote project root: `/home/junxuanliu/projects/RNA-Folding-Frontier` (empty).
- Login-node Python is 3.6.8; Conda modules are available. A dedicated environment
  must be created before any experiment.
- E0 first pass: ArchiveII overlaps bpRNA-1m TR0 by 20.47% exact sequence
  duplicates and 38.25% near-duplicates (k=15 containment >= 0.5). TS0, VL0 and
  bpRNA-new are sequence-clean at this threshold.
- **E0e:** against the union of BPfold's two trainers (RNAStrAlign ∪ bpRNA TR0),
  ArchiveII is **61.20% exact duplicates** and **88.40% near-duplicates**. The
  long families (16S, 23S, group II) are clean only because the trainers cap
  length at 600 nt. ArchiveII also contains 313 duplicate sequences internally.
- **E0f:** ArchiveII vs RNAStrAlign labels agree at mean pair F1 **0.9946**
  (91.2% identical), while ArchiveII vs bpRNA TR0 agree at **0.9191** (46.6%
  identical). Overlap plus label identity is the mechanism that turns leakage
  into F1.
- **E0c:** Rfam 15.0 has 4,178 families, 147 clans, 463 clan memberships. ArchiveII
  spans 8 clans; tRNA and tmRNA share CL00001, so leave-one-family-out leaks on
  that fold.
- **bpRNA-1m mirror:** the HuggingFace `rouskinlab/bpRNA-1m` mirror (66,715
  sequences, sha256 `d36c62ee...`) contains 98.7-98.9% of TR0/VL0/TS0 sequences
  and is a strict superset of their bpseq labels (0 missing pairs in 13,419
  sequences). Filtering it to canonical+GU reproduces bpseq exactly for 91.2% of
  TR0. A 90% subset (CD-HIT 4.8.1) has 26,701 representatives.
- Infrastructure: Infernal 1.1.5 built from source; Rfam 15.0 CMs pressed;
  cmscan running over 37,343 sequences (bpRNA-90 reps + ArchiveII + TS0 +
  bpRNA-new) on the login node because the CPU partition is fully allocated.
- Environment: the rff env has transformers 5.20 which requires torch >= 2.6;
  torch 2.5.1 is being replaced with 2.6.0+cu124 (driver 550.144.03 supports it).

## Open questions

- Which of the PI's three candidate directions is most novel and feasible?
- Which datasets and splits give a defensible, leakage-audited benchmark?
- What is the strongest reproducible baseline achievable under the compute cap?

## Current blockers

- CPU partition `c-80-255` is fully allocated by other users, so cmscan runs on
  the login node at low priority with 14 workers.
- `tools/mmseqs-c/bin/mmseqs` segfaults against this host's glibc; CD-HIT 4.8.1
  was built instead and is used for the 90% clustering.

## Do not claim

- Any performance number, SOTA statement, or novelty claim that is not backed by a
  committed artifact in this repository.
- That BPfold calls ArchiveII out-of-distribution: they call it a *sequence-wise*
  dataset and reserve "cross-family" for Rfam 12.3-14.10.

