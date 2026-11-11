# test_novel evaluation set — frozen, and the clan gate it fails

Status: **frozen 2026-11-11 22:25 (+0800)**, before any model was scored on it.
Builder: `scripts/build_test_novel_eval.py` (seed 17). Protocol: PREREGISTRATION
Amendment 1.6.

## Files

| file | sha256 | content |
|---|---|---|
| `experiments/001_e1prime/test_novel_eval.jsonl` | `9d1ab36ca8d07a47b67b5972bf69b2be6ebd34127b777eabbc552b0706ee209b` | 2,372 sequences, 310 families |
| `experiments/001_e1prime/test_novel_catalogue.json` | (derived, same run) | per-family clan, type, consensus pairs, multiloop flag, n_in_split, n_selected |

Selection: novel clans only; SEED `#=GF SS` source = *published*; SS_cons has
≥ 4 pairs; within a family, greedy 90%-identity dedup over a seeded random order,
cap 10 sequences. `n_in_split` records how many sequences the family had, so the
cap's effect is visible per family.

## Composition

- 310 families, 2,372 sequences (mean 7.7, median 8 per family).
- Type: sRNA 186, cis-regulatory/riboswitch 108, ribozyme 14, other 2.
- 70 families contain at least one multiloop.

## The gate fails: 6 clans, not 150

Amendment 1.6 says: *"If the published-structure subset is smaller than ~150
clans, stop and re-decide rather than relaxing the rule silently."* It is 6.

The reason is structural, not a bug. A clan groups families that share a
structure; a novel family with no bpRNA support is almost always the only member
of its clan, and `scripts/parse_cmscan.py` deliberately maps singleton clans to
their own accession rather than labelling them `SINGLETON` (that invariant makes
clan identity meaningful for the bpRNA side). So 310 novel *families* collapse
into 6 multi-member clans plus 304 singletons.

Consequences, stated plainly:

- **Clan-macro F1 on test_novel is not a usable primary metric.** Bootstrapping
  by clan over 6 units gives an interval whose width is set by six numbers, and
  it would silently weight the few multi-member clans far above everything else.
- **Family-macro F1 with a family bootstrap is usable**, but it is a weaker unit
  than the pre-registered one and must be labelled as such.
- The number to compare against BPfold is also affected: BPfold reports per-family
  results, so a family-level comparison is at least like-for-like on this set.

This is escalated to the PI rather than decided here. It does not block the P2
mechanism test, which runs on `id_test` / `test_withheld` where clans are real
groups.
