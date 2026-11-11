# Round 4g reply — Claude Opus 5.5 (High), PI

Conversation: `claude.ai/chat/cc2ad40a-406e-481e-8e7b-89a64b78f2f4` (tab 959262454)
Archived: 2026-11-12, verbatim.

---

## Round 4g: good catches. Redirecting the outage work away from a new decoder

The NaN catch is exactly why the dry run was worth doing, and the slip fix and the restored Amendment 4 close those items. Two new issues from your report, then the redirect.

### 1. Check which bytes the locked hashes cover (urgent, about 10 minutes)

With `core.autocrlf=true`, a file can be CRLF in your working tree and LF in git. If the quoted hashes (Amendments 2–6, the eval subset) were computed on working-tree bytes, a Linux checkout of the same commit will hash differently. The `.gitattributes` change can also renormalise files on the next commit. So:

- Recompute every locked hash from the git blob (`git show <commit>:<path> | sha256sum`) at the commit where each hash was first quoted, and compare.
- If any differ, record both values, and from now on define every lock as a hash of the git blob bytes. The guard should verify blob bytes, not working-tree bytes.
- Run the check from a fresh clone on a Linux box (the cluster, when it returns) as well.

### 2. Labels: the pseudoknot gap matters less than it looks, but fix it without touching membership

Because Rfam marks most pseudoknots with letter pairs, dropping them makes the label sources asymmetric: bpRNA keeps crossing pairs (23.7% of sequences), while seeds keep almost none. The primary F1 is scored against the maximum nested subset for both sources, so the asymmetry mostly affects secondary F1 on seeds. Rulings:

- **Membership stays frozen.** Don't rebuild the splits.
- **Make the corrected seed labels a separate artifact.** Build `labels_v2` deterministically from the local `Rfam.seed.gz`, using the split IDs and a WUSS parser that handles letter pairs. Pin its hash in Amendment 7. Primary F1 uses the existing labels; secondary F1 on seeds uses `labels_v2`. Report how many pairs v2 adds per split.
- **Replace the ≤16-pair exact solver.** The maximum non-crossing subset of a *fixed* set of pairs is a polynomial DP: Nussinov restricted to the given arcs, about O(n²) cells. Use it with the same deterministic tie-break, and check it against your current solver on every component of 16 pairs or fewer. Also tell me what currently happens to components larger than 16.
- **Fix the manifest's "pseudoknots dropped" sentence**, with a dated correction note. It isn't a hashed file.
- **One fact to report:** were the Round 4 heads trained with crossing pairs as positives? Either is fine. It just needs to be recorded.

### 3. Don't build the stacking DP now

H2.4's confirmatory arm is A-energy-fixed (ViennaRNA with soft constraints and the fitted s, b), which is already implemented and parity-tested. Its comparator, greedy nested matching, is the production decoder. A max-weight nested DP with stacking is a third decoder. Amendment 3 classes it as exploratory (A-struct), and building it now would spend the outage on the arm that can't confirm anything.

Better use of the outage: rehearse the confirmatory path on real sequences.

1. **Stand-in scores.** Use local ArchiveII with stand-in head scores. EternaFold or CONTRAfold base-pair probabilities are best if you can install either locally. ViennaRNA probabilities are the fallback, but ViennaRNA probabilities fed back into ViennaRNA are circular, so either way this is non-protocol and a pipeline test only.
2. **Run the confirmatory chain end to end on it:** `fit_energy_calibration`, then greedy and A-energy-fixed, then `h24_decoder_compare`. Confirm the verdict basis (recall-matched, pair-budget or null) comes out, and that the guard still refuses any test_novel path.
3. **Measure the cost of the (s, b) grid.** Time ViennaRNA soft-constrained MFE per sequence by length bin. A 7×7 grid × 5 folds × about 17.7k training sequences is roughly 4M folds, which could be well over 100 CPU-h at the median length. Pre-register in Amendment 7, before any fit:
   - a coarse-to-fine grid, for example 5×5 coarse, then 5×5 around the best point;
   - a family-stratified fitting subsample (for example 3,000 sequences, fixed seed) if the full set is too slow;
   - and that the fitted (s, b) gets frozen and hashed before scoring id_test or withheld.

A-struct can be built after the confirmatory arms have run, labelled exploratory.

### Outage priorities, in order

1. Lock-hash blob check (§1).
2. Polynomial maximum-nested-subset DP (§2).
3. `labels_v2` for seeds (§2).
4. Confirmatory-path rehearsal plus grid timing, then Amendment 7 (§3).
5. A-struct, only if all of the above are done.

The cluster queue is unchanged. Lead the next report with the hash-basis result.
