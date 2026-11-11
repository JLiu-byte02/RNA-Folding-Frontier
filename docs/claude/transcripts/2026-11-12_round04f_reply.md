# Round 4f reply — Claude Opus 5.5 (High), PI

Conversation: `claude.ai/chat/cc2ad40a-406e-481e-8e7b-89a64b78f2f4` (tab 959262454)
Archived: 2026-11-12, verbatim.

---

## Round 4f: accepted, with three fixes before the dry run

The provenance account is accepted. The eval-subset hash being real and committed is what protects the confirmatory set, and you recorded the source of the contradiction correctly. What's still open is the full test_novel split's membership and hash; that waits for the A4.6 procedure.

### Fix 1: never edit a file after its hash has been quoted, even to append

Amendment 4's hash moved from 08d6771a… to c6471002…. A frozen amendment is a document whose hash we have quoted to each other; if appending is allowed, the hash stops proving anything.

- Restore AMENDMENT_4.md to its original bytes, so it hashes to 08d6771a… again.
- Move A4.6 into a new Amendment 6.
- Keep the git history so the episode is visible, and add a rule to `docs/decisions.md`: amendments are append-only as files; every change goes in a new amendment.
- Check that Amendment 5 wasn't edited after you reported its hash.

### Fix 2: a ±1 slip F1 of 0.0055 looks wrong

Allowing a one-position shift can only add matches, so slip-allowed F1 must be at least the exact F1 (0.54 micro here). If 0.0055 is the increment over exact, rename the field `slip_increment`. If it's the slip-allowed F1 itself, it's a bug in the slip path. The independent scorer only checked the exact path, so extend it to the slip variant too.

### Fix 3: freeze the pair and pseudoknot policy in Amendment 6, before any real evaluation

Round 1b set the label pair policy: canonical plus GU pairs, applied to all labels, with pseudoknots recorded but excluded from primary F1. The scorer you describe does something different (crossing true pairs scored as given). Please first report what the label pipeline actually does, for bpRNA labels and for projected Rfam seeds separately:

- Are non-canonical pairs removed?
- Does the WUSS parser handle `<>`, `[]`, `{}` and the letter-pair (`Aa`) pseudoknot notation?

Then freeze this:

- **Labels:** canonical plus GU pairs only, for every source. This is a scoring-time transform, so split membership and hashes are unaffected.
- **Primary F1:** against the nested true structure. Where a structure contains pseudoknots, take the maximum nested subset by pair count, ties broken deterministically. This matches Round 1b, and every decoder in our comparisons (greedy nested, decoder A, RNAfold) produces nested output only.
- **Secondary F1:** against the full canonical structure including pseudoknotted pairs, so baselines that can predict crossing pairs get credit. Also report per-split pseudoknot prevalence.
- **Non-canonical prediction candidates:** report whether the greedy head masks non-canonical cells. If it doesn't, keep it unmasked for H2.4 (it's the condition that was measured) and log how many of its predicted pairs are non-canonical.

On the baseline table: the MXfold2 correction is accepted. Its default training sets (TrainSetA and TrainSetB) are Rfam-derived, so their overlap with our withheld and novel clans should be checked by cmscan when the cluster returns. Recording the licences closes the ViennaRNA check.

### Next local task: the synthetic dry run, yes, with two additions

1. **Adversarial synthetic dumps:** empty predictions, a dense all-pairs map, NaN or inf logits, L = 1024, a family with one sequence, a family absent from predictions, and a pseudoknotted true structure. Every one should either produce a defined value or fail loudly; none should produce a silent wrong number.
2. **An embargo guard:** the confirmatory scripts (H2.4, P5) refuse any path under test_novel unless a `--confirmatory` flag is given **and** the hashes of Amendments 2–6 and the eval subset match the values recorded in a committed lockfile. Test the guard in both directions: it refuses without the flag, and it refuses when an amendment's hash doesn't match.

### Queue when the cluster returns (unchanged)

1. Provenance resolution under A6.
2. D1/D2, then K1.
3. P2-block cache.
4. A-energy-fixed: fit (s, b), then H2.4.
5. Train-and-test heads.
6. Positive control.
7. Baseline retraining plus clan-overlap checks.
8. A-energy-learned.

The next report should start with the A4 restore, the slip answer, and the implemented label policy.
