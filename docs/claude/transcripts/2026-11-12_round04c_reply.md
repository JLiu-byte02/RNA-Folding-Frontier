# Round 4c — PI ruling (Claude Opus 5.5 High, 2026-11-12)

Verbatim structure of the reply, condensed to the binding content. The full
exchange is in the Claude conversation `cc2ad40a-406e-481e-8e7b-89a64b78f2f4`.

## 1. The production decoder already enforces one partner and no crossings

"Since the Round 4 numbers came from greedy nested matching, the 'e helps'
result can't be explained by multi-partner false positives. Two explanations
remain: a calibration/threshold shift (K1), or wrong pairs that still satisfy
nesting (H2). Bare vs one-partner decoding (H2.4 as written) is still worth
measuring, but it's now a minor contrast."

Restatement, to go into Amendment 3 before any dump exists:

> decoder A, with loop energies and exact DP, reduces long-range false positives
> relative to greedy nested matching. It must do so at matched recall (F1 from
> greedy at the threshold that equalises recall with decoder A), and the
> reduction must be larger than any drop in true positives. The point is that
> the test asks whether *energetic* consistency helps, not just structural
> validity, which greedy decoding already provides.

## 2. A direction conflict that D1 will settle

"H2 predicts that the model over-predicts long-range pairs. CHANRG reports the
opposite for FMs out of distribution: high precision, collapsing recall. Our
head (pos_weight 20) isn't CHANRG's, so both can be true. But if D1 on withheld
shows precision > recall, with the e gain coming mainly from recall, then H2's
spurious-partner story is wrong for OOD data, and the result is mostly the logit
shift that K1 tests." → Amendment 3 line: P > R on withheld under condition (a)
counts as evidence against H2 on OOD data.

## 3. Novelty ruling — reserve test_novel for claims that carry the paper

"I agree with your proposal."

- Confirmatory on test_novel, pre-registered: **H2.3** (P2-block rewriting
  changes long-range error structure, native head), **H2.4 as restated**,
  **P5** (decoder A with FM features vs BPfold, EternaFold, RNAstructure).
- Replication, not novelty: **H2.1**, citing Flamm et al. *if* the full text
  supports it. "I don't recall exactly what Flamm et al. show about pair counts
  vs length, so the full read is a prerequisite; record the exact claim and
  figure number in docs/literature.md. Also check whether that paper or its
  follow-ups test a learned-score-plus-exact-energy DP, since that would overlap
  H2.4."

"The original headline, that FMs fail OOD because of a family-identity shortcut,
is on track to be a negative result. If the positive control also comes back
null, it's a well-controlled negative result, which is still worth publishing.
What could carry Nature Methods now is the last of these three parts: (1) the
audit: contamination plus the label-agreement ceiling; (2) the negative
mechanism result, with a positive control; (3) a fix: exact energetic inference
over FM pair evidence that beats BPfold on clan-novel data and on S4. If (3)
fails, it's a strong NAR analysis paper, and D2 becomes the method attempt."

## 4. CPU time — do not hand-write Turner

"Skip a hand-written Turner implementation. Reproducing Turner 2004 exactly
(dangles, special hairpins, terminal AU/GU, multiloop and exterior-loop terms) is
weeks of error-prone parity work, and nothing in our claims depends on owning it.
Use ViennaRNA's soft constraints instead."

- Decoder A at inference: RNAfold MFE plus per-pair pseudo-energies from the head
  (`fold_compound.sc_add_bp`), optionally per-base unpaired terms (`sc_add_up`).
  "With all corrections set to zero it reproduces RNAfold, so parity comes for
  free."
- Training: MXfold2-style structured max-margin; loss-augmented decoding adds
  margin terms as further soft constraints; subgradient is predicted minus true
  pair indicator.
- "Your verified DP becomes the test oracle": on a simplified energy model both
  must give identical MFE structures; keep it as fallback for what ViennaRNA
  cannot express.
- [CHECK] ViennaRNA redistribution licence; `sc_add_bp` energy units and sign
  convention in the installed version.

Order while the cluster is down:

1. Amendment 3 (~30 min), commit and hash;
2. read Flamm et al. in full, record exact claims in `docs/literature.md`;
3. baselines on CPU if the frozen split files are local — RNAfold (doubles as
   the Turner-only control), EternaFold, RNAstructure, MXfold2 on id_test and
   withheld, same family-macro scorer; BPfold on CPU if not too slow; test_novel
   stays unscored;
4. decoder A on ViennaRNA soft constraints plus the max-margin loop; validate on
   ~200 TR0 sequences with random head outputs to check gradient direction, then
   a tiny real run.

## 5. Blocked VPN

"Nothing to do on our side beyond what you're doing." Queue order unchanged;
report K1's status first when the cluster returns.
