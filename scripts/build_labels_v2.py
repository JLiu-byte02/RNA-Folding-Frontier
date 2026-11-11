#!/usr/bin/env python3
"""Build the corrected Rfam seed labels (``labels_v2``) without touching splits.

PI Round 4g §2: membership is frozen, so the fix for dropped letter-pair
pseudoknots cannot be a rebuild. Instead this writes a *delta* artifact — for
every seed whose v2 pair set differs from v1, the full v2 canonical pair list.
Primary F1 keeps using v1 labels; secondary F1 on seeds can use v2. Applying the
delta to a split is then exact and reversible:

    pairs_v2 = v1_pairs if id not in delta else delta[id]

The same seed sequence appears once per family record, keyed ``<RF accession>|<seed
name>`` exactly as ``scripts/build_rfam_index.py`` keys ``seed_structures.json``.

Per-split counts need the membership lists, which are server-side (Amendment 6
§A6.1). Pass ``--split-dir`` on the cluster and this prints them; locally it
reports the global numbers and says so.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.util
import json
import re
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rff.data.io import pairs_from_dbn  # noqa: E402
from rff.data.wuss import canonical_pairs, pairs_from_wuss, project_wuss  # noqa: E402

_spec = importlib.util.spec_from_file_location(
    "build_rfam_index", ROOT / "scripts" / "build_rfam_index.py"
)
rfam_builder = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(rfam_builder)


def iter_seeds(path: Path, max_per_family: int):
    with rfam_builder.open_maybe_gzip(path) as handle:
        for meta, seqs, ss_cons in rfam_builder.parse_stockholm(handle):
            acc = meta.get("AC", "").split()[0] if meta.get("AC") else ""
            if not acc or not ss_cons:
                continue
            items = list(seqs.items())
            if max_per_family:
                items = items[:max_per_family]
            for seq_name, aligned in items:
                ungapped = re.sub(r"[-._~]", "", aligned).upper().replace("T", "U")
                if not ungapped:
                    continue
                dbn_v1 = rfam_builder.project_ss(ss_cons, aligned)
                dbn_v2 = project_wuss(ss_cons, aligned)
                if len(dbn_v1) != len(ungapped) or len(dbn_v2) != len(ungapped):
                    continue
                yield acc, seq_name, ungapped, dbn_v1, dbn_v2


def load_split_ids(split_dir: Path) -> dict[str, set[str]]:
    """Read every jsonl/txt in ``split_dir`` as a set of sequence ids."""
    out: dict[str, set[str]] = {}
    for path in sorted(split_dir.iterdir()):
        if path.suffix not in {".jsonl", ".json", ".txt", ".tsv"}:
            continue
        ids: set[str] = set()
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                if line.startswith("{"):
                    rec = json.loads(line)
                    value = rec.get("id") or rec.get("name")
                    if value:
                        ids.add(str(value))
                elif path.suffix in {".txt", ".tsv"}:
                    ids.add(line.split("\t")[0])
        if ids:
            out[path.name] = ids
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", default="data/raw/rfam/Rfam.seed.gz")
    ap.add_argument("--max-per-family", type=int, default=0)
    ap.add_argument("--split-dir", default="", help="optional: membership lists (server)")
    ap.add_argument("--out-dir", default="experiments/008_labels_v2")
    args = ap.parse_args()

    out_dir = ROOT / args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    delta: dict[str, list[list[int]]] = {}
    per_family: dict[str, dict[str, int]] = {}
    n_seqs = 0
    n_with_pairs_v1 = 0
    n_changed = 0
    pairs_v1 = pairs_v2 = 0
    letters_seen = Counter()
    n_v1_unbalanced = 0
    n_v2_unbalanced = 0
    n_v1_failed_v2_ok = 0

    for acc, name, seq, dbn_v1, dbn_v2 in iter_seeds(ROOT / args.seed, args.max_per_family):
        n_seqs += 1
        for ch in dbn_v2:
            if ch.isalpha():
                letters_seen[ch] += 1
        try:
            p1 = canonical_pairs(seq, pairs_from_dbn(dbn_v1))
        except ValueError:
            # the v1 projection dropped a bracket level the row still closes
            n_v1_unbalanced += 1
            p1 = None
        try:
            p2 = canonical_pairs(seq, pairs_from_wuss(dbn_v2))
        except ValueError:
            n_v2_unbalanced += 1
            p2 = None
        if p2 is None:
            continue
        if p1 is None:
            # build_s1prime leaves these rows in the pool with an EMPTY truth:
            # `except ValueError: pairs = []`. So v1's effective label is the
            # empty set, and v2 is a pure gain.
            n_v1_failed_v2_ok += 1
            p1 = set()
        pairs_v1 += len(p1)
        pairs_v2 += len(p2)
        if p1:
            n_with_pairs_v1 += 1
        key = f"{acc}|{name}"
        fam = per_family.setdefault(acc, {"n_sequences": 0, "n_changed": 0, "added_pairs": 0})
        fam["n_sequences"] += 1
        if p2 != p1:
            n_changed += 1
            fam["n_changed"] += 1
            fam["added_pairs"] += len(p2 - p1)
            delta[key] = [list(p) for p in sorted(p2)]

    summary = {
        "purpose": "recover WUSS letter-pair pseudoknots in Rfam seed labels without touching splits",
        "seed_source": "data/raw/rfam/Rfam.seed.gz",
        "max_per_family": args.max_per_family,
        "n_sequences": n_seqs,
        "n_sequences_with_pairs_v1": n_with_pairs_v1,
        "n_sequences_changed": n_changed,
        "pairs_v1": pairs_v1,
        "pairs_v2": pairs_v2,
        "pairs_added": pairs_v2 - pairs_v1,
        "letters_present": dict(letters_seen.most_common(10)),
        "n_v1_unbalanced_rows": n_v1_unbalanced,
        "n_v2_unbalanced_rows": n_v2_unbalanced,
        "n_v1_failed_v2_ok": n_v1_failed_v2_ok,
        "n_families": len(per_family),
        "families_changed": {k: v for k, v in per_family.items() if v["n_changed"]},
    }

    if args.split_dir:
        splits = load_split_ids(ROOT / args.split_dir)
        per_split = {}
        for name, ids in splits.items():
            touched = [k for k in ids if k in delta]
            per_split[name] = {
                "n_ids": len(ids),
                "n_ids_with_v1_and_v2_difference": len(touched),
                "added_pairs": sum(len(delta[k]) for k in touched),
            }
        summary["per_split"] = per_split
    else:
        summary["per_split"] = {
            "status": "not computed: membership lists are server-side (A6.1)",
            "what_to_run_on_the_cluster": "scripts/build_labels_v2.py --split-dir splits/s1prime",
        }

    # gzip with mtime=0 so the artifact's bytes — and therefore its hash — are
    # reproducible on any machine
    payload = json.dumps(delta, sort_keys=True).encode("utf-8")
    blob = gzip.compress(payload, mtime=0)
    delta_path = out_dir / "labels_v2_delta.json.gz"
    delta_path.write_bytes(blob)
    summary["delta_file"] = str(delta_path.relative_to(ROOT))
    summary["delta_bytes"] = delta_path.stat().st_size
    summary["delta_bytes_uncompressed"] = len(payload)
    summary["delta_sha256"] = hashlib.sha256(blob).hexdigest()
    (out_dir / "labels_v2_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({k: summary[k] for k in ("n_sequences", "n_sequences_changed",
                                              "pairs_v1", "pairs_v2", "pairs_added",
                                              "delta_bytes")}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
