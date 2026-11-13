#!/usr/bin/env python3
"""Run RNAstructure ``Fold`` (Turner 2004 MFE) on a split.

Emits one ``{"id", "dbn"}`` line per sequence for
``scripts/score_external_baseline.py score`` -- the shared v3 label policy,
``--min-true-pairs 4`` exclusion, and nested projection are reused, not
reimplemented here.

Two environment details are load-bearing on this cluster (verified 2026-11-13):

* the shipped ``exe/Fold`` needs ``GLIBCXX_3.4.21``, which CentOS 7 lacks, so
  ``LD_LIBRARY_PATH`` must point at the conda env's ``libstdc++.so.6``;
* ``Fold`` reads its thermodynamic tables from ``DATAPATH``.

Both are set in the child environment rather than assumed from the parent
shell.  ``Fold`` accepts a one-record FASTA and writes a CT file; the CT is
parsed with the project's existing ``load_ct`` so there is exactly one CT
reader in the codebase.

Failure policy: by default any sequence without a CT, or whose CT length
disagrees with the input, aborts the whole run -- a silently dropped row would
bias the family-macro average.  ``Fold`` rejects IUPAC ambiguity codes
(``RYSWKMBDHVN``) outright, so ``--skip-ambiguous`` drops those rows *by
declared rule* and records them in ``<out>.skipped.jsonl``; the scorer then
reports them as missing predictions and the comparison table uses an explicit
common cohort (see ``experiments/002_e1prime_runs/baseline_cohort_*.jsonl``).
``--cache`` makes the run resumable: every result is appended as it completes,
so a timeout on one long sequence does not discard the rest.  Paths containing
``test_novel`` are refused (embargo).
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rff.data.io import load_ct  # noqa: E402


def refuse_embargoed(path: Path) -> None:
    if "test_novel" in str(path).lower():
        raise SystemExit(f"refusing {path}: test_novel is embargoed")


def read_jsonl(path: Path) -> list[dict]:
    with open(path, encoding="utf-8") as fh:
        return [json.loads(line) for line in fh if line.strip()]


def predict_one(payload: tuple[str, str, str, str, str, int]) -> tuple[str, str | None, str]:
    rid, seq, fold_bin, datapath, libdir, timeout = payload
    seq = seq.upper().replace("T", "U")
    env = dict(os.environ)
    env["DATAPATH"] = datapath
    env["LD_LIBRARY_PATH"] = libdir + os.pathsep + env.get("LD_LIBRARY_PATH", "")
    with tempfile.TemporaryDirectory(prefix="rs_") as tmp:
        seqfile = Path(tmp) / "seq.fa"
        ctfile = Path(tmp) / "out.ct"
        seqfile.write_text(f">{rid}\n{seq}\n", encoding="ascii")
        try:
            proc = subprocess.run(
                [fold_bin, str(seqfile), str(ctfile)],
                capture_output=True,
                text=True,
                env=env,
                timeout=timeout,
            )
        except subprocess.TimeoutExpired:
            return rid, None, "timeout"
        if not ctfile.exists():
            return rid, None, (proc.stderr or proc.stdout)[-300:]
        try:
            rec = load_ct(ctfile)
        except Exception as exc:  # noqa: BLE001 - surface the parser message
            return rid, None, f"ct parse: {exc}"
        if len(rec.seq) != len(seq):
            return rid, None, f"ct length {len(rec.seq)} != {len(seq)}"
        if rec.seq != seq:
            return rid, None, "ct sequence disagrees with input"
        chars = ["."] * len(seq)
        for i, j in rec.pairs:
            chars[i] = "("
            chars[j] = ")"
        return rid, "".join(chars), ""


def _load_cache(path: Path) -> dict[str, tuple[str | None, str]]:
    done: dict[str, tuple[str | None, str]] = {}
    if path.exists():
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                if not line.strip():
                    continue
                rec = json.loads(line)
                done[rec["id"]] = (rec.get("dbn"), rec.get("err") or "")
    return done


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--split", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--fold", default="/home/junxuanliu/tools/RNAstructure/exe/Fold")
    ap.add_argument("--datapath", default="/home/junxuanliu/tools/RNAstructure/data_tables")
    ap.add_argument("--libdir", default=os.path.expanduser("~/.conda/envs/rff/lib"))
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--timeout", type=int, default=1800)
    ap.add_argument(
        "--max-length",
        type=int,
        default=0,
        help="drop rows longer than this many nt (0 = no cap) and log them",
    )
    ap.add_argument(
        "--skip-ambiguous",
        action="store_true",
        help="drop rows containing non-ACGU residues (Fold rejects them) and log them",
    )
    ap.add_argument(
        "--cache",
        default=None,
        help="append-only per-id result cache; reused on restart (default <out>.partial.jsonl)",
    )
    ap.add_argument(
        "--allow-partial",
        action="store_true",
        help="write predictions for the sequences that succeeded even if some failed",
    )
    args = ap.parse_args()

    split = Path(args.split)
    if not split.is_absolute():
        split = ROOT / split
    refuse_embargoed(split)
    out_path = Path(args.out)
    if not out_path.is_absolute():
        out_path = ROOT / out_path
    cache_path = Path(args.cache) if args.cache else out_path.with_suffix(".partial.jsonl")
    if not cache_path.is_absolute():
        cache_path = ROOT / cache_path

    rows = read_jsonl(split)
    ids = [r["id"] for r in rows]
    if len(set(ids)) != len(ids):
        raise SystemExit("duplicate ids in split")

    to_fold: list[dict] = []
    skipped: list[dict] = []
    for r in rows:
        extra = sorted(set(r["seq"].upper()) - set("ACGU"))
        if args.skip_ambiguous and extra:
            skipped.append({"id": r["id"], "reason": "non-ACGU", "codes": extra})
        else:
            to_fold.append(r)
    if skipped:
        skip_path = out_path.with_suffix(".skipped.jsonl")
        skip_path.parent.mkdir(parents=True, exist_ok=True)
        with open(skip_path, "w", encoding="utf-8") as fh:
            for rec in skipped:
                fh.write(json.dumps(rec, separators=(",", ":")) + "\n")
        print(f"[rnastructure] skipped {len(skipped)} ambiguous rows -> {skip_path}", flush=True)

    too_long: list[dict] = []
    if args.max_length:
        kept: list[dict] = []
        for r in to_fold:
            if len(r["seq"]) > args.max_length:
                too_long.append({"id": r["id"], "length": len(r["seq"])})
            else:
                kept.append(r)
        to_fold = kept
    if too_long:
        len_path = out_path.with_suffix(".excluded_length.jsonl")
        len_path.parent.mkdir(parents=True, exist_ok=True)
        with open(len_path, "w", encoding="utf-8") as fh:
            for rec in too_long:
                fh.write(json.dumps(rec, separators=(",", ":")) + "\n")
        print(
            f"[rnastructure] excluded {len(too_long)} rows longer than "
            f"{args.max_length} nt -> {len_path}",
            flush=True,
        )

    done = _load_cache(cache_path)
    payloads = [
        (r["id"], r["seq"], args.fold, args.datapath, args.libdir, args.timeout)
        for r in to_fold
        if r["id"] not in done
    ]
    # Shortest first: ProcessPoolExecutor.map returns results in submission
    # order, so a single 4,000 nt fold at the head of the queue stalls every
    # worker behind it (job 15649 froze at 5,872/17,060 this way).
    payloads.sort(key=lambda p: len(p[1]))
    if done:
        print(f"[rnastructure] cache hit {len(done)}; running {len(payloads)}", flush=True)

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    run_ids = {p[0] for p in payloads}
    n_done = 0
    with open(cache_path, "a", encoding="utf-8") as cfh, cf.ProcessPoolExecutor(
        max_workers=max(1, args.workers)
    ) as ex:
        for rid, dbn, err in ex.map(predict_one, payloads, chunksize=4):
            done[rid] = (dbn, err)
            cfh.write(
                json.dumps({"id": rid, "dbn": dbn, "err": err or None}, separators=(",", ":"))
                + "\n"
            )
            cfh.flush()
            n_done += 1
            if n_done % 200 == 0 or n_done == len(payloads):
                n_err = sum(1 for k in run_ids if k in done and done[k][0] is None)
                print(f"[rnastructure] {n_done}/{len(payloads)} errors={n_err}", flush=True)

    errors = [(rid, dbn_or_err[1]) for rid, dbn_or_err in done.items() if dbn_or_err[0] is None]
    found = {rid: v[0] for rid, v in done.items() if v[0] is not None}
    if errors:
        err_path = out_path.with_suffix(".errors.jsonl")
        with open(err_path, "w", encoding="utf-8") as fh:
            for rid, err in errors:
                fh.write(json.dumps({"id": rid, "err": err}, separators=(",", ":")) + "\n")
        for rid, err in errors[:10]:
            print(f"[error] {rid}: {err[:200]}", file=sys.stderr)
        print(f"[rnastructure] {len(errors)} failures -> {err_path}", flush=True)
        if not args.allow_partial:
            raise SystemExit(
                f"{len(errors)} sequences failed; refusing partial output (cache kept at {cache_path})"
            )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    n_written = 0
    with open(out_path, "w", encoding="utf-8") as fh:
        for rid in ids:
            if rid in found:
                fh.write(json.dumps({"id": rid, "dbn": found[rid]}, separators=(",", ":")) + "\n")
                n_written += 1
    print(
        json.dumps(
            {
                "tool": "rnastructure",
                "split": str(split),
                "n_split": len(ids),
                "n_skipped_ambiguous": len(skipped),
                "n_written": n_written,
                "n_errors": len(errors),
                "out": str(out_path),
                "cache": str(cache_path),
                "workers": args.workers,
            }
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
