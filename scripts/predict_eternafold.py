#!/usr/bin/env python3
"""Run EternaFold (CONTRAfold-SE, EternaFoldParams.v1) on a split.

Emits one ``{"id", "dbn"}`` line per sequence so the result can be scored with
``scripts/score_external_baseline.py score`` -- the same v3 label policy,
``--min-true-pairs 4`` exclusion and nested projection used for every other
baseline in this project.  The load-bearing detail is that the scorer is
shared, not reimplemented here.

Why one sequence per invocation: CONTRAfold's multi-record input is an
*alignment* reader and aborts with "Not all sequences have the same length"
when a split FASTA mixes lengths.  Single-structure prediction is the
documented path ``contrafold predict <seqfile> --params <params>``; we run it
per sequence in a process pool, which keeps the semantics exactly as
documented and still finishes a 16k split in tens of minutes on 8 cores.

Failure policy: any sequence without a parsed structure aborts the run.  A
silently dropped sequence would bias the family-macro average.  The script
refuses paths containing ``test_novel`` (embargo, Amendment 14 pending).
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


def refuse_embargoed(path: Path) -> None:
    if "test_novel" in str(path).lower():
        raise SystemExit(f"refusing {path}: test_novel is embargoed")


def read_jsonl(path: Path) -> list[dict]:
    with open(path, encoding="utf-8") as fh:
        return [json.loads(line) for line in fh if line.strip()]


def parse_structure(text: str) -> str | None:
    """Extract the dot-bracket that follows the final '>structure' tag."""
    lines = text.splitlines()
    for i, line in enumerate(lines):
        if line.strip() == ">structure":
            j = i + 1
            while j < len(lines) and not lines[j].strip():
                j += 1
            if j < len(lines):
                return lines[j].strip()
    return None


def predict_one(payload: tuple[str, str, str, str, int]) -> tuple[str, str | None, str]:
    rid, seq, contrafold, params, timeout = payload
    seq = seq.upper().replace("T", "U")
    with tempfile.TemporaryDirectory(prefix="ef_") as tmp:
        seqfile = Path(tmp) / "seq.fa"
        seqfile.write_text(f">{rid}\n{seq}\n", encoding="ascii")
        try:
            proc = subprocess.run(
                [contrafold, "predict", str(seqfile), "--params", params],
                capture_output=True,
                text=True,
                timeout=timeout,
            )
        except subprocess.TimeoutExpired:
            return rid, None, "timeout"
        if proc.returncode != 0:
            return rid, None, (proc.stderr or proc.stdout)[-300:]
        dbn = parse_structure(proc.stdout)
        if dbn is None:
            return rid, None, "no structure block"
        if len(dbn) != len(seq):
            return rid, None, f"length {len(dbn)} != {len(seq)}"
        return rid, dbn, ""


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--split", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument(
        "--contrafold", default="/home/junxuanliu/tools/EternaFold/src/contrafold"
    )
    ap.add_argument(
        "--params",
        default="/home/junxuanliu/tools/EternaFold/parameters/EternaFoldParams.v1",
    )
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--timeout", type=int, default=3600)
    args = ap.parse_args()

    split = Path(args.split)
    if not split.is_absolute():
        split = ROOT / split
    refuse_embargoed(split)
    out_path = Path(args.out)
    if not out_path.is_absolute():
        out_path = ROOT / out_path

    rows = read_jsonl(split)
    ids = [r["id"] for r in rows]
    if len(set(ids)) != len(ids):
        raise SystemExit("duplicate ids in split")

    payloads = [
        (r["id"], r["seq"], args.contrafold, args.params, args.timeout) for r in rows
    ]
    found: dict[str, str] = {}
    errors: list[tuple[str, str]] = []
    n_done = 0
    with cf.ProcessPoolExecutor(max_workers=max(1, args.workers)) as ex:
        for rid, dbn, err in ex.map(predict_one, payloads, chunksize=8):
            n_done += 1
            if dbn is None:
                errors.append((rid, err))
            else:
                found[rid] = dbn
            if n_done % 200 == 0 or n_done == len(payloads):
                print(f"[eternafold] {n_done}/{len(payloads)} errors={len(errors)}", flush=True)

    if errors:
        for rid, err in errors[:5]:
            print(f"[error] {rid}: {err}", file=sys.stderr)
        raise SystemExit(f"{len(errors)} sequences failed; refusing partial output")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as fh:
        for rid in ids:
            fh.write(
                json.dumps({"id": rid, "dbn": found[rid]}, separators=(",", ":"))
                + "\n"
            )
    print(
        json.dumps(
            {
                "tool": "eternafold",
                "split": str(split),
                "n": len(ids),
                "out": str(out_path),
                "workers": args.workers,
            }
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())