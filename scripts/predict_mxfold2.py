"""Fold a split with a retrained MXfold2 model and emit ``{"id", "dbn"}`` JSONL.

MXfold2's released CLI predicts a whole FASTA at once, but has no shard logic,
no resume, and no per-row provenance.  This wrapper adds the three things the
external-baseline protocol needs:

* an explicit length cap with the excluded ids written to disk, so the common
  cohort can account for them instead of silently scoring a smaller set;
* a chunked, resumable append-only cache -- a killed job resumes where it
  stopped rather than refolding 17k sequences;
* a hard failure on any chunk that returns the wrong number of rows or a
  structure whose length disagrees with its sequence.  A silently dropped row
  would bias the family-macro average.

The folding command is the upstream CLI itself (``mxfold2 predict``), invoked
with ``@model.conf`` so the trained hyperparameters are read from the file the
training job wrote, not re-specified here.  Paths containing ``test_novel`` are
refused (Amendment 8 embargo).
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


def _parse_stdout(text: str) -> list[tuple[str, str]]:
    """Parse ``>header / SEQ / PRED (score)`` triples from mxfold2 predict."""
    lines = [ln.rstrip("\n") for ln in text.splitlines() if ln.strip()]
    if len(lines) % 3 != 0:
        raise ValueError(
            f"output has {len(lines)} non-empty lines, not a multiple of 3; "
            "a truncated record would silently drop a sequence"
        )
    out: list[tuple[str, str]] = []
    for i in range(0, len(lines), 3):
        if not lines[i].startswith(">"):
            raise ValueError(f"expected '>' header at line {i}: {lines[i][:60]!r}")
        header = lines[i][1:].strip()
        seq = lines[i + 1].strip()
        dbn = lines[i + 2].split()[0]
        if not dbn:
            raise ValueError(f"empty structure for {header!r}")
        out.append((header, dbn))
    return out


def predict_chunk(payload: tuple[list[tuple[str, str]], str, str, str, str | None, str, int, int]) -> tuple[list[tuple[str, str, str | None]], str]:
    chunk, mx_python, mx_root, param, conf, model, gpu, timeout = payload
    with tempfile.TemporaryDirectory(prefix="mx_") as tmp:
        fasta = Path(tmp) / "chunk.fa"
        with open(fasta, "w", encoding="ascii", newline="\n") as fh:
            for header, seq in chunk:
                fh.write(f">{header}\n{seq}\n")
        cmd = [mx_python, "-m", "mxfold2", "predict", str(fasta), "--param", param, "--model", model]
        if conf:
            cmd.append(f"@{conf}")
        cmd += ["--gpu", str(gpu)]
        env = dict(os.environ)
        # The vendor checkout carries the compiled ``interface`` extension; the
        # site-packages install does not.  Training used PYTHONPATH the same way.
        env["PYTHONPATH"] = mx_root + (os.pathsep + env["PYTHONPATH"] if env.get("PYTHONPATH") else "")
        try:
            proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout, env=env)
        except subprocess.TimeoutExpired:
            return [(h, "", "timeout") for h, _ in chunk], "timeout"
        if proc.returncode != 0:
            return [(h, "", (proc.stderr or proc.stdout)[-300:]) for h, _ in chunk], "nonzero"
        try:
            parsed = _parse_stdout(proc.stdout)
        except ValueError as exc:
            return [(h, "", f"parse: {exc}") for h, _ in chunk], "parse"
    if len(parsed) != len(chunk):
        return [(h, "", f"row count {len(parsed)} != {len(chunk)}") for h, _ in chunk], "count"
    results = []
    for (header, seq), (got_header, dbn) in zip(chunk, parsed):
        if len(dbn) != len(seq):
            results.append((header, "", f"dbn length {len(dbn)} != {len(seq)}"))
        else:
            results.append((header, dbn, None))
    return results, ""


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
    ap.add_argument("--param", required=True, help="model.pth written by mxfold2 train")
    ap.add_argument("--conf", default=None, help="model.conf written by mxfold2 train")
    ap.add_argument("--mx-python", default="/home/junxuanliu/.conda/envs/mxfold2/bin/python")
    ap.add_argument("--mx-root", default="/home/junxuanliu/projects/vendor/MXfold2")
    ap.add_argument("--model", default="MixC", help="folding model name; must match training")
    ap.add_argument("--gpu", type=int, default=-1)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--chunk-size", type=int, default=64)
    ap.add_argument("--timeout", type=int, default=7200)
    ap.add_argument("--max-length", type=int, default=600)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--cache", default=None)
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
    too_long: list[dict] = []
    for r in rows:
        if len(r["seq"]) > args.max_length:
            too_long.append({"id": r["id"], "length": len(r["seq"])})
        else:
            to_fold.append(r)
    if args.limit:
        to_fold = to_fold[: args.limit]
    if too_long:
        len_path = out_path.with_suffix(".excluded_length.jsonl")
        len_path.parent.mkdir(parents=True, exist_ok=True)
        with open(len_path, "w", encoding="utf-8") as fh:
            for rec in too_long:
                fh.write(json.dumps(rec, separators=(",", ":")) + "\n")
        print(f"[mxfold2] excluded {len(too_long)} rows longer than {args.max_length} nt -> {len_path}", flush=True)
    if not to_fold:
        raise SystemExit("nothing to fold after the length cap")

    done = _load_cache(cache_path)
    todo = [r for r in to_fold if r["id"] not in done]
    print(f"[mxfold2] {len(to_fold)} rows under cap; cache hit {len(done)}; running {len(todo)}", flush=True)

    chunks = [todo[i : i + args.chunk_size] for i in range(0, len(todo), args.chunk_size)]
    payloads = [([(r["id"], r["seq"].upper().replace("T", "U")) for r in c], args.mx_python, args.mx_root, args.param, args.conf, args.model, args.gpu, args.timeout) for c in chunks]
    workers = 1 if args.gpu >= 0 else max(1, args.workers)

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    n_done = 0
    with open(cache_path, "a", encoding="utf-8") as cfh, cf.ProcessPoolExecutor(max_workers=workers) as ex:
        for results, chunk_err in ex.map(predict_chunk, payloads):
            for rid, dbn, err in results:
                done[rid] = (dbn or None, err or chunk_err)
                cfh.write(json.dumps({"id": rid, "dbn": dbn or None, "err": err or chunk_err or None}, separators=(",", ":")) + "\n")
            cfh.flush()
            n_done += len(results)
            if n_done % (args.chunk_size * 10) < args.chunk_size or n_done == len(todo):
                n_err = sum(1 for r in todo if (done.get(r["id"], (None, ""))[0] is None))
                print(f"[mxfold2] {n_done}/{len(todo)} errors={n_err}", flush=True)

    errors = [(rid, v[1]) for rid, v in done.items() if v[0] is None]
    found = {rid: v[0] for rid, v in done.items() if v[0] is not None}
    if errors:
        err_path = out_path.with_suffix(".errors.jsonl")
        with open(err_path, "w", encoding="utf-8") as fh:
            for rid, err in errors:
                fh.write(json.dumps({"id": rid, "err": err}, separators=(",", ":")) + "\n")
        for rid, err in errors[:10]:
            print(f"[error] {rid}: {str(err)[:200]}", file=sys.stderr)
        print(f"[mxfold2] {len(errors)} failures -> {err_path}", flush=True)
        raise SystemExit(f"{len(errors)} sequences failed; refusing partial output (cache kept at {cache_path})")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    n_written = 0
    with open(out_path, "w", encoding="utf-8", newline="\n") as fh:
        for rid in ids:
            if rid in found:
                fh.write(json.dumps({"id": rid, "dbn": found[rid]}, separators=(",", ":")) + "\n")
                n_written += 1
    print(json.dumps({
        "tool": "mxfold2",
        "model": "MixC-retrained",
        "split": str(split),
        "param": str(args.param),
        "conf": str(args.conf) if args.conf else None,
        "max_length": args.max_length,
        "n_split": len(ids),
        "n_excluded_too_long": len(too_long),
        "n_written": n_written,
        "out": str(out_path),
        "cache": str(cache_path),
        "workers": workers,
    }))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
