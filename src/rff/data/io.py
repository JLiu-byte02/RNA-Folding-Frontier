"""Loaders for common RNA secondary-structure formats.

All loaders return :class:`Record` objects with:
- a normalized sequence (upper-case, T -> U),
- base pairs as 0-based ``(i, j)`` tuples with ``i < j``,
- optional family / clan / split metadata.
"""

from __future__ import annotations

import gzip
import json
import os
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Iterator


@dataclass
class Record:
    rid: str
    seq: str
    pairs: list[tuple[int, int]] = field(default_factory=list)
    family: str | None = None
    clan: str | None = None
    source: str | None = None
    split: str | None = None

    @property
    def length(self) -> int:
        return len(self.seq)


def normalize_seq(seq: str) -> str:
    seq = seq.strip().upper().replace("T", "U")
    return re.sub(r"[^ACGU N]", "N", seq).replace(" ", "")


def gzip_deterministic(payload: bytes) -> bytes:
    """Gzip bytes that hash identically on Windows and Linux.

    ``gzip.compress(payload, mtime=0)`` still lets zlib choose the OS byte in
    the gzip header: this machine writes ``0x0a`` (NTFS), the Linux server writes
    ``0xff`` (unknown). The compressed stream is identical, but the artifact
    sha256 is not. Pin the metadata byte so a frozen hash is reproducible across
    platforms. The value 0x0a is retained because the committed v2/v3 artifacts
    were produced with it; gzip readers treat the field as informational.
    """
    if len(payload) == 0:
        raise ValueError("refusing to write an empty artifact")
    blob = bytearray(gzip.compress(payload, mtime=0))
    blob[9] = 0x0A
    return bytes(blob)


def pairs_from_dbn(dbn: str) -> list[tuple[int, int]]:
    stacks: dict[str, list[int]] = {"(": [], "[": [], "{": [], "<": []}
    closing = {")": "(", "]": "[", "}": "{", ">": "<"}
    pairs: list[tuple[int, int]] = []
    for i, ch in enumerate(dbn.strip()):
        if ch in stacks:
            stacks[ch].append(i)
        elif ch in closing:
            opener = closing[ch]
            if not stacks[opener]:
                raise ValueError(f"unbalanced '{ch}' at position {i}")
            j = stacks[opener].pop()
            pairs.append((j, i))
    if any(stacks.values()):
        raise ValueError("unbalanced opening bracket(s)")
    return sorted(pairs)


def load_fasta_dbn(path: str | os.PathLike) -> list[Record]:
    records: list[Record] = []
    name: str | None = None
    seq_parts: list[str] = []
    dbn_parts: list[str] = []
    mode = "seq"

    def flush() -> None:
        if name is None:
            return
        seq = normalize_seq("".join(seq_parts))
        dbn = "".join(dbn_parts).strip()
        if dbn and len(dbn) != len(seq):
            raise ValueError(f"{name}: sequence and structure lengths differ")
        records.append(
            Record(rid=name, seq=seq, pairs=pairs_from_dbn(dbn) if dbn else [])
        )

    with open(path, encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith(">"):
                flush()
                name = line[1:].strip() or f"seq{len(records)}"
                seq_parts = []
                dbn_parts = []
                mode = "seq"
            elif not line.strip():
                mode = "dbn"
            elif mode == "seq":
                seq_parts.append(line.strip())
            else:
                dbn_parts.append(line.strip())
    flush()
    return records


def load_ct(path: str | os.PathLike) -> Record:
    with open(path, encoding="utf-8", errors="replace") as handle:
        lines = [line.rstrip("\n") for line in handle]
    header = lines[0].split()
    if len(header) < 2:
        raise ValueError(f"{path}: invalid CT header")
    length = int(header[0])
    bases: list[str] = []
    partner: dict[int, int] = {}
    for line in lines[1 : 1 + length]:
        fields = line.split()
        if len(fields) < 6:
            raise ValueError(f"{path}: invalid CT row: {line!r}")
        idx = int(fields[0])
        base = fields[1]
        j = int(fields[4])
        bases.append(base)
        if j > 0:
            partner[idx] = j
    pairs = sorted({tuple(sorted((i - 1, j - 1))) for i, j in partner.items()})
    rid = Path(path).stem
    return Record(rid=rid, seq=normalize_seq("".join(bases)), pairs=pairs)


def load_bpseq(path: str | os.PathLike) -> Record:
    seq_chars: list[str] = []
    pairs: list[tuple[int, int]] = []
    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            fields = line.split()
            if len(fields) < 3:
                continue
            idx = int(fields[0]) - 1
            base = fields[1]
            j = int(fields[2])
            if idx >= len(seq_chars):
                seq_chars.extend(["N"] * (idx + 1 - len(seq_chars)))
            seq_chars[idx] = base
            if j > 0:
                pairs.append(tuple(sorted((idx, j - 1))))
    return Record(
        rid=Path(path).stem,
        seq=normalize_seq("".join(seq_chars)),
        pairs=sorted(set(pairs)),
    )


def load_archiveii_json(path: str | os.PathLike) -> list[Record]:
    with open(path, encoding="utf-8") as handle:
        data = json.load(handle)
    records: list[Record] = []
    for rid, payload in data.items():
        seq = normalize_seq(payload["sequence"])
        pairs = sorted({tuple(sorted((int(i), int(j)))) for i, j in payload["structure"]})
        family = rid.split("_", 1)[0]
        records.append(
            Record(rid=rid, seq=seq, pairs=pairs, family=family, source=str(path))
        )
    return records


def load_path(path: str | os.PathLike) -> list[Record]:
    path = Path(path)
    if path.is_dir():
        records: list[Record] = []
        for file in sorted(path.rglob("*")):
            if file.suffix.lower() in {".ct", ".bpseq", ".dbn", ".fasta", ".fa"}:
                records.extend(load_path(file))
        return records
    suffix = path.suffix.lower()
    if suffix == ".json":
        return load_archiveii_json(path)
    if suffix == ".ct":
        return [load_ct(path)]
    if suffix == ".bpseq":
        return [load_bpseq(path)]
    if suffix in {".dbn", ".fasta", ".fa"}:
        return load_fasta_dbn(path)
    raise ValueError(f"unsupported format: {path}")


def iter_families(records: Iterable[Record]) -> Iterator[str]:
    for rec in records:
        yield rec.family or "unknown"

