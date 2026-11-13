#!/usr/bin/env python3
"""RiNALMo-giga masked pseudo-perplexity for E-C1 (Amendment 12)."""

from __future__ import annotations

import argparse
import json
import math
import zlib
from pathlib import Path

import numpy as np
import torch
import torch.nn.functional as F

ROOT = Path(__file__).resolve().parents[1]


def read_jsonl(path: Path) -> list[dict]:
    with open(path, encoding="utf-8") as fh:
        return [json.loads(line) for line in fh if line.strip()]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--families", default="experiments/002_e1prime_runs/ec1_families.json"
    )
    ap.add_argument(
        "--all-split",
        action="store_true",
        help="compute every record in --split instead of an E-C1 family map",
    )
    ap.add_argument("--split", default="splits/s1prime_filtered/test_withheld.jsonl")
    ap.add_argument("--checkpoint", default="data/models/rinalmo-giga")
    ap.add_argument("--masks", type=int, default=4)
    ap.add_argument("--mask-frac", type=float, default=0.15)
    ap.add_argument("--seeds", nargs="+", type=int, default=[17, 23, 41, 53])
    ap.add_argument("--batch-size", type=int, default=8)
    ap.add_argument("--device", default="cuda")
    ap.add_argument(
        "--out", default="experiments/002_e1prime_runs/ec1_p2_pppl.json"
    )
    args = ap.parse_args()

    from multimolecule.models.rinalmo import RiNALMoForPreTraining
    from multimolecule.tokenisers.rna import RnaTokenizer

    split_rows = read_jsonl(ROOT / args.split)
    if args.all_split:
        wanted = {rec["id"] for rec in split_rows}
        id_to_family = {rec["id"]: rec.get("family", "?") for rec in split_rows}
    else:
        families = json.loads((ROOT / args.families).read_text(encoding="utf-8"))
        wanted = {
            rid for entry in families["families"].values() for rid in entry["ids"]
        }
        id_to_family = {
            rid: family
            for family, entry in families["families"].items()
            for rid in entry["ids"]
        }
    seqs = {
        rec["id"]: rec["seq"]
        for rec in split_rows
        if rec["id"] in wanted
    }
    missing = wanted - set(seqs)
    if missing:
        raise SystemExit(f"{len(missing)} family ids missing from split")

    device = args.device if torch.cuda.is_available() else "cpu"
    tokenizer = RnaTokenizer.from_pretrained(ROOT / args.checkpoint)
    model = (
        RiNALMoForPreTraining.from_pretrained(ROOT / args.checkpoint)
        .to(device)
        .eval()
    )
    for param in model.parameters():
        param.requires_grad_(False)

    by_length: dict[int, list[str]] = {}
    for rid in sorted(wanted):
        by_length.setdefault(len(seqs[rid]), []).append(rid)

    results: dict[str, dict] = {}
    mask_id = int(tokenizer.mask_token_id)
    for length, ids in sorted(by_length.items()):
        for start in range(0, len(ids), args.batch_size):
            chunk = ids[start : start + args.batch_size]
            token_ids = [
                tokenizer(
                    seqs[rid], return_tensors="pt", add_special_tokens=False
                )["input_ids"][0]
                for rid in chunk
            ]
            for seed in args.seeds:
                inputs = []
                labels = []
                for rid, ids_tensor in zip(chunk, token_ids):
                    n = len(ids_tensor)
                    n_mask = max(1, int(round(args.mask_frac * n)))
                    rng = np.random.default_rng([seed, zlib.crc32(rid.encode())])
                    pos = torch.tensor(
                        rng.choice(n, size=n_mask, replace=False), dtype=torch.long
                    )
                    masked = ids_tensor.clone()
                    target = torch.full_like(ids_tensor, -100)
                    masked[pos] = mask_id
                    target[pos] = ids_tensor[pos]
                    inputs.append(masked)
                    labels.append(target)
                input_ids = torch.stack(inputs).to(device)
                label_ids = torch.stack(labels).to(device)
                attention = torch.ones_like(input_ids)
                with torch.no_grad():
                    logits = model(
                        input_ids=input_ids, attention_mask=attention
                    ).logits.float()
                    log_probs = F.log_softmax(logits, dim=-1)
                    nll = -log_probs.gather(-1, label_ids.clamp_min(0).unsqueeze(-1)).squeeze(-1)
                    valid = label_ids.ne(-100)
                    per_row = (nll * valid).sum(dim=1) / valid.sum(dim=1).clamp_min(1)
                for row, rid in enumerate(chunk):
                    value = float(per_row[row].item())
                    results.setdefault(
                        rid,
                        {
                            "family": id_to_family[rid],
                            "length": len(seqs[rid]),
                            "nll": [],
                        },
                    )["nll"].append(value)

    for row in results.values():
        row["pppl"] = round(
            float(np.mean([math.exp(v) for v in row["nll"]])), 6
        )
        row["nll"] = [round(float(v), 6) for v in row["nll"]]

    out = ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(
        json.dumps(
            {
                "amendment": "AMENDMENT_12.md",
                "mask_fraction": args.mask_frac,
                "masks": args.masks,
                "seeds": args.seeds,
                "n_sequences": len(results),
                "values": results,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    print(json.dumps({"out": args.out, "n_sequences": len(results)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
