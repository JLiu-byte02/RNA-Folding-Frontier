"""Frozen RNA foundation-model wrappers with a controlled receptive field.

The E1' design needs the *same* position embedded twice: once with its full
sequence as context and once with only +/-k nucleotides around it. That is a
property of the input, not of the model, so a single frozen checkpoint serves
both conditions.

Everything here is inference-only. Weights are never updated; the head that maps
embeddings to energy corrections is trained separately.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Iterable, Sequence

import torch


@dataclass
class BackboneSpec:
    name: str
    hf_id: str
    license: str
    max_positions: int


BACKBONES = {
    "rinalmo-giga": BackboneSpec("rinalmo-giga", "multimolecule/rinalmo-giga", "AGPL-3.0", 1024),
    "rinalmo-micro": BackboneSpec("rinalmo-micro", "multimolecule/rinalmo-micro", "AGPL-3.0", 1024),
    "rna-fm": BackboneSpec("rna-fm", "multimolecule/rnafm", "MIT", 1024),
}


class FrozenBackbone:
    """Wraps a HuggingFace RNA model and never computes gradients."""

    def __init__(
        self,
        spec: BackboneSpec,
        device: str = "cuda",
        dtype: torch.dtype = torch.float16,
        path: str | None = None,
    ) -> None:
        # Import the modelling submodules directly: multimolecule's top-level
        # package pulls in `datasets`, which is not installed here.
        from multimolecule.models.rinalmo import RiNALMoModel
        from multimolecule.tokenisers.rna import RnaTokenizer

        self.spec = spec
        self.device = device
        self.dtype = dtype
        local = path or os.environ.get(
            f"RFF_MODEL_{spec.name.upper().replace('-', '_')}", ""
        )
        source = local if local else spec.hf_id
        self.tokenizer = RnaTokenizer.from_pretrained(source)
        self.model = (
            RiNALMoModel.from_pretrained(source)
            .to(device)
            .to(dtype)
            .eval()
        )
        for param in self.model.parameters():
            param.requires_grad_(False)

    @torch.inference_mode()
    def embed_strings(
        self, strings: Sequence[str], batch_size: int = 96
    ) -> list[torch.Tensor]:
        """Embed each string; returns one [len(s), d] tensor per input."""
        outputs: list[torch.Tensor] = []
        for start in range(0, len(strings), batch_size):
            chunk = list(strings[start : start + batch_size])
            batch = self.tokenizer(
                chunk,
                return_tensors="pt",
                padding=True,
                add_special_tokens=False,
            ).to(self.device)
            hidden = self.model(**batch).last_hidden_state.to(torch.float32)
            mask = batch["attention_mask"].to(torch.bool)
            for row, text in enumerate(chunk):
                valid = int(mask[row].sum())
                assert valid == len(text), (
                    f"tokenizer produced {valid} tokens for a {len(text)}-nt string; "
                    "special tokens would shift the crop centre"
                )
                outputs.append(hidden[row, :valid].cpu())
        return outputs

    def embed_full(self, seqs: Sequence[str]) -> list[torch.Tensor]:
        """Full-context embeddings, truncated to the model's position limit."""
        limit = self.spec.max_positions
        return self.embed_strings([s[:limit] for s in seqs])

    def embed_crops(
        self,
        seqs: Sequence[str],
        k: int,
        positions: Iterable[int] | None = None,
    ) -> list[torch.Tensor]:
        """Per-position embeddings computed from a +/-k window only.

        For each sequence the window at position ``i`` is ``seq[max(0,i-k):i+k+1]``
        and the returned embedding is the activation at that position. Windows are
        clipped at the sequence ends (no padding tokens), so edge positions see
        shorter context, which is the same convention for every k.
        """
        limit = self.spec.max_positions
        results: list[torch.Tensor] = []
        for seq in seqs:
            seq = seq[:limit]
            idxs = list(range(len(seq))) if positions is None else list(positions)
            crops: list[str] = []
            centres: list[int] = []
            for i in idxs:
                lo = max(0, i - k)
                hi = min(len(seq), i + k + 1)
                crops.append(seq[lo:hi])
                centres.append(i - lo)
            emb = self.embed_strings(crops)
            stacked = torch.stack([e[c] for e, c in zip(emb, centres)])
            results.append(stacked)
        return results


def load_backbone(
    name: str, device: str = "cuda", path: str | None = None
) -> FrozenBackbone:
    if name not in BACKBONES:
        raise KeyError(f"unknown backbone {name!r}; choices: {sorted(BACKBONES)}")
    return FrozenBackbone(BACKBONES[name], device=device, path=path)
