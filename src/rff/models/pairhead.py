"""Pair-map decoder: frozen embeddings in, base-pair logits out.

This is decoder B of E1' (the matched-capacity learned 2D network) and the
reference implementation against which the structured DP decoder is compared.
It deliberately contains no structure prior, so any gain from the DP decoder is
attributable to the decoder rather than to extra capacity.

Pair (i, j) is described by a projection of the two residue embeddings, their
elementwise product and a learned relative-distance embedding. A stack of
dilated 3x3 convolutions then mixes local pair context. The output is
symmetrised so the model cannot depend on pair order.
"""

from __future__ import annotations

import torch
import torch.nn as nn
import torch.nn.functional as F


class ResBlock(nn.Module):
    def __init__(self, channels: int, dilation: int, dropout: float = 0.1) -> None:
        super().__init__()
        self.conv1 = nn.Conv2d(
            channels, channels, 3, padding=dilation, dilation=dilation
        )
        self.conv2 = nn.Conv2d(
            channels, channels, 3, padding=dilation, dilation=dilation
        )
        self.norm1 = nn.GroupNorm(8, channels)
        self.norm2 = nn.GroupNorm(8, channels)
        self.drop = nn.Dropout2d(dropout)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        y = F.gelu(self.norm1(self.conv1(x)))
        y = self.drop(y)
        y = self.norm2(self.conv2(y))
        return F.gelu(x + y)


class PairMapHead(nn.Module):
    def __init__(
        self,
        d_in: int = 1280,
        d: int = 32,
        n_blocks: int = 8,
        max_len: int = 1024,
        dropout: float = 0.1,
    ) -> None:
        super().__init__()
        self.proj = nn.Sequential(nn.LayerNorm(d_in), nn.Linear(d_in, d), nn.GELU())
        self.rel = nn.Embedding(max_len, d)
        nn.init.normal_(self.rel.weight, std=0.02)
        self.in_conv = nn.Conv2d(3 * d, d, 1)
        dilations = [1, 2, 4, 8] * (n_blocks // 4 + 1)
        self.blocks = nn.ModuleList(
            [ResBlock(d, dil, dropout) for dil in dilations[:n_blocks]]
        )
        self.out = nn.Conv2d(d, 1, 1)
        self.max_len = max_len

    def forward(self, h: torch.Tensor, mask: torch.Tensor | None = None) -> torch.Tensor:
        """h: [B, L, d_in] embeddings. Returns [B, L, L] logits."""
        b, n, _ = h.shape
        z = self.proj(h)
        zi = z.unsqueeze(2)  # [B, L, 1, d]
        zj = z.unsqueeze(1)  # [B, 1, L, d]
        feats = torch.cat(
            [zi.expand(b, n, n, -1), zj.expand(b, n, n, -1), zi * zj], dim=-1
        )
        idx = torch.arange(n, device=h.device)
        dist = (idx[None, :] - idx[:, None]).abs().clamp(max=self.max_len - 1)
        feats = feats + self.rel(dist).unsqueeze(0)
        x = self.in_conv(feats.permute(0, 3, 1, 2))
        for block in self.blocks:
            x = block(x)
        logits = self.out(x).squeeze(1)
        logits = 0.5 * (logits + logits.transpose(1, 2))
        if mask is not None:
            valid = (mask[:, :, None] * mask[:, None, :]).bool()
            logits = logits.masked_fill(~valid, -1e4)
        return logits


def pair_mask(lengths: torch.Tensor, n: int, min_loop: int = 3) -> torch.Tensor:
    """Boolean [B, L, L]: cells that are syntactic partners (i + min_loop < j)."""
    idx = torch.arange(n, device=lengths.device)
    valid = idx[None, :] < lengths[:, None]
    sep = (idx[None, :] - idx[:, None]) > min_loop
    return valid[:, :, None] & valid[:, None, :] & sep


def pair_bce(
    logits: torch.Tensor,
    target: torch.Tensor,
    lengths: torch.Tensor,
    pos_weight: float = 20.0,
) -> torch.Tensor:
    mask = pair_mask(lengths, logits.shape[-1], 3)
    t = target.bool()
    weight = torch.ones_like(logits) * mask
    weight = weight * torch.where(t, pos_weight, 1.0)
    loss = F.binary_cross_entropy_with_logits(
        logits, t.to(logits.dtype), weight=weight, reduction="sum"
    )
    return loss / weight.sum().clamp(min=1.0)


def pairs_to_matrix(pairs, n: int, device=None) -> torch.Tensor:
    m = torch.zeros(n, n, dtype=torch.float32, device=device)
    for i, j in pairs:
        if 0 <= i < n and 0 <= j < n:
            m[i, j] = 1.0
            m[j, i] = 1.0
    return m
