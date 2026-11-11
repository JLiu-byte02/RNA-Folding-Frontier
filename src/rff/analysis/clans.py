"""Canonical family -> clan mapping, shared by every leakage analysis.

Rfam 15.0 has 4,178 families, of which only 463 belong to a clan. The release
file marks the rest ``SINGLETON``; collapsing them into one bucket would make a
leave-one-clan-out split remove 3,715 families at once, so a clanless family is
normalised to its own accession here and every script goes through this module.

The same normalisation is applied to the cmscan parser (which reads the clan
membership file directly), so a family's clan identifier is identical no matter
which table produced it.
"""

from __future__ import annotations

from pathlib import Path

SINGLETON_MARKERS = {"", "-", "SINGLETON", "NA", "none", "None"}


def normalise_clan(clan: str, accession: str) -> str:
    clan = (clan or "").strip()
    return accession if clan in SINGLETON_MARKERS else clan


def load_family_table(path: str | Path) -> dict[str, dict[str, str]]:
    """accession -> {name, clan, description, n_seed} with clans normalised."""
    path = Path(path)
    fams: dict[str, dict[str, str]] = {}
    if not path.exists():
        return fams
    lines = path.read_text().splitlines()
    if not lines:
        return fams
    header = lines[0].split("\t")
    for line in lines[1:]:
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        row = dict(zip(header, parts))
        acc = row.get("accession") or parts[0]
        fams[acc] = {
            "name": row.get("name", ""),
            "clan": normalise_clan(row.get("clan", ""), acc),
            "description": row.get("description", ""),
            "n_seed": row.get("n_seed", ""),
        }
    return fams


def load_clan_membership(path: str | Path) -> dict[str, str]:
    """clan_membership.txt (clan, accession) -> accession -> clan."""
    path = Path(path)
    clans: dict[str, str] = {}
    if not path.exists():
        return clans
    for line in path.read_text().splitlines():
        parts = line.split("\t")
        if len(parts) >= 2 and parts[0].startswith("CL"):
            clans[parts[1].strip()] = parts[0].strip()
    return clans


def clan_for(clans: dict[str, str], accession: str) -> str:
    return clans.get(accession) or accession
