#!/usr/bin/env python3
"""Extract how method papers describe ArchiveII and whether they deduplicated.

Europe PMC is used because it exposes full text as XML for open-access papers,
which makes the quotes verifiable rather than paraphrased. For each query we
print the matched sentences, so the curation step is reading, not searching.
"""

from __future__ import annotations

import argparse
import html
import json
import re
import sys
import urllib.parse
import urllib.request

SEARCH = (
    "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
    "?query={q}&format=json&resultType=core&pageSize=3"
)
FULLTEXT = (
    "https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML"
)


def fetch(url: str, timeout: int = 60) -> str:
    req = urllib.request.Request(url, headers={"User-Agent": "rff-lit-tally/0.1"})
    with urllib.request.urlopen(req, timeout=timeout) as fh:
        return fh.read().decode("utf-8", errors="replace")


def sentences(text: str, pattern: str, window: int = 2) -> list[str]:
    parts = re.split(r"(?<=[.!?])\s+", text)
    hits = []
    for idx, part in enumerate(parts):
        if re.search(pattern, part, flags=re.I):
            lo = max(0, idx - window)
            hi = min(len(parts), idx + window + 1)
            hits.append(" ".join(parts[lo:hi]))
    return hits


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("queries", nargs="+")
    ap.add_argument(
        "--pattern", default=r"ArchiveII|RNAStralign|duplicat|redundan|non-redundant"
    )
    args = ap.parse_args()

    for query in args.queries:
        print("=" * 78)
        print("QUERY:", query)
        try:
            meta = json.loads(fetch(SEARCH.format(q=urllib.parse.quote(query))))
        except Exception as exc:  # noqa: BLE001
            print("  search failed:", exc)
            continue
        results = meta.get("resultList", {}).get("result", [])
        if not results:
            print("  no results")
            continue
        hit = next((r for r in results if r.get("pmcid")), results[0])
        print(
            "  ",
            hit.get("title", "")[:120],
            "|",
            hit.get("journalTitle", ""),
            hit.get("pubYear", ""),
            "| PMCID:",
            hit.get("pmcid", "-"),
        )
        pmcid = hit.get("pmcid")
        if not pmcid:
            print("  no open full text")
            continue
        try:
            xml = fetch(FULLTEXT.format(pmcid=pmcid))
        except Exception as exc:  # noqa: BLE001
            print("  fulltext failed:", exc)
            continue
        text = re.sub(r"<[^>]+>", " ", xml)
        text = html.unescape(text)
        text = re.sub(r"\s+", " ", text)
        for sent in sentences(text, args.pattern)[:4]:
            print("  >", sent[:700])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
