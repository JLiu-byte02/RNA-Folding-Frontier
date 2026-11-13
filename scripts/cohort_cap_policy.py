"""Cohort length-cap policy: rows/families lost under candidate caps.

Read-only. Reports, per split and per cap, how many sequences exceed the cap,
how many families lose at least one member, and how many families vanish
entirely. This informs the common-cohort policy for external baselines whose
cost grows as L^3 (EternaFold, RNAstructure).
"""

import json
import sys
from collections import Counter, defaultdict

CAPS = [800, 1000, 1200, 1500, 2000, 2500, 3000]


def load(path):
    rows = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if line:
                rows.append(json.loads(line))
    return rows


def main(split_paths):
    for path in split_paths:
        rows = load(path)
        lengths = [len(r["seq"]) for r in rows]
        families = defaultdict(list)
        for row, length in zip(rows, lengths):
            families[row["family"]].append(length)
        n_rows = len(rows)
        n_fam = len(families)
        print(f"\n=== {path} ===")
        print(f"rows={n_rows} families={n_fam} "
              f"median_len={sorted(lengths)[n_rows // 2]} max_len={max(lengths)}")
        for cap in CAPS:
            dropped = [x for x in lengths if x > cap]
            fam_any = sum(1 for v in families.values() if any(x > cap for x in v))
            fam_all = sum(1 for v in families.values() if all(x > cap for x in v))
            print(f"cap={cap:5d} dropped_rows={len(dropped):6d} "
                  f"({100.0 * len(dropped) / n_rows:5.2f}%) "
                  f"fam_any_loss={fam_any:4d} ({100.0 * fam_any / n_fam:5.1f}%) "
                  f"fam_lost_all={fam_all:3d} "
                  f"rows_left={n_rows - len(dropped):6d}")


if __name__ == "__main__":
    main(sys.argv[1:])
