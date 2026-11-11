#!/usr/bin/env python3
r"""Poll the aTrust SDP line and, when it answers, verify SSH and say what to do.

The client's "all line out" state is sticky: after the 01:21 failure it ran no
further detection, so waiting for the tray to notice a restored line is not a
recovery plan. This watcher is the cheap half of one — it only probes, logs, and
on success runs the read-only checks that must happen *before* any queued job is
submitted.

It never changes network configuration, never restarts aTrust, and never submits
anything to Slurm. Registering it with Task Scheduler is a separate, explicit
step (see `docs/network_diagnosis_2026-11-12.md`):

    schtasks /Create /SC MINUTE /MO 15 /TN RFF-tunnel-watch ^
        /TR "E:\python 3.12.9\python.exe H:\RSSP - CC\scripts\tunnel_watch.py --once"

Exit codes: 0 = line reachable and SSH works, 2 = line reachable but SSH failed,
3 = line still closed (the normal case right now).
"""

from __future__ import annotations

import argparse
import json
import socket
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LOG = ROOT / "state" / "tunnel_watch.log"
DEFAULT_TARGETS = [("218.18.235.188", 441), ("218.18.235.189", 441)]


def probe(host: str, port: int, timeout: float) -> tuple[bool, float]:
    start = time.perf_counter()
    try:
        with socket.create_connection((host, port), timeout=timeout):
            return True, round((time.perf_counter() - start) * 1000, 1)
    except OSError:
        return False, round((time.perf_counter() - start) * 1000, 1)


def log(entry: dict) -> None:
    LOG.parent.mkdir(parents=True, exist_ok=True)
    entry["ts"] = datetime.now().astimezone().isoformat(timespec="seconds")
    with open(LOG, "a") as fh:
        fh.write(json.dumps(entry) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--once", action="store_true", help="single check (default for Task Scheduler)")
    ap.add_argument("--interval", type=float, default=900, help="seconds between checks")
    ap.add_argument("--timeout", type=float, default=8.0)
    ap.add_argument("--target", action="append", default=[], help="host:port, repeatable")
    args = ap.parse_args()

    targets = (
        [tuple([t.rsplit(":", 1)[0], int(t.rsplit(":", 1)[1])]) for t in args.target]
        if args.target
        else DEFAULT_TARGETS
    )

    while True:
        results = []
        reachable = False
        for host, port in targets:
            ok, ms = probe(host, port, args.timeout)
            results.append({"target": f"{host}:{port}", "open": ok, "ms": ms})
            reachable = reachable or ok

        entry: dict = {"probe": results, "line_open": reachable}
        if reachable:
            doctor = subprocess.run(
                [sys.executable, str(ROOT / "scripts" / "remote.py"), "doctor"],
                capture_output=True,
                text=True,
                cwd=ROOT,
            )
            entry["ssh_ok"] = doctor.returncode == 0
            entry["doctor"] = (doctor.stdout or doctor.stderr)[-1500:]
            entry["next"] = (
                "line is back: run `python scripts/remote.py doctor`, then "
                "`squeue -u $USER` and sacct reconciliation before any submission; "
                "the aTrust tray may still show all-line-out until a disconnect/reconnect"
            )
        log(entry)

        if args.once or reachable:
            print(json.dumps(entry))
            if not reachable:
                return 3
            return 0 if entry.get("ssh_ok") else 2
        time.sleep(args.interval)


if __name__ == "__main__":
    raise SystemExit(main())
