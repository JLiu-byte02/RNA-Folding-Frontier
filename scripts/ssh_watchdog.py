"""Probe the cluster's SSH banner and append one JSON line per check.

The 2026-11-14 outage had a failure mode that a TCP check cannot see: the local
VPN interceptor accepted connections in under a millisecond and then forwarded
nothing, so `Test-NetConnection` said "open" while no SSH banner ever arrived.
This watchdog therefore classifies the *banner*, not the socket, and records
whether the connect was suspiciously fast (local interception) or slow (a real
remote peer that stayed silent).

It is deliberately passive: it never restarts the VPN client, never changes a
route, and never retries inside one run. A scheduled task can call it every few
minutes; repeated failures are logged, not hammered.

Configuration comes from the environment so the repository stays free of
internal addresses:
    RFF_HOST       host to probe (required unless --host is given)
    RFF_SSH_PORT   port, default 22
    RFF_WATCH_LOG  log path, default %LOCALAPPDATA%\\rff\\ssh_watchdog.log

Exit codes: 0 by default (the log carries the verdict), 2 with
--exit-on-down when no banner arrives, 3 on a configuration error.
"""

from __future__ import annotations

import argparse
import json
import os
import socket
import sys
import time
from datetime import datetime, timezone
from pathlib import Path


def default_log() -> Path:
    base = os.environ.get("LOCALAPPDATA") or os.environ.get("XDG_STATE_HOME") or str(Path.home())
    return Path(base) / "rff" / "ssh_watchdog.log"


def classify(connect_ms: float, banner: bytes, error: str | None) -> str:
    if banner.startswith(b"SSH-"):
        return "up"
    if connect_ms < 5.0:
        return "intercepted"
    if error:
        return "error"
    return "silent_peer"


def probe(host: str, port: int, timeout: float) -> dict:
    started = time.perf_counter()
    banner = b""
    error = None
    connect_ms = float("nan")
    try:
        with socket.create_connection((host, port), timeout=timeout) as sock:
            connect_ms = (time.perf_counter() - started) * 1000.0
            sock.settimeout(timeout)
            try:
                banner = sock.recv(255)
            except (socket.timeout, TimeoutError):
                error = "banner timeout"
            except OSError as exc:
                error = f"banner read: {exc.__class__.__name__}"
    except OSError as exc:
        error = f"connect: {exc.__class__.__name__}: {exc}"
    verdict = classify(connect_ms, banner, error)
    return {
        "at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "host": host,
        "port": port,
        "connect_ms": None if connect_ms != connect_ms else round(connect_ms, 2),
        "banner_bytes": len(banner),
        "banner_prefix": banner[:16].decode("ascii", "replace"),
        "error": error,
        "verdict": verdict,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--host", default=os.environ.get("RFF_HOST"))
    ap.add_argument("--port", type=int, default=int(os.environ.get("RFF_SSH_PORT", "22")))
    ap.add_argument("--timeout", type=float, default=8.0)
    ap.add_argument("--log", default=None)
    ap.add_argument("--quiet", action="store_true")
    ap.add_argument(
        "--exit-on-down",
        action="store_true",
        help="exit 2 when no banner arrives; default is 0 so a scheduler does "
             "not record an expected outage as a task failure",
    )
    args = ap.parse_args()

    if not args.host:
        print("RFF_HOST is not set and --host was not given", file=sys.stderr)
        return 3

    record = probe(args.host, args.port, args.timeout)
    log_path = Path(args.log) if args.log else default_log()
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with open(log_path, "a", encoding="utf-8") as fh:
        fh.write(json.dumps(record, sort_keys=True) + "\n")

    if not args.quiet:
        print(json.dumps(record, indent=1, sort_keys=True))
    if args.exit_on_down and record["verdict"] != "up":
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
