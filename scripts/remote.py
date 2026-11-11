#!/usr/bin/env python3
"""SSH/SFTP helper with a connection manager for the compute server.

Usage:
    python scripts/remote.py exec "squeue -u $USER"
    python scripts/remote.py batch commands.txt
    python scripts/remote.py put local/file remote/path
    python scripts/remote.py get remote/path local/file
    python scripts/remote.py submit path/to/job.sbatch
    python scripts/remote.py sync scripts src tests
    python scripts/remote.py doctor
    python scripts/remote.py status

Password is read from the SSHPASS environment variable. Prefer an SSH key once
one is provisioned; this helper exists so the project is never blocked on
interactive password entry.

Connection policy (see docs/runbook_ssh.md):

* one login per invocation, never parallel logins from one agent;
  ``batch`` runs a list of commands over a single connection;
* ``keepalive`` is enabled so a long command does not die of an idle transport;
* the host key is pinned on first success and verified afterwards;
* failures are counted in ``state/ssh_status.json``. After two consecutive
  failures the helper refuses to connect for 30 minutes unless ``--force`` is
  passed, so a dead server cannot be hammered by a retry loop.
"""

import argparse
import json
import os
import posixpath
import socket
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import paramiko

HOST = os.environ.get("RFF_HOST", "10.22.19.39")
USER = os.environ.get("RFF_USER", "junxuanliu")
REMOTE_ROOT = os.environ.get(
    "RFF_REMOTE_ROOT", "/home/junxuanliu/projects/RNA-Folding-Frontier"
)
ROOT = Path(__file__).resolve().parents[1]
STATE = ROOT / "state" / "ssh_status.json"
COOLDOWN_SECONDS = 30 * 60


def now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def load_state() -> dict:
    if STATE.exists():
        try:
            return json.loads(STATE.read_text())
        except json.JSONDecodeError:
            pass
    return {"status": "UNKNOWN", "failures": 0, "cooldown_until": 0}


def save_state(state: dict) -> None:
    STATE.parent.mkdir(parents=True, exist_ok=True)
    STATE.write_text(json.dumps(state, indent=2, sort_keys=True) + "\n")


def gate(state: dict, force: bool) -> None:
    if force:
        return
    until = float(state.get("cooldown_until") or 0)
    if until > time.time():
        left = int(until - time.time())
        sys.exit(
            f"SSH marked {state.get('status')} with {state.get('failures')} "
            f"consecutive failures; cooling down for {left // 60} more min. "
            f"Last error: {state.get('last_error')}. Use --force to override."
        )


def record_failure(state: dict, reason: str) -> dict:
    state["failures"] = int(state.get("failures", 0)) + 1
    state["last_error"] = reason
    state["last_attempt"] = now()
    if state["failures"] >= 2:
        state["status"] = "SSH_UNAVAILABLE"
        state["cooldown_until"] = time.time() + COOLDOWN_SECONDS
    else:
        state["status"] = "RETRY"
        state["cooldown_until"] = time.time() + 60
    save_state(state)
    return state


def pin_host_key(client: paramiko.SSHClient, state: dict) -> None:
    """Trust-on-first-use pinning, verified on every later connection.

    paramiko has no portable known_hosts entry for an arbitrary host without
    writing a file, so the pinned key is stored here and compared after the
    handshake. A mismatch aborts the session before any command runs.
    """
    key = client.get_transport().get_remote_server_key()
    fingerprint = f"{key.get_name()} {key.get_base64()}"
    known = state.get("host_fingerprint")
    if known and known != fingerprint:
        client.close()
        sys.exit(
            "host key changed: pinned %r but server offered %r. Refusing to "
            "connect; investigate before trusting this host." % (known, fingerprint)
        )
    if not known:
        state["host_fingerprint"] = fingerprint
        save_state(state)


def connect(state: dict):
    password = os.environ.get("SSHPASS")
    key_path = os.environ.get("RFF_KEY")
    client = paramiko.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    kwargs = dict(
        hostname=HOST,
        username=USER,
        timeout=30,
        banner_timeout=30,
        auth_timeout=30,
    )
    if key_path:
        kwargs["key_filename"] = key_path
    else:
        if not password:
            sys.exit("Set SSHPASS or RFF_KEY before running this script.")
        kwargs["password"] = password
    client.connect(**kwargs)
    pin_host_key(client, state)
    client.get_transport().set_keepalive(30)
    state.update({"status": "OK", "failures": 0, "cooldown_until": 0, "last_success": now()})
    save_state(state)
    return client


def resolve(path):
    if path.startswith("/"):
        return path
    return posixpath.join(REMOTE_ROOT, path)


def cmd_exec(client, command, cwd=None):
    full = f"cd {resolve(cwd)} && {command}" if cwd else command
    _, stdout, stderr = client.exec_command(full, timeout=None)
    out = stdout.read().decode("utf-8", errors="replace")
    err = stderr.read().decode("utf-8", errors="replace")
    code = stdout.channel.recv_exit_status()
    if out:
        print(out, end="")
    if err:
        print(err, end="", file=sys.stderr)
    return code


def cmd_sync(client, local_dirs, remote_root, dry_run=False):
    """Upload every file under ``local_dirs`` whose size or mtime differs.

    Skipping unchanged files matters: a full push of scripts/ and src/ is slow
    over SFTP and, more importantly, it is easy to forget one file by hand and
    then debug a stale copy on the server.
    """
    sftp = client.open_sftp()
    skipped = uploaded = 0
    try:
        for local_dir in local_dirs:
            root = os.path.abspath(local_dir)
            for dirpath, dirnames, filenames in os.walk(root):
                dirnames[:] = [
                    d
                    for d in dirnames
                    if d not in {"__pycache__", ".git", ".ipynb_checkpoints"}
                ]
                rel_dir = os.path.relpath(dirpath, root)
                for name in filenames:
                    if name.endswith((".pyc", ".pyo", ".swp")):
                        continue
                    local = os.path.join(dirpath, name)
                    rel = name if rel_dir == "." else posixpath.join(
                        rel_dir.replace(os.sep, "/"), name
                    )
                    remote = posixpath.join(remote_root, os.path.basename(root), rel)
                    try:
                        st = sftp.stat(remote)
                        same_size = st.st_size == os.path.getsize(local)
                        same_time = abs(st.st_mtime - os.path.getmtime(local)) < 1
                        if same_size and same_time:
                            skipped += 1
                            continue
                    except OSError:
                        remote_dir = posixpath.dirname(remote)
                        try:
                            sftp.stat(remote_dir)
                        except OSError:
                            sftp.mkdir(remote_dir)
                    if dry_run:
                        print(f"would upload {local} -> {remote}")
                        uploaded += 1
                        continue
                    sftp.put(local, remote)
                    sftp.utime(remote, (os.path.getatime(local), os.path.getmtime(local)))
                    uploaded += 1
                    print(f"uploaded {local} -> {remote}")
    finally:
        sftp.close()
    print(f"sync: {uploaded} uploaded, {skipped} unchanged")
    return 0


def cmd_put(client, local, remote):
    sftp = client.open_sftp()
    target = resolve(remote)
    sftp.put(local, target)
    sftp.close()
    print(f"uploaded {local} -> {target}")


def cmd_get(client, remote, local):
    sftp = client.open_sftp()
    sftp.get(resolve(remote), local)
    sftp.close()
    print(f"downloaded {resolve(remote)} -> {local}")


def cmd_batch(client, path: str) -> int:
    """Run a list of commands over one connection.

    Each non-empty, non-comment line is one command and gets a banner so the
    output stays readable. This exists so a multi-step check costs one login
    instead of ten.
    """
    code = 0
    for line in Path(path).read_text().splitlines():
        command = line.strip()
        if not command or command.startswith("#"):
            continue
        print(f"\n$ {command}", flush=True)
        code = cmd_exec(client, command, REMOTE_ROOT)
        if code != 0:
            print(f"[exit {code}]", flush=True)
    return code


def cmd_doctor() -> int:
    """Protocol-level diagnosis; writes a redacted report under logs/."""
    report: dict = {"host": HOST, "port": 22, "at": now()}
    sock = socket.socket()
    sock.settimeout(10)
    started = time.time()
    try:
        sock.connect((HOST, 22))
        report["tcp_connect"] = True
        report["tcp_connect_seconds"] = round(time.time() - started, 3)
        try:
            banner = sock.recv(256)
            report["ssh_banner_bytes"] = len(banner)
            report["ssh_banner"] = banner.decode("utf-8", "replace").strip()
        except OSError as exc:
            report["ssh_banner_bytes"] = 0
            report["ssh_banner_error"] = f"{type(exc).__name__}: {exc}"
    except OSError as exc:
        report["tcp_connect"] = False
        report["tcp_connect_error"] = f"{type(exc).__name__}: {exc}"
    finally:
        sock.close()

    # A handshake that completes essentially instantly is answered by a local
    # interceptor, not by the far end. On this host the Sangfor aTrust client
    # accepts 10.22.19.39:22 for the WFP redirect rule and then fails upstream
    # with "no proxy line", which is indistinguishable from a wedged sshd
    # unless the timing is looked at. See docs/network_diagnosis_2026-11-12.md.
    local_accept = (
        report.get("tcp_connect_seconds") is not None
        and report["tcp_connect_seconds"] < 0.005
    )
    report["suspected_local_interception"] = bool(
        report.get("tcp_connect") and local_accept
    )
    if not report.get("tcp_connect"):
        report["diagnosis"] = "tcp refused or unreachable"
    elif report.get("ssh_banner_bytes"):
        report["diagnosis"] = "sshd reachable"
    elif report["suspected_local_interception"]:
        report["diagnosis"] = (
            "local interceptor accepted the connection in <5 ms and forwarded "
            "nothing: aTrust SDP tunnel is down (upstream: no proxy line). "
            "sshd was never contacted; retrying SSH cannot help."
        )
    else:
        report["diagnosis"] = (
            "remote end accepted TCP but sent no SSH identification: sshd-side "
            "refusal, connection limit, or an inline filter"
        )
    out = ROOT / "logs" / f"ssh_doctor_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return 0 if report.get("ssh_banner_bytes") else 1


def main():
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="action", required=True)

    p_exec = sub.add_parser("exec")
    p_exec.add_argument("command")
    p_exec.add_argument("--cwd", default=None)
    p_exec.add_argument("--force", action="store_true")

    p_put = sub.add_parser("put")
    p_put.add_argument("local")
    p_put.add_argument("remote")

    p_get = sub.add_parser("get")
    p_get.add_argument("remote")
    p_get.add_argument("local")

    p_submit = sub.add_parser("submit")
    p_submit.add_argument("script")

    p_sync = sub.add_parser("sync")
    p_sync.add_argument("dirs", nargs="+")
    p_sync.add_argument("--dry-run", action="store_true")
    p_sync.add_argument("--force", action="store_true")

    p_batch = sub.add_parser("batch")
    p_batch.add_argument("path")
    p_batch.add_argument("--force", action="store_true")

    sub.add_parser("doctor")
    sub.add_parser("status")

    args = parser.parse_args()
    if args.action == "doctor":
        sys.exit(cmd_doctor())
    state = load_state()
    if args.action == "status":
        print(json.dumps(state, indent=2))
        return
    gate(state, getattr(args, "force", False))
    try:
        client = connect(state)
    except Exception as exc:  # noqa: BLE001 - any failure must be recorded
        record_failure(state, f"{type(exc).__name__}: {exc}")
        print(
            f"connect failed: {type(exc).__name__}: {exc}",
            file=sys.stderr,
        )
        print(
            f"marked {state['status']}; failures={state['failures']}; "
            f"see state/ssh_status.json",
            file=sys.stderr,
        )
        sys.exit(1)
    try:
        if args.action == "exec":
            sys.exit(cmd_exec(client, args.command, args.cwd or REMOTE_ROOT))
        if args.action == "put":
            cmd_put(client, args.local, args.remote)
        if args.action == "get":
            cmd_get(client, args.remote, args.local)
        if args.action == "submit":
            sys.exit(cmd_exec(client, f"sbatch {args.script}", REMOTE_ROOT))
        if args.action == "sync":
            sys.exit(cmd_sync(client, args.dirs, REMOTE_ROOT, args.dry_run))
        if args.action == "batch":
            sys.exit(cmd_batch(client, args.path))
    finally:
        client.close()


if __name__ == "__main__":
    main()

