#!/usr/bin/env python3
"""Minimal SSH/SFTP helper for the RNA-Folding-Frontier compute server.

Usage:
    python scripts/remote.py exec "squeue -u $USER"
    python scripts/remote.py put local/file remote/path
    python scripts/remote.py get remote/path local/file
    python scripts/remote.py submit path/to/job.sbatch

Password is read from the SSHPASS environment variable. Prefer an SSH key once
one is provisioned; this helper exists so the project is never blocked on
interactive password entry.
"""

import argparse
import os
import posixpath
import sys

import paramiko

HOST = os.environ.get("RFF_HOST", "10.22.19.39")
USER = os.environ.get("RFF_USER", "junxuanliu")
REMOTE_ROOT = os.environ.get(
    "RFF_REMOTE_ROOT", "/home/junxuanliu/projects/RNA-Folding-Frontier"
)


def connect():
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


def main():
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="action", required=True)

    p_exec = sub.add_parser("exec")
    p_exec.add_argument("command")
    p_exec.add_argument("--cwd", default=None)

    p_put = sub.add_parser("put")
    p_put.add_argument("local")
    p_put.add_argument("remote")

    p_get = sub.add_parser("get")
    p_get.add_argument("remote")
    p_get.add_argument("local")

    p_submit = sub.add_parser("submit")
    p_submit.add_argument("script")

    args = parser.parse_args()
    client = connect()
    try:
        if args.action == "exec":
            sys.exit(cmd_exec(client, args.command, args.cwd))
        if args.action == "put":
            cmd_put(client, args.local, args.remote)
        if args.action == "get":
            cmd_get(client, args.remote, args.local)
        if args.action == "submit":
            sys.exit(cmd_exec(client, f"sbatch {args.script}", REMOTE_ROOT))
    finally:
        client.close()


if __name__ == "__main__":
    main()

