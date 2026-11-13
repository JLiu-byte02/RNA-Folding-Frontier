# SSH runbook — connection, failure, and recovery

The compute server is an internal-IP host reached over the campus/VPN path. Access
is treated as a scarce, shared resource: one login per action, no parallel
channels, no retry storms.

## 1. Before the first connection of a session

```powershell
python scripts/remote.py doctor
python scripts/remote.py status
```

`doctor` records TCP reachability, whether an SSH identification string arrives,
and a diagnosis, into `logs/ssh_doctor_<timestamp>.json`. `status` prints the
failure counter and any cooldown.

Address the three failure classes differently:

| symptom | meaning | action |
|---|---|---|
| TCP connect fails / times out | network path, VPN, or host down | fix the path; do not retry in a loop |
| TCP connects, no SSH banner | sshd-side refusal (pre-auth connection limit) or inline filter | cool down, then a single check |
| banner arrives, auth fails | credential/key problem | fix the key; never retry wrong passwords |

## 2. Connection policy

- One login per invocation. For several checks, use
  `python scripts/remote.py batch commands.txt` — a list of commands over a
  single connection.
- `keepalive` 30 s is set on every transport, so long jobs do not die of an idle
  socket.
- Retries: at most two, with 60 s and then 5 min of backoff.
- After **two consecutive failures** the helper sets `status = SSH_UNAVAILABLE`
  and refuses to connect for **30 minutes** unless `--force` is passed. This is
  deliberate: the server's pre-auth limit is only made worse by hammering it.
- Prefer an authorised key (`RFF_KEY`) over a password. The host key is pinned on
  first success and every later connection is checked against it.
- Never change server security configuration, never switch source IP to evade a
  limit, never create accounts to bypass quota.

## 3. While SSH is unavailable (SSH_UNAVAILABLE)

1. `python scripts/remote.py doctor` once, and keep the JSON report.
2. Pause **submissions**. Do not cancel or re-queue existing Slurm jobs — they
   keep running without us.
3. Continue work that does not need the cluster: literature, PI discussion,
   model/decoder design, local development, unit tests, GitHub updates.
4. Keep the pending work in `state/task_queue.json` with the exact command for
   each item, so resumption is mechanical.
5. Health-check at low frequency (>= 30 min apart), `doctor` only.

## 4. On recovery

1. `doctor` succeeds; `status` shows `OK` and a pinned host key.
2. `pwd` and `hostname` — confirm identity and working directory.
3. `squeue -u $USER` and `sacct` — reconcile Job IDs so nothing is submitted
   twice.
4. Check `logs/*.out`, results files, and checkpoints for anything that finished
   while we were away.
5. Resume from the last valid checkpoint; if a task cannot be resumed, record why
   before deciding to re-run it.
6. Sync results into the repo and report to the PI.

## 5. Known-good access

- Host `10.22.19.39`, user `junxuanliu`, port 22.
- Project root `/home/junxuanliu/projects/RNA-Folding-Frontier`.
- Env `/home/junxuanliu/.conda/envs/rff/bin/python` (3.10, torch 2.6.0+cu124).
- Slurm: `g-3090-2` (project cap 6 GPUs); CPU partition `c-80-255` is saturated
  by other users, so CPU work runs on the login node at `nice 19`.

## 6. Incident log

| when (local) | symptom | diagnosis | resolution |
|---|---|---|---|
| 2026-11-11 ~23:45 – 2026-11-12 00:30 | TCP connects, `recv` returns 0 bytes, `ssh -vvv` reports `kex_exchange_identification: read: Unknown error` | server-side pre-auth refusal, not auth and not "connection refused" | cooled down, switched to `doctor` + single checks |
| 2026-11-14 00:33 onward | TCP accepts in a few ms, zero-byte SSH banner, `kex_exchange_identification` timeout; the VPN client still reports online | local tunnel interceptor completes the handshake and forwards nothing (no proxy line); `doctor` alternates between "local interception" and "remote silent" | pending network repair; per-user watchdog installed to timestamp every probe |

## 7. Banner watchdog

`scripts/ssh_watchdog.py` reads the SSH identification string, not just the
socket, because the failure above looks "open" to a TCP check. Each run appends
one JSON line to `%LOCALAPPDATA%\rff\ssh_watchdog.log` with the connect
latency, banner size, and a verdict: `up`, `intercepted` (connect under 5 ms and
no banner), `silent_peer` (a real peer that stayed silent), or `error`.

It is passive by design: no VPN restart, no route change, no retry inside a run.
Register a per-user task (no administrator rights) with:

```
powershell -File scripts/install_ssh_watchdog.ps1 -SshHost <host> -EveryMinutes 5
powershell -File scripts/install_ssh_watchdog.ps1 -Uninstall
```

The host is supplied at install time and lives only in the local task
definition. The script exits 0 even when the host is down -- the verdict is in
the log -- so a scheduler does not record an expected outage as a task failure;
add `--exit-on-down` when a caller wants the exit code instead.
