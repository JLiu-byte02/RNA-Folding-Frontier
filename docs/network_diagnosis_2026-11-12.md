# Why SSH to 10.22.19.39:22 fails (2026-11-12, ~00:00 onward)

Status: **resolved to a cause, not to a working connection.** The break is on the
school's side of the aTrust SDP tunnel; nothing on this laptop can restore it
without a fresh login that re-fetches the line list.

## Symptom

TCP to `10.22.19.39:22` succeeds in about 0 ms, then zero bytes arrive: no SSH
identification string, then the connection closes.

## What it is not

| hypothesis | evidence | verdict |
|---|---|---|
| sshd is down or refusing | packets never leave the laptop's aTrust tunnel; see below | **no** |
| wrong route (traffic escaping to the internet) | see the aTrust redirect evidence; `Find-NetRoute` is misleading here because the redirect happens in the Windows Filtering Platform, below routing | **no** |
| wrong password / auth | no banner is exchanged at all; auth is never reached | **no** |
| IP ban / fail2ban | would still show a banner or a reset from the real host | **no** |
| ZoogVPN breaking aTrust | the failing leg is the school's WAN line, and third-party vantage points fail identically (below) | **no** |

## Evidence

1. **aTrust owns the connection.** `aTrustTray`'s `xtunnel/tcpAccess.log` records
   every attempt as `flowSrc: wfpRedirect`, `appName: 唐继军老师服务器-10.22.19.39`,
   `resource: 10.22.19.39:22`, `use0RTT: true`. So the 0 ms accept is the local
   proxy, not the server.

2. **The upstream leg fails with no line.**

   ```
   "upstream":{"ec":1,"ecat":"proxyErrorCategory","ecMsg":"no proxy line",
               "tcpConnTime":9223372036854775807, ...}
   ```

   `9223372036854775807` is `int64` max: "never connected". Zero bytes in both
   directions, which is exactly the missing banner.

3. **The line itself died at 23:59.** `xtunnel/xtunnel.log`:

   ```
   18:34:37  choose current line 218.18.235.189:441 ... delayMS 46
   21:42:50  ... delayMS 1055          (working, but degrading)
   23:59:29  ... error: dial tcp 218.18.235.189:441: i/o timeout
   23:59:32  enter all line out status, nodeGroup Default state success
   ```

   The node group has exactly one WAN line. When it times out, every resource on
   that node group is dead. Two live sessions (opened 22:40 and 23:50) were then
   torn down at 00:01:18, which is when the break became visible.

4. **The school's line port is closed to the whole internet, not just to us.**
   `check-host.net` TCP probes to `218.18.235.189:441` from Indonesia, Iceland
   and two US nodes: **all "Connection timed out"**. The same four-node test on
   `218.18.235.189:443` (the aTrust portal): **all succeed, 236-288 ms**. The host
   is up; port 441 is the part that is gone.

5. **Local corroboration.** From this laptop, `:441` and `:444` on that host time
   out (~15 s), while `:80` and `:443` connect. `aTrustClientTools.exe` refuses to
   run for a non-Sangfor parent, so the client cannot be driven from the shell.

## Where the workaround has to come from

The client caches the line list obtained at login (18:34). If the school moved
the SDP line to another port or address, only a re-login fetches the new list;
the current process will keep dialling the dead `:441` line.

So, in order of preference:

1. **Re-login to aTrust** (tray menu -> 退出登录, then log in). This re-fetches
   the line list. Needs the account password and possibly the SMS factor, which
   is why it has not been done automatically.
2. If re-login still shows all lines out, the SDP gateway is down or its line
   port changed server-side; ask the school's network administrator whether the
   aTrust line endpoint is still `218.18.235.189:441`.
3. Anything that changes the network path (campus Wi-Fi, phone hotspot) is worth
   one test, because it also refreshes the line negotiation.

## What not to do

## Update 02:12 — the re-login did fetch new lines, and one is new

The 01:20 re-login **did** refresh the line list, which corrects the earlier
"controller publishes zero lines" reading. At 01:21:30 the config reload carries
three lines:

```
10.20.4.21:441          (lan)
vpn.suat-sz.edu.cn:441  (wan)  -> 218.18.235.188
218.18.235.189:441      (wan)  second address, same node group
```

So the school has a second WAN address (`218.18.235.188`, published under the
hostname `vpn.suat-sz.edu.cn`) that the earlier diagnosis never tested. It is
usable in principle: DNS resolves the hostname to `.188`, and `.188:443`
answers (as does `.189:443`). The SDP line port, however, is closed on both:

| probe (01:35–02:12) | result |
|---|---|
| `218.18.235.188:441` | closed (15 s timeout) |
| `218.18.235.189:441` | closed (15 s timeout) |
| `218.18.235.188:443` | open |
| `218.18.235.189:443` | open |
| `https://vpn.suat-sz.edu.cn/` | HTTP 200, serves the aTrust login page |

That combination is the useful part: the gateway is *healthy* (portal 200 on both
addresses), and only the line service on `:441` is unreachable from the public
internet on both addresses. A dead appliance, an expired certificate or a
blocked client would not leave the portal answering.

The client then entered `all line out` at 01:21:36, six seconds after the reload,
and has run no further detection since: the log's last selectLine event is that
transition, and detection in the healthy period was triggered by network change
or latency degradation, not by a free-running timer. So the "all line out" state
is sticky — when the school restores the service, a tray disconnect/reconnect
(not necessarily a full re-login) is needed to force a fresh detection round.

**Question for the school's network administrator** (now specific enough to
answer in one reply): is the aTrust SDP line port still 441, is `218.18.235.188`
the current line address, and is there a maintenance window in progress? The
portal is up and both candidate addresses answer on 443, so the line port is the
only thing to explain.

- Do not keep retrying SSH: each attempt is accepted locally and dropped, so it
  cannot succeed and only adds noise.
- Do not change default routes or delete the ZoogVPN default route. The redirect
  happens in WFP; routing changes will not move this traffic and could break the
  Claude path.
- Do not disable the aTrust VNIC or the aTrust service. The tunnel is down, but
  reinstalling the data plane would not create a line that the WAN cannot reach.

## After it is restored

`python scripts/remote.py doctor` now distinguishes the two failure modes:
`suspected_local_interception` is true when the handshake completes in under
5 ms, which means a local proxy accepted it and the far end was never contacted.

## Update 01:35 — after a manual aTrust logout/login

The user logged out and back in at about 01:20. Nothing improved, and the state
changed in one informative way:

| | before re-login | after re-login |
|---|---|---|
| node group `Default` | one line, `218.18.235.189:441`, `delayMS 46` -> `1055` -> `i/o timeout` | `{"state":"selectline:no current line"}` — **zero lines** |
| `:441` from this laptop | timeout | timeout |
| `:443` from this laptop | open | open |
| `remote.py doctor` | 0 ms accept, 0-byte banner | 32 ms accept, 0-byte banner |
| `ssh -vvv` | `kex_exchange_identification: read: Unknown error` | identical |

So the line address and port did **not** change; what changed is that the SDP
controller no longer offers any line for this node group. A client cannot select
a line that the controller does not publish, and the previous line's port is
still closed to the internet (checked again at 01:35; four external nodes had
already failed at 00:35 and the control on `:443` succeeded).

The break began at 23:59-00:00 exactly, which is the shape of a scheduled
maintenance window rather than a configuration drift. There is no client-side
action left: the next check should be after the school's service window, and the
question for their administrator is whether the aTrust SDP gateway is down or
its line endpoint changed.
