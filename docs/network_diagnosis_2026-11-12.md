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
