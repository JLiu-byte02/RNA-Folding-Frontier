# Public-repository hygiene

Audit: 2026-11-14. The repository is public, so anything committed is permanent.

## What was checked

- every tracked text file, for private-key blocks and for credentials written as
  literal assignments (an SSH password, a login password, an API key, a token);
- every tracked filename, for key-like names (`id_rsa`, `*.pem`, `.env`,
  `secrets.*`);
- the full git history, for the cluster password literal and for OpenSSH/RSA
  private-key blocks (`git log --all -S`).

## Result

**No credentials are present, in the working tree or anywhere in history.** The
two checks are now CI tests (`tests/test_public_repo_hygiene.py`), so the
fresh-clone job fails if one is introduced later.

## Known non-secret identifiers

The repository documents the compute environment, so it contains infrastructure
identifiers: the cluster host and account name (`docs/environment.md`,
`docs/corrections.md`), the campus VPN endpoints and line diagnostics
(`docs/network_diagnosis_2026-11-12.md`, `docs/decision_log.md`), and the SSH
runbook. These are not credentials:

- the cluster address is a private address, unreachable from the public
  internet without the institution's VPN;
- the VPN endpoints are published by the institution's DNS;
- the account name identifies the author, who is also the repository owner.

They are kept because the runbook and the reproduction instructions need them.
They were committed before this audit, so they also remain in history; removing
them from the current tree would not remove them from clones of older commits.

## Policy

1. No credential, token, or private key ever enters the repository, including
   in notebooks, logs, or workflow files.
2. Secrets stay in environment variables on the machine that needs them.
3. If a credential is ever committed: rotate it first, then decide with the
   owner whether a history rewrite is warranted. Rotation is the fix; rewriting
   history is not a substitute for it.
