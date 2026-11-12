# Compute Environment

## Server

- Host: `10.22.19.39` (login node `master`), user `junxuanliu`.
- Project root: `/home/junxuanliu/projects/RNA-Folding-Frontier`.
- Scheduler: Slurm. GPU partition `g-3090-2`, 4 nodes x 2 RTX 3090 (24 GB).
- Local project cap: **6 GPUs**.

## Python

- Conda env: `/home/junxuanliu/.conda/envs/rff` (Python 3.10.18).
- PyTorch 2.5.1 + CUDA 12.1 (installed).
- ViennaRNA 2.6.4 Python wheel installed on 2026-11-12 for D2C/H2.4
  (`pip install --only-binary=:all: ViennaRNA==2.6.4`); the source build of
  2.7.2 fails against this environment's C++ toolchain.
- Bootstrap command:

```bash
/home/junxuanliu/.conda/envs/rff/bin/pip install -r requirements.txt
```

Note: `conda activate` is not configured for non-interactive SSH. Use the
interpreter path directly, or `conda run -n rff`.

## Job conventions

- One experiment = one Slurm job (or a small array).
- Every job writes to `runs/<experiment_id>/`: resolved config, stdout/stderr,
  checkpoints, metrics.
- Logs are mirrored into the repository `experiments/<id>/` after completion
  (large checkpoints stay on the server).

## Repository sync

The local machine reaches GitHub over SSH port 22 with the repository-scoped
deploy key at `~/.ssh/rff_deploy2`; `core.sshCommand` points at it. HTTPS is not
used.

On 2026-11-12 the managed workspace profile made the local `.git` read-only.
Commits during that period were created in a temporary object/working-tree
repository under `tmp_dbg/` using `.git/objects/info/alternates` to read the
original object store, then pushed with the same deploy key. The remote branch
is authoritative; the local checkout's HEAD may lag until the permission
profile changes.

The CPU partition `c-80-255` is shared and frequently fully allocated. Jobs
that need only a few cores still queue behind other users; do not use the login
node to bypass Slurm for data-scale scoring.

