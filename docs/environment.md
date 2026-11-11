# Compute Environment

## Server

- Host: `10.22.19.39` (login node `master`), user `junxuanliu`.
- Project root: `/home/junxuanliu/projects/RNA-Folding-Frontier`.
- Scheduler: Slurm. GPU partition `g-3090-2`, 4 nodes x 2 RTX 3090 (24 GB).
- Local project cap: **6 GPUs**.

## Python

- Conda env: `/home/junxuanliu/.conda/envs/rff` (Python 3.10.18).
- PyTorch 2.5.1 + CUDA 12.1 (installed).
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

The local machine cannot reach `github.com:443`, but SSH on port 22 works via a
repository-scoped deploy key stored at `~/.ssh/rff_deploy2`. Git is configured
with `core.sshCommand` accordingly.

