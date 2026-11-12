# Remote provenance — 2026-11-12

This records the A6.1 recovery after SSH returned, and the one provenance defect
the recovery exposed.

## Connectivity

At 14:44 +0800 a non-sandboxed `scripts/remote.py doctor` reported TCP connect
in 12 ms and `SSH-2.0-OpenSSH_7.4`. A real login succeeded on `master` as
`junxuanliu`, with `squeue` empty. The public `:441` probe still reports closed,
so it is not a reliable liveness test for the aTrust redirect; doctor plus an
authenticated command is.

## The project directory was not a git repository

`/home/junxuanliu/projects/RNA-Folding-Frontier` had 55 GB of data and all four
split files, but no `.git`. Every earlier “server copy” was therefore a file
copy, not a traceable checkout. The recovery procedure was:

1. Back up every non-data directory to
   `/home/junxuanliu/projects/_rff_provenance_backup/code_20261112_1448.tgz`
   (395 MB). `data/`, `runs/`, `checkpoints/`, and `logs/` were excluded from
   the archive and left untouched.
2. `git init`, add the public GitHub remote, fetch `main`, and check out
   `origin/main` with `-f -B main`.
3. Record the resulting HEAD in the project state and verify `git status`.

After the checkout the only untracked paths are data/results that are
intentionally server-side:

```
experiments/000_e0_leakage/...
experiments/001_e1prime/crop_sanity_giga.json
experiments/002_e1prime_runs/diag_*_perseq.jsonl
splits/
```

`data/`, `runs/`, and `checkpoints/` are ignored. No server data was deleted.

## Frozen split hashes

Computed on the server after the checkout:

| file | sha256 |
|---|---|
| `splits/s1prime/train.jsonl` | `0e3711f65c6eec001c82a3c5110561afaa8cd5c3033d30e5ad4e19a2b64bade2` |
| `splits/s1prime/id_test.jsonl` | `f40324ec033543a99ff4cc10dc4f5637c20b6ebb29c2c7d6c53c51a8db88e137` |
| `splits/s1prime/test_withheld.jsonl` | `bca8ae19baee5f6c9141b973b88e3f6092a83603bdca912837cb522e8acb9b9b` |
| `splits/s1prime/test_novel.jsonl` | `ab8eaeee608055127c7f8c8fbc27c95064dbdfb870e752f8006d319cee68d09d` |
| `splits/s1prime/manifest.json` | `a9ceb336aef1f0a6b10127138ff35751350ebbe86518d0589aaed0e60b7440fe` |
| `splits/s1prime_filtered/train.jsonl` | `b8af23357fc6c36aab0d4ebae1647bde452582ea8a0872590efd747fa1f94ecd` |
| `splits/s1prime_filtered/id_test.jsonl` | `4c280160357559070e7f44efd40ef3eff5be373fe8f83e36ca4f7c4b49d7cccf` |

The filtered `test_novel` and `test_withheld` files are byte-identical to their
`s1prime` counterparts. Membership was not rebuilt or changed.

## Two artifact defects found and fixed

1. **Git treated gzip as text.** `.gitattributes` had a global `text=auto
   eol=lf`. A gzip stream may contain the bytes `0d0a`; git’s clean filter then
   strips `0d` and reports a phantom modification, and a later `git add` would
   store the altered bytes. Gzip, `.pt`, `.h5`, `.npz`, `.pkl`, `.safetensors`
   and checkpoint suffixes are now explicitly `binary`.
2. **Gzip OS byte was platform-dependent.** The compressed stream was identical
   on Windows and Linux, but the header OS byte was `0x0a` here and `0xff` on
   the server, so the artifact sha256 differed. `gzip_deterministic` pins that
   metadata byte to `0x0a`; the server now reproduces both frozen hashes:

```
labels_v2_delta.json.gz   27bbaa32d11837c77b95207c8eac00dbc81890d55eb72e4f1f509f0626deb0af
labels_v3_pairs.json.gz   b8fe71abfab30c97298bc90140bd82fc7b787f2cc439b838d86a1b0b6b895ec1
```

No frozen label content changed. The v3 summary now carries server-side
per-split counts; see `experiments/011_labels_v3/README.md`.
