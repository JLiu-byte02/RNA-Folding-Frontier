#!/usr/bin/env bash
#SBATCH --job-name=rff_baselines
#SBATCH --partition=g-3090-2
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/home/junxuanliu/projects/RNA-Folding-Frontier/logs/baselines_%j.out
# CPU-only: no --gres, so this makes no claim on a GPU.
#
# Compute nodes have NO outbound network (curl fails with "Could not resolve
# host"), so every download happens on the login node and this job only builds
# or unpacks from local tarballs in ${TOOLS}.  Pinned digests (login node):
#   EternaFold.tar.gz    5cf3a40feae73a8e79acd9b21f5c3324ae22a4a511d32ad276e3039c73a90e84
#   RNAstructure.tar.gz  4b419a3ec13e107311ea99050da374f338aedaf48b3508d75a5b0455015ecb13
# Layout facts (verified 2026-11-13):
#   * EternaFold has NO top-level Makefile; the build lives in src/ and the
#     generated Defaults.ipp ships in the tarball.
#   * RNAstructure.tar.gz already contains compiled exe/ binaries under a
#     top-level RNAstructure/ directory, so it is unpacked, not built.
set -euo pipefail

ROOT="/home/junxuanliu/projects/RNA-Folding-Frontier"
TOOLS="/home/junxuanliu/tools"
mkdir -p "${TOOLS}"
cd "${TOOLS}"

check_tar() {
  local f="$1"
  if [ ! -s "$f" ]; then echo "[FATAL] missing $f"; exit 1; fi
  if ! tar tzf "$f" >/dev/null 2>&1; then echo "[FATAL] corrupt $f"; exit 1; fi
}

echo "[start] $(date -Iseconds)"
sha256sum EternaFold.tar.gz RNAstructure.tar.gz || true

if [ ! -x "${TOOLS}/EternaFold/src/contrafold" ]; then
  check_tar EternaFold.tar.gz
  rm -rf EternaFold EternaFold-master
  tar xzf EternaFold.tar.gz
  mv EternaFold-master EternaFold
  echo "[eternafold build] $(date -Iseconds)"
  make -C "${TOOLS}/EternaFold/src" -j "${SLURM_CPUS_PER_TASK}" contrafold
fi
echo "[eternafold done] $(date -Iseconds)"
ls -la "${TOOLS}/EternaFold/src/contrafold"

if [ ! -x "${TOOLS}/RNAstructure/exe/Fold" ]; then
  check_tar RNAstructure.tar.gz
  rm -rf RNAstructure
  tar xzf RNAstructure.tar.gz -C "${TOOLS}"
fi
echo "[rnastructure done] $(date -Iseconds)"
ls -la "${TOOLS}/RNAstructure/exe/Fold"
"${TOOLS}/RNAstructure/exe/Fold" --version 2>&1 | head -3 || true
echo "[all done] $(date -Iseconds)"