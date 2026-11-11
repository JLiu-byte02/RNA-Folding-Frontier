#!/usr/bin/env bash
# Build Infernal (cmsearch/cmscan/cmpress) from source into <root>/tools/infernal.
# Idempotent: does nothing if the binaries already exist.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
mkdir -p "${ROOT}/tools"
cd "${ROOT}/tools"

if [ -x "${ROOT}/tools/infernal/bin/cmscan" ]; then
  echo "infernal already installed"
  exit 0
fi

VER=1.1.5
if [ ! -f "infernal-${VER}.tar.gz" ]; then
  wget -c --tries=20 --waitretry=10 --timeout=120 \
    -O "infernal-${VER}.tar.gz" \
    "http://eddylab.org/software/infernal/infernal-${VER}.tar.gz" ||
  wget -c --tries=20 --waitretry=10 --timeout=120 \
    -O "infernal-${VER}.tar.gz" \
    "https://github.com/EddyRivasLab/infernal/archive/refs/tags/infernal-${VER}.tar.gz"
fi

rm -rf "infernal-${VER}"
tar xzf "infernal-${VER}.tar.gz"
cd "infernal-${VER}"
./configure --prefix="${ROOT}/tools/infernal"
make -j 16
make install

echo "installed:"
"${ROOT}/tools/infernal/bin/cmscan" -h | head -3
