#!/usr/bin/env bash
# Download Rfam 15.0 artefacts (seed alignments, covariance models, clan map),
# verify md5 where the release ships a checksum file, and press the CMs.
set -uo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DEST="${ROOT}/data/raw/rfam"
BASE="https://ftp.ebi.ac.uk/pub/databases/Rfam/15.0"
mkdir -p "${DEST}"
cd "${DEST}"

fetch() {
  local name="$1"
  if [ -s "${name}" ]; then
    echo "have ${name}"
    return 0
  fi
  echo "fetching ${name}"
  wget -c --tries=50 --waitretry=30 --timeout=120 -O "${name}.part" "${BASE}/${name}" \
    && mv "${name}.part" "${name}"
}

fetch Rfam.seed.gz
fetch Rfam.cm.gz
fetch Rfam.clanin
fetch md5.txt

if [ -s md5.txt ]; then
  grep -E 'Rfam\.(seed|cm)\.gz|Rfam\.clanin' md5.txt > md5.subset || true
  echo "--- md5 check ---"
  md5sum -c md5.subset || echo "WARNING: checksum mismatch"
fi

echo "--- press CMs ---"
if [ ! -s Rfam.cm ]; then
  gzip -dc Rfam.cm.gz > Rfam.cm
fi
CMPRESS="${ROOT}/tools/infernal/bin/cmpress"
if [ -x "${CMPRESS}" ]; then
  "${CMPRESS}" -F Rfam.cm && echo "pressed"
else
  echo "cmpress not found; run scripts/server_build_infernal.sh first"
fi

ls -la "${DEST}"
