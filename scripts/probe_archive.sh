#!/usr/bin/env bash
set -u
for u in \
  "http://rna.urmc.rochester.edu/pub/archiveII.tar.gz" \
  "https://rna.urmc.rochester.edu/pub/archiveII.tar.gz" \
  "https://rna.urmc.rochester.edu/pub/" \
  "https://rna.urmc.rochester.edu/RNAstructureWeb/Download/ArchiveII.html"; do
  echo "==== $u"
  curl -sIL -m 20 -A 'Mozilla/5.0' "$u" | tr -d '\r' \
    | grep -iE '^HTTP/|^content-length|^location' | tail -4
done
