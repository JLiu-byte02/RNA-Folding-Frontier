#!/usr/bin/env bash
# Probe public data sources for accessibility and size.
set -u

urls=(
  "https://rna.urmc.rochester.edu/pub/archiveII/"
  "https://rna.urmc.rochester.edu/pub/archiveII/archiveII.tar.gz"
  "https://rna.urmc.rochester.edu/pub/archiveII/ArchiveII_sequences2.ct"
  "https://ftp.ebi.ac.uk/pub/databases/Rfam/15.0/family.txt.gz"
  "https://ftp.ebi.ac.uk/pub/databases/Rfam/15.0/Rfam.cm.gz"
  "https://ftp.ebi.ac.uk/pub/databases/Rfam/15.0/Rfam.clanin"
  "https://bprna.cgrb.oregonstate.edu/"
)

for u in "${urls[@]}"; do
  code=$(curl -sIL -m 20 -o /dev/null -w '%{http_code}' "$u" || true)
  len=$(curl -sIL -m 20 "$u" 2>/dev/null | tr -d '\r' | awk 'tolower($1)=="content-length:"{v=$2} END{print v}')
  echo "$code  ${len:-?}  $u"
done
