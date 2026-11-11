#!/usr/bin/env bash
set -u
echo "== archiveII =="
curl -sIL -m 25 https://rna.urmc.rochester.edu/pub/archiveII.tar.gz \
  | tr -d '\r' | grep -iE '^HTTP/|^content-length|^location' | tail -4

echo "== bpRNA download page links =="
curl -sL -m 25 https://bprna.cgrb.oregonstate.edu/download.php \
  | grep -ioE 'href="[^"]*(zip|tar|gz|tgz)[^"]*"' | head -30

echo "== Rfam dir (raw grep) =="
curl -sL -m 25 https://ftp.ebi.ac.uk/pub/databases/Rfam/15.0/ | head -c 1200
echo
