#!/usr/bin/env bash
set -u
echo "== Rfam 15.0 listing =="
curl -sL -m 25 https://ftp.ebi.ac.uk/pub/databases/Rfam/15.0/ \
  | grep -oE 'href="[^"]+"' | sed 's/href="//; s/"$//' | head -40

echo "== URMC root archive mentions =="
curl -sL -m 25 https://rna.urmc.rochester.edu/ \
  | grep -ioE 'href="[^"]*archive[^"]*"' | head -20

echo "== bpRNA downloads =="
curl -sL -m 25 https://bprna.cgrb.oregonstate.edu/ \
  | grep -ioE 'href="[^"]*"' | head -40
