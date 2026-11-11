#!/usr/bin/env bash
set -euo pipefail
cd /home/junxuanliu/projects/RNA-Folding-Frontier
mkdir -p data/raw/rfam

curl -sSL -m 60 -o data/raw/rfam/Rfam.clanin https://ftp.ebi.ac.uk/pub/databases/Rfam/15.0/Rfam.clanin
curl -sSL -m 300 -o data/raw/rfam/Rfam.seed.gz https://ftp.ebi.ac.uk/pub/databases/Rfam/15.0/Rfam.seed.gz
ls -lh data/raw/rfam/

echo "== HF archiveII tree =="
curl -sSL -m 60 "https://huggingface.co/api/datasets/rouskinlab/archiveII/tree/main" | head -c 1500
echo
