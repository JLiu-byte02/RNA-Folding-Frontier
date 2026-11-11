#!/usr/bin/env bash
set -u
cd /home/junxuanliu/projects/RNA-Folding-Frontier/data/raw/bpfold/BPfold_data

for d in bpRNA/TR0 bpRNA/VL0 bpRNA/TS0 archiveII/archiveII bpRNAnew/bpRNAnew.nr500.canonicals RNAStrAlign_bpseq; do
  echo "==== $d  (files: $(find "$d" -type f | wc -l))"
  find "$d" -type f | head -3
done

echo "==== sample .ct"
sample=$(find bpRNA/TR0 -type f | head -1)
echo "$sample"
head -8 "$sample"

echo "==== sample .bpseq"
sample=$(find archiveII/archiveII -type f | head -1)
echo "$sample"
head -8 "$sample"
