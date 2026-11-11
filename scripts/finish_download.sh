#!/usr/bin/env bash
set -euo pipefail
cd /home/junxuanliu/projects/RNA-Folding-Frontier

URL="https://github.com/heqin-zhu/BPfold/releases/download/v0.2/BPfold_data.tar.gz"
OUT="data/raw/BPfold_data.tar.gz"
SIZE=166902534
CHUNKS=24
CHUNK=$(( (SIZE + CHUNKS - 1) / CHUNKS ))

DIR=$(ls -d "${OUT}".parts.* | head -1)
echo "using $DIR"

for i in $(seq 0 $((CHUNKS - 1))); do
  start=$(( i * CHUNK ))
  end=$(( (i + 1) * CHUNK - 1 ))
  [ "$end" -ge "$SIZE" ] && end=$(( SIZE - 1 ))
  want=$(( end - start + 1 ))
  f="$DIR/part.$(printf '%03d' "$i")"
  have=0
  [ -f "$f" ] && have=$(stat -c%s "$f")
  if [ "$have" -ne "$want" ]; then
    echo "refetch $i: $have/$want"
    curl -sL --retry 5 -r "${start}-${end}" -o "$f" "$URL"
  fi
done

cat "$DIR"/part.* > "$OUT"
actual=$(stat -c%s "$OUT")
echo "size: $actual (expected $SIZE)"
[ "$actual" -eq "$SIZE" ]
sha256sum "$OUT"
rm -rf "$DIR"
