#!/usr/bin/env bash
# Parallel ranged downloader for throttled single-connection links.
set -euo pipefail

URL="$1"
OUT="$2"
SIZE="$3"
CHUNKS="${4:-16}"

DIR="${OUT}.parts"
mkdir -p "$DIR"

CHUNK=$(( (SIZE + CHUNKS - 1) / CHUNKS ))

fetch_part() {
  local i=$1 start=$2 end=$3 want=$4 f=$5
  local have=0
  [ -f "$f" ] && have=$(stat -c%s "$f")
  [ "$have" -eq "$want" ] && return 0
  curl -sL --max-time 600 --retry 3 -r "${start}-${end}" -o "$f" "$URL" || true
  have=0
  [ -f "$f" ] && have=$(stat -c%s "$f")
  [ "$have" -eq "$want" ]
}

for attempt in 1 2 3; do
  pids=()
  for i in $(seq 0 $((CHUNKS - 1))); do
    start=$(( i * CHUNK ))
    end=$(( (i + 1) * CHUNK - 1 ))
    [ "$end" -ge "$SIZE" ] && end=$(( SIZE - 1 ))
    [ "$start" -gt "$end" ] && continue
    want=$(( end - start + 1 ))
    f="$DIR/part.$(printf '%03d' "$i")"
    fetch_part "$i" "$start" "$end" "$want" "$f" &
    pids+=($!)
  done
  for p in "${pids[@]}"; do wait "$p" || true; done
  # Verify; break when all parts are complete.
  ok=1
  for i in $(seq 0 $((CHUNKS - 1))); do
    start=$(( i * CHUNK ))
    end=$(( (i + 1) * CHUNK - 1 ))
    [ "$end" -ge "$SIZE" ] && end=$(( SIZE - 1 ))
    want=$(( end - start + 1 ))
    f="$DIR/part.$(printf '%03d' "$i")"
    have=0
    [ -f "$f" ] && have=$(stat -c%s "$f")
    [ "$have" -eq "$want" ] || { ok=0; echo "incomplete $i: $have/$want"; }
  done
  [ "$ok" -eq 1 ] && break
  echo "attempt $attempt incomplete; retrying"
done

[ "$ok" -eq 1 ] || { echo "giving up; partial parts kept in $DIR" >&2; exit 1; }

cat "$DIR"/part.* > "$OUT"
actual=$(stat -c%s "$OUT")
echo "size: $actual (expected $SIZE)"
[ "$actual" -eq "$SIZE" ] || { echo "size mismatch" >&2; exit 1; }
sha256sum "$OUT"
rm -rf "$DIR"
