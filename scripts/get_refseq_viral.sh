#!/usr/bin/env bash
# Download NCBI RefSeq viral proteins (+ GenPept records carrying the taxonomy lineage).
# Usage: bash scripts/get_refseq_viral.sh [outdir]      (default: db/refseq_viral)
# Output: viral.protein.faa.gz  -> config `viral_faa`
#         viral.protein.gpff.gz -> config `viral_gpff`
set -euo pipefail
OUT=${1:-db/refseq_viral}
BASE=https://ftp.ncbi.nlm.nih.gov/refseq/release/viral
mkdir -p "$OUT" && cd "$OUT"

# resume a partial download, but skip straight past a file that's already complete -- asking
# curl to resume a file the server considers finished gets a 416 and (with set -e) kills the script
fetch() {
  local f=$1 url=$2 remote_size local_size
  remote_size=$(curl -sSIf "$url" | tr -d '\r' | awk -F': ' 'tolower($1)=="content-length"{s=$2} END{print s}')
  local_size=$([ -f "$f" ] && wc -c < "$f" || echo 0)
  if [ -n "$remote_size" ] && [ "$local_size" = "$remote_size" ]; then
    echo "  $f already complete ($remote_size bytes), skipping" >&2
    return 0
  fi
  curl -sSf --retry 3 -C - -o "$f" "$url"
}

# the release is split into numbered parts (viral.1, viral.2, ...); fetch whatever exists
FILES=$(curl -sSf "$BASE/" | grep -oE 'viral\.[0-9]+\.protein\.(faa|gpff)\.gz' | sort -u)
[ -n "$FILES" ] || { echo "could not list $BASE" >&2; exit 1; }
for f in $FILES; do
  fetch "$f" "$BASE/$f"
done
curl -sSf -o RELEASE_NUMBER "https://ftp.ncbi.nlm.nih.gov/refseq/release/RELEASE_NUMBER" \
  || echo "?" > RELEASE_NUMBER        # lives one level up from release/viral/; non-fatal if missing

# concatenated gzip files are valid gzip
cat viral.*.protein.faa.gz  > viral.protein.faa.gz
cat viral.*.protein.gpff.gz > viral.protein.gpff.gz
rm -f viral.[0-9]*.protein.*.gz

echo "RefSeq viral release $(cat RELEASE_NUMBER 2>/dev/null || echo '?'):" \
     "$(zgrep -c '^>' viral.protein.faa.gz) proteins -> $PWD"
