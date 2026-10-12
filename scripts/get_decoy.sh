#!/usr/bin/env bash
# Build a NON-viral protein decoy for the reciprocal (reverse) search.
# Usage: bash scripts/get_decoy.sh swissprot [outdir]   ~570k reviewed proteins, ~90 MB, minutes
#        bash scripts/get_decoy.sh uniref50  [outdir]   ~60M+ clusters, ~12 GB download, hours
# Output: <outdir>/decoy_<mode>_nonviral.faa.gz -> config `cellular_decoy_faa`
#
# Viruses are removed by taxonomy (NCBI taxon 10239 and all descendants), so EVEs are only
# rejected when their best match is a cellular protein.
set -euo pipefail
MODE=${1:-swissprot}
OUT=${2:-db/decoy}
mkdir -p "$OUT" && cd "$OUT"
REST=https://rest.uniprot.org

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

case "$MODE" in
  swissprot)
    # Reviewed UniProtKB entries, excluding viruses; streamed and gzip-compressed by UniProt
    curl -sSf --retry 3 -o decoy_swissprot_nonviral.faa.gz \
      "$REST/uniprotkb/stream?compressed=true&format=fasta&query=%28reviewed%3Atrue%29%20NOT%20%28taxonomy_id%3A10239%29"
    ;;
  uniref50)
    # 1) IDs of UniRef50 clusters whose common taxon is viral (small list)
    curl -sSf --retry 3 -o viral_uniref50.ids.gz \
      "$REST/uniref/stream?compressed=true&format=list&query=%28identity%3A0.5%29%20AND%20%28taxonomy_id%3A10239%29"
    # 2) full UniRef50 from the FTP mirror (much faster than streaming 60M records over REST)
    fetch uniref50.fasta.gz https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz
    # 3) drop the viral clusters
    zcat viral_uniref50.ids.gz > viral_uniref50.ids
    seqkit grep -v -f viral_uniref50.ids uniref50.fasta.gz -o decoy_uniref50_nonviral.faa.gz
    rm -f viral_uniref50.ids
    ;;
  *) echo "mode must be swissprot or uniref50" >&2; exit 1 ;;
esac

F=decoy_${MODE}_nonviral.faa.gz
echo "$(zgrep -c '^>' "$F") non-viral proteins -> $PWD/$F"
