#!/usr/bin/env bash
# Download GSA read files with wget -c: one connection per file, resumed after every drop, MD5-verified.
# Usage: wget_files.sh FILE_LIST [N_PARALLEL]
#   FILE_LIST  TSV: species, CRR, url, md5, bytes (metadata/four_species_files.tsv or metadata/split/files_NN.tsv)
# Chunks left by fetch_chunked.sh (<file>.parts/) are joined into the start of the file so they are not refetched.
set -uo pipefail
source /path/to/TemporalGenomicsEvaluation/scripts/safe_mkdir.sh   # race-safe mkdir -p (uutils mkdir fails on NFS races)
B=/path/to/TemporalGenomicsEvaluation
LIST=$1
NPAR=${2:-3}

fetch() {
    local sp=$1 crr=$2 url=$3 md5=$4
    local dir=$B/raw_reads/$sp/$crr f i drops attempt
    f=$dir/$(basename "$url")
    [[ -f $f.md5ok ]] && { echo "SKIP $f"; return 0; }
    mkdir -p "$dir"
    # finished chunks all have the full chunk size, so the leading run of them is a valid partial file
    if [[ -d $f.parts && ! -f $f ]]; then
        i=0
        while [[ -f $(printf '%s.parts/%06d' "$f" $i) ]]; do cat "$(printf '%s.parts/%06d' "$f" $i)" >> "$f"; i=$((i + 1)); done
        rm -rf "$f.parts"; echo "REUSED $i chunks $f"
    fi
    for attempt in 1 2 3; do
        # the server drops and refuses connections in bursts: keep resuming until wget reports
        # the transfer complete (gives up after 500 consecutive failed sessions)
        drops=0
        until wget -q -c --tries=20 --retry-connrefused --waitretry=30 --read-timeout=120 -O "$f" "$url"; do
            drops=$((drops + 1))
            [[ $drops -ge 500 ]] && { echo "FAIL (no progress) $f"; return 1; }
            sleep 30
        done
        if [[ $(md5sum "$f" | cut -d' ' -f1) == "$md5" ]]; then
            touch "$f.md5ok"; echo "OK   $f (resumed $drops times) $(date)"; return 0
        fi
        # complete transfer but wrong checksum: the file is corrupt, start over
        echo "BADMD5 attempt $attempt $f"; rm -f "$f"
    done
    echo "FAIL $f"; return 1
}
export -f fetch; export B

tail -n +2 "$LIST" | cut -f1-4 |
    xargs -P "$NPAR" -d '\n' -I{} bash -c 'IFS=$'"'"'\t'"'"' read -r sp crr url md5 <<< "{}"; fetch "$sp" "$crr" "$url" "$md5"'
echo "finished $(date)"
