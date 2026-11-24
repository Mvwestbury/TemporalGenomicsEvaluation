#!/usr/bin/env bash
# Segmented downloader for GSA: the server throttles each connection and drops long ones,
# so each file is fetched as many short HTTP range requests in parallel, then joined.
#
# Usage: fetch_chunked.sh FILE_LIST MODE [N_PARALLEL] [CRR_FILTER]
#   FILE_LIST  TSV: species, CRR, url, md5, bytes (metadata/four_species_files.tsv)
#   MODE       full  -> raw_reads/<species>/<CRR>/<file>, MD5-verified (.md5ok marker)
#              head  -> pilot_heads/<species>/<CRR>/<file>: first HEAD_READS reads of each
#                       mate, cut to the same read count in R1 and R2 and re-gzipped
#   CRR_FILTER optional file with one CRR per line; only those runs are fetched
set -uo pipefail
source /path/to/TemporalGenomicsEvaluation/scripts/safe_mkdir.sh   # race-safe mkdir -p (uutils mkdir fails on NFS races)
B=/path/to/TemporalGenomicsEvaluation
LIST=$1
MODE=$2
NPAR=${3:-12}
FILTER=${4:-}
CHUNK=$((64 * 1000 * 1000))
[[ $MODE == head ]] && CHUNK=$((10 * 1000 * 1000))
HEAD_BYTES=$((200 * 1000 * 1000))   # comfortably more than HEAD_READS reads at ~100 bytes/read
HEAD_READS=1100000

# fetch one byte range until it has exactly the expected length
fetch_chunk() {
    local url=$1 a=$2 b=$3 out=$4
    local want=$((b - a + 1)) try
    [[ -f $out && $(stat -c %s "$out") -eq $want ]] && return 0
    for try in $(seq 1 60); do
        curl -s -r "$a-$b" --connect-timeout 30 --max-time 900 -o "$out.tmp" "$url" &&
            [[ $(stat -c %s "$out.tmp") -eq $want ]] && { mv "$out.tmp" "$out"; return 0; }
        sleep $((10 + RANDOM % 30))
    done
    echo "CHUNKFAIL $out"; return 1
}
export -f fetch_chunk

fetch_file() {
    local sp=$1 crr=$2 url=$3 md5=$4 size=$5
    local name dir dest parts end n i
    name=$(basename "$url")
    if [[ $MODE == full ]]; then dir=$B/raw_reads/$sp/$crr; end=$size
    else dir=$B/pilot_heads/$sp/$crr; end=$(( HEAD_BYTES < size ? HEAD_BYTES : size )); fi
    dest=$dir/$name; parts=$dest.parts
    [[ -f $dest.md5ok || -f $dest.headok ]] && { echo "SKIP $dest"; return 0; }
    mkdir -p "$parts"
    n=$(( (end + CHUNK - 1) / CHUNK ))
    for ((i = 0; i < n; i++)); do
        printf '%s %d %d %s/%06d\n' "$url" $((i * CHUNK)) $(( (i + 1) * CHUNK < end ? (i + 1) * CHUNK - 1 : end - 1 )) "$parts" $i
    done | xargs -P "$NPAR" -L 1 bash -c 'fetch_chunk "$@"' _ || { echo "FAIL $dest"; return 1; }
    cat "$parts"/[0-9]* > "$dest.tmp" && rm -rf "$parts"
    if [[ $MODE == full ]]; then
        if [[ $(md5sum "$dest.tmp" | cut -d' ' -f1) == "$md5" ]]; then
            mv "$dest.tmp" "$dest"; touch "$dest.md5ok"; echo "OK   $dest"
        else
            rm -f "$dest.tmp"; echo "BADMD5 $dest"; return 1
        fi
    else
        mv "$dest.tmp" "$dest.trunc"; echo "HEAD $dest.trunc"
    fi
}

# head mode: turn the truncated gz pair into valid, equal-length FASTQs
finish_heads() {
    local dir=$1 crr r1 r2 n1 n2 n
    crr=$(basename "$dir")
    r1=$dir/${crr}_r1.fq.gz; r2=$dir/${crr}_r2.fq.gz
    [[ -f $r1.headok && -f $r2.headok ]] && return 0
    [[ -f $r1.trunc && -f $r2.trunc ]] || return 0
    n1=$( (gzip -dc "$r1.trunc" 2>/dev/null || true) | head -n $((HEAD_READS * 4)) | awk 'END{print int(NR/4)}')
    n2=$( (gzip -dc "$r2.trunc" 2>/dev/null || true) | head -n $((HEAD_READS * 4)) | awk 'END{print int(NR/4)}')
    n=$(( n1 < n2 ? n1 : n2 ))
    for m in r1 r2; do
        (gzip -dc "$dir/${crr}_$m.fq.gz.trunc" 2>/dev/null || true) | head -n $((n * 4)) | gzip -1 > "$dir/${crr}_$m.fq.gz"
    done
    # mates must pair up record by record
    if paste <(gzip -dc "$r1" | awk 'NR%4==1{sub(/\/1$/,"",$1); print $1}') \
             <(gzip -dc "$r2" | awk 'NR%4==1{sub(/\/2$/,"",$1); print $1}') | awk '$1!=$2{bad++} END{exit bad>0}'; then
        touch "$r1.headok" "$r2.headok"; rm -f "$r1.trunc" "$r2.trunc"; echo "HEADOK $dir $n reads"
    else
        echo "HEADPAIRFAIL $dir"
    fi
}

export B MODE NPAR CHUNK HEAD_BYTES
tail -n +2 "$LIST" | while IFS=$'\t' read -r sp crr url md5 size; do
    [[ -n $FILTER ]] && ! grep -qx "$crr" "$FILTER" && continue
    fetch_file "$sp" "$crr" "$url" "$md5" "$size"
done
if [[ $MODE == head ]]; then
    for d in "$B"/pilot_heads/*/*/; do finish_heads "${d%/}"; done
fi
echo "finished $MODE $(date)"
