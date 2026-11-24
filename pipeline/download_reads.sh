#!/usr/bin/env bash
# Download GSA reads listed in a run table (species, voucher, period, CRX, CRR, gsa_name, r1_url, r1_md5, r2_url, r2_md5, bytes)
# Uses HTTPS mirror (FTP is often overloaded), resumes partial files, verifies MD5.
# Usage: download_reads.sh metadata/four_species_runs.tsv [n_parallel]
set -uo pipefail
source /path/to/TemporalGenomicsEvaluation/scripts/safe_mkdir.sh   # race-safe mkdir -p (uutils mkdir fails on NFS races)
B=/path/to/TemporalGenomicsEvaluation
TABLE=$1
NPAR=${2:-3}
LOG=$B/logs/download_$(date +%Y%m%d_%H%M).log

fetch() {
    local sp=$1 crr=$2 url=$3 md5=$4
    local dir=$B/raw_reads/${sp// /_}/$crr
    local f=$dir/$(basename "$url")
    url=${url/ftp:\/\/download.big.ac.cn/https:\/\/download.cncb.ac.cn}
    mkdir -p "$dir"
    if [[ -f $f.md5ok ]]; then echo "SKIP $f"; return 0; fi
    for attempt in 1 2 3; do
        # the server drops connections often: keep resuming the partial file until wget
        # reports the transfer complete (gives up after 200 consecutive failed sessions)
        local drops=0
        until wget -q -c --tries=20 --retry-connrefused --waitretry=30 --read-timeout=120 -O "$f" "$url"; do
            drops=$((drops + 1))
            [[ $drops -ge 200 ]] && { echo "FAIL (no progress) $f"; return 1; }
            sleep 30
        done
        if [[ "$(md5sum "$f" | cut -d' ' -f1)" == "${md5,,}" ]]; then
            touch "$f.md5ok"; echo "OK   $f (resumed $drops times)"; return 0
        fi
        # complete transfer but wrong checksum: the file is corrupt, start over
        echo "BADMD5 attempt $attempt $f"; rm -f "$f"
    done
    echo "FAIL $f"; return 1
}
export -f fetch; export B

tail -n +2 "$TABLE" | awk -F'\t' '{print $1"\t"$5"\t"$7"\t"$8; print $1"\t"$5"\t"$9"\t"$10}' |
    xargs -P "$NPAR" -d '\n' -I{} bash -c 'IFS=$'"'"'\t'"'"' read -r sp crr url md5 <<< "{}"; fetch "$sp" "$crr" "$url" "$md5"' \
    >> "$LOG" 2>&1
echo "done; $(grep -c '^OK' "$LOG") ok, $(grep -c '^FAIL' "$LOG") failed" >> "$LOG"
