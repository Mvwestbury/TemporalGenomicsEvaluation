#!/usr/bin/env bash
#SBATCH --job-name=rb_reads
#SBATCH --partition=all
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --output=/path/to/TemporalGenomicsEvaluation/logs/rb_reads_%j.out
# Reads for the 10 robins, in the Liu layout raw_reads/Petroica_traversi/<run>/ with a .md5ok flag per file
# (MD5 from the ENA file report, robin/robin_files.tsv):
#   modern (paired)   : downloaded from ENA over HTTPS (resumed, MD5-checked), 2 files at a time -> <run>_r1/_r2.fq.gz
#   historical (SE)   : already downloaded in GT_damage/robin/robin_historical/raw_data -> symlinked as <run>.fq.gz,
#                       MD5-checked against ENA (re-downloaded from ENA if the local copy is missing or differs)
set -uo pipefail
source /path/to/TemporalGenomicsEvaluation/scripts/safe_mkdir.sh   # race-safe mkdir -p (uutils mkdir fails on NFS races)
B=/path/to/TemporalGenomicsEvaluation
H=/path/to/GT_damage/robin/robin_historical/raw_data
D=$B/raw_reads/Petroica_traversi
echo "host $(hostname) start $(date)"

fetch() {  # run read url md5
    local run=$1 rd=$2 url=$3 md5=$4 dir=$D/$1
    local f=$dir/$run.fq.gz; [[ $rd == se ]] || f=$dir/${run}_$rd.fq.gz
    mkdir -p $dir
    [[ -f $f.md5ok ]] && { echo "SKIP $f"; return 0; }
    if [[ $rd == se && -s $H/$run.fastq.gz ]]; then
        ln -sfn $H/$run.fastq.gz $f
        [[ "$(md5sum < $f | cut -d' ' -f1)" == "$md5" ]] && { touch $f.md5ok; echo "OK   $f (local copy)"; return 0; }
        echo "BADMD5 local copy $f, downloading"; rm -f $f
    fi
    for attempt in 1 2 3; do
        local drops=0
        until wget -q -c --tries=20 --retry-connrefused --waitretry=30 --read-timeout=120 -O $f "$url"; do
            drops=$((drops + 1)); [[ $drops -ge 50 ]] && { echo "FAIL (no progress) $f"; return 1; }; sleep 30
        done
        [[ "$(md5sum < $f | cut -d' ' -f1)" == "$md5" ]] && { touch $f.md5ok; echo "OK   $f (resumed $drops times)"; return 0; }
        echo "BADMD5 attempt $attempt $f"; rm -f $f
    done
    echo "FAIL $f"; return 1
}
export -f fetch; export D H

# historical: local MD5 checks, 4 at a time; modern: downloads, 2 at a time (run concurrently)
awk -F'\t' 'NR>1 && $4=="se"{print $3"\t"$4"\t"$5"\t"$6}' $B/robin/robin_files.tsv |
    xargs -P 4 -d '\n' -I{} bash -c 'IFS=$'"'"'\t'"'"' read -r a b c d <<< "{}"; fetch $a $b $c $d' &
awk -F'\t' 'NR>1 && $4!="se"{print $3"\t"$4"\t"$5"\t"$6}' $B/robin/robin_files.tsv |
    xargs -P 2 -d '\n' -I{} bash -c 'IFS=$'"'"'\t'"'"' read -r a b c d <<< "{}"; fetch $a $b $c $d' &
wait
n=$(find $D -name '*.md5ok' | wc -l); echo "verified files: $n of $(( $(wc -l < $B/robin/robin_files.tsv) - 1 ))"
echo "end $(date)"
[[ $n -eq $(( $(wc -l < $B/robin/robin_files.tsv) - 1 )) ]]
