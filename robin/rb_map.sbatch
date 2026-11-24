#!/usr/bin/env bash
#SBATCH --job-name=rb_map
#SBATCH --partition=all
#SBATCH --cpus-per-task=10
#SBATCH --mem=32G
#SBATCH --output=/path/to/TemporalGenomicsEvaluation/logs/rb_map_%A_%a.out
# Robin copy of full_map.sbatch: PlainMap in modern mode (bwa mem), row TASK+1 of robin/robin_samples.tsv.
# Modern: one paired run (r1, r2). Historical: the 10 single-end runs of the bird (reads are already
# merged in the ENA files; PlainMap classifies them as SE and maps them with bwa mem as SE).
#   sbatch --array=1-10 rb_map.sbatch
set -euo pipefail
source /path/to/TemporalGenomicsEvaluation/scripts/safe_mkdir.sh   # race-safe mkdir -p (uutils mkdir fails on NFS races)
B=/path/to/TemporalGenomicsEvaluation
PLAINMAP=/path/to/software/PlainMap/plainmap.sh
T=${SLURM_CPUS_PER_TASK:-10}

IFS=$'\t' read -r species voucher crr period ref <<< "$(awk -F'\t' -v n=$((SLURM_ARRAY_TASK_ID + 1)) 'NR==n' $B/robin/robin_samples.tsv)"
sp=${species// /_}
prefix=${voucher}_${crr}
outdir=$B/mapped/modern/$sp
echo "task $SLURM_ARRAY_TASK_ID: $species $voucher $crr $period ref=$ref host=$(hostname) $(date)"
[[ -f $ref.bwt ]] || { echo "ERROR: reference $ref not indexed"; exit 1; }

files=$(awk -F'\t' -v b=$voucher -v D=$B/raw_reads/$sp 'NR>1 && $1==b{
    f = D"/"$3"/"$3 ($4=="se" ? "" : "_"$4) ".fq.gz"; print f}' $B/robin/robin_files.tsv)
for f in $files; do [[ -f $f.md5ok ]] || { echo "ERROR: $f not downloaded/verified"; exit 1; }; done

mkdir -p $B/mapped/manifests $outdir
manifest=$B/mapped/manifests/$prefix.txt
printf '%s\n' $files > $manifest
cat $manifest

bash $PLAINMAP \
    --manifest $manifest \
    --prefix $prefix \
    --ref $ref \
    --outdir $outdir \
    --library-type modern \
    --threads $T
echo "done $(date)"
