#!/usr/bin/env bash
#SBATCH --job-name=full_map
#SBATCH --partition=all
#SBATCH --cpus-per-task=10
#SBATCH --mem=32G
#SBATCH --output=/path/to/TemporalGenomicsEvaluation/logs/full_map_%A_%a.out
# Full-data mapping with PlainMap in modern mode (bwa mem). One array task per data line of
# metadata/full_map_samples.tsv (species, voucher, CRR, period, ref):
#   sbatch --array=1-<nrows>%10 full_map.sbatch
# Runs whose reads are not downloaded and MD5-verified yet are skipped; resubmit their task IDs later.
set -euo pipefail
source /path/to/TemporalGenomicsEvaluation/scripts/safe_mkdir.sh   # race-safe mkdir -p (uutils mkdir fails on NFS races)
B=/path/to/TemporalGenomicsEvaluation
PLAINMAP=/path/to/software/PlainMap/plainmap.sh
T=${SLURM_CPUS_PER_TASK:-10}

IFS=$'\t' read -r species voucher crr period ref <<< "$(awk -F'\t' -v n=$((SLURM_ARRAY_TASK_ID + 1)) 'NR==n' $B/metadata/full_map_samples.tsv)"
sp=${species// /_}
reads=$B/raw_reads/$sp/$crr
prefix=${voucher}_${crr}
outdir=$B/mapped/modern/$sp
echo "task $SLURM_ARRAY_TASK_ID: $species $voucher $crr $period ref=$ref host=$(hostname) $(date)"

for f in $reads/${crr}_r1.fq.gz $reads/${crr}_r2.fq.gz; do
    [[ -f $f.md5ok ]] || { echo "SKIP: $f not downloaded/verified yet"; exit 0; }
done
[[ -f $ref.bwt ]] || { echo "SKIP: reference $ref not indexed yet"; exit 0; }

mkdir -p $B/mapped/manifests $outdir
manifest=$B/mapped/manifests/$prefix.txt
printf '%s\n' $reads/${crr}_r1.fq.gz $reads/${crr}_r2.fq.gz > $manifest

bash $PLAINMAP \
    --manifest $manifest \
    --prefix $prefix \
    --ref $ref \
    --outdir $outdir \
    --library-type modern \
    --threads $T
echo "done $(date)"
