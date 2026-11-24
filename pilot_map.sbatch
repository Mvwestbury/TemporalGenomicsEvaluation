#!/usr/bin/env bash
#SBATCH --job-name=pilot_map
#SBATCH --partition=all
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --output=/path/to/TemporalGenomicsEvaluation/logs/pilot_map_%A_%a.out
# Pilot mapping: 1M fragments per sample with PlainMap, then mapDamage (via PlainMap) and TvCount.
# Historical samples run in both modes (modern = bwa mem, ancient-all = bwa aln);
# modern controls run in modern mode only. One array task per data line of
# metadata/pilot_samples.tsv (species, voucher, CRR, period, ref, mode):
#   sbatch --array=1-<nrows> pilot_map.sbatch
set -euo pipefail
source /path/to/TemporalGenomicsEvaluation/scripts/safe_mkdir.sh   # race-safe mkdir -p (uutils mkdir fails on NFS races)
B=/path/to/TemporalGenomicsEvaluation
PLAINMAP=/path/to/software/PlainMap/plainmap.sh
MAPDAMAGE=/path/to/software/miniconda3/bin/mapDamage
PY=/path/to/software/miniconda3/bin/python
TVCOUNT=$B/tools/TvCount/TvCount.py
export PATH=$PATH:/path/to/software/miniconda3/envs/r_phylo_env/bin   # Rscript for TvCount PDFs
PILOT=1000000
T=${SLURM_CPUS_PER_TASK:-8}

IFS=$'\t' read -r species voucher crr period ref MODE <<< "$(awk -F'\t' -v n=$((SLURM_ARRAY_TASK_ID + 1)) 'NR==n' $B/metadata/pilot_samples.tsv)"
sp=${species// /_}
# prefer the pre-cut 1.1M-read heads (fetch_chunked.sh head); fall back to the full download
reads=$B/pilot_heads/$sp/$crr
ok=headok
[[ -f $reads/${crr}_r1.fq.gz.headok ]] || { reads=$B/raw_reads/$sp/$crr; ok=md5ok; }
prefix=${voucher}_${crr}
outdir=$B/pilot/$MODE/$sp
echo "task $SLURM_ARRAY_TASK_ID: $species $voucher $crr $period mode=$MODE ref=$ref host=$(hostname)"

for f in $reads/${crr}_r1.fq.gz $reads/${crr}_r2.fq.gz; do
    [[ -f $f.$ok ]] || { echo "SKIP: $f not downloaded/verified yet"; exit 0; }
done
[[ -f $ref.bwt ]] || { echo "SKIP: reference $ref not indexed yet"; exit 0; }

mkdir -p $B/pilot/manifests $outdir
manifest=$B/pilot/manifests/$prefix.txt
printf '%s\n' $reads/${crr}_r1.fq.gz $reads/${crr}_r2.fq.gz > $manifest

bash $PLAINMAP \
    --manifest $manifest \
    --prefix $prefix \
    --ref $ref \
    --outdir $outdir \
    --library-type $MODE \
    --threads $T \
    --pilot-fragments $PILOT \
    --run-mapdamage --mapdamage-serial --mapdamage $MAPDAMAGE

bam=$(find $outdir -name "$prefix.final.RG.bam" | head -1)
[[ -n $bam ]] || { echo "ERROR: final BAM not found for $prefix"; exit 1; }
mkdir -p $outdir/tvcount
$PY $TVCOUNT \
    --bam $bam \
    --reference $ref \
    --sample-label "${prefix}_$MODE" \
    --summary-tsv-out $outdir/tvcount/$prefix.summary.tsv \
    --pos-tsv-out $outdir/tvcount/$prefix.positions.tsv \
    --plot-pdf-out $outdir/tvcount/$prefix.pdf \
    --plot-log-y \
    --max-reads 0 \
    --threads $T
echo "done $(date)"
