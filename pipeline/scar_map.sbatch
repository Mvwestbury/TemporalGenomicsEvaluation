#!/usr/bin/env bash
#SBATCH --job-name=scar_map
#SBATCH --partition=all
#SBATCH --cpus-per-task=10
#SBATCH --mem=32G
#SBATCH --output=/path/to/TemporalGenomicsEvaluation/logs/scar_map_%A_%a.out
# Map one SCAR arm with PlainMap (modern mode, SE), then damage QC on 1M reads (mapDamage + TvCount).
# usage: sbatch --array=1-(4*N) scar_map.sbatch <samples.tsv>
#   task -> sample row ((i-1)/4)+1, arm (i-1)%4 of: original frag frag_deam frag_deam_tvex
set -euo pipefail
source /path/to/TemporalGenomicsEvaluation/scripts/safe_mkdir.sh   # race-safe mkdir -p (uutils mkdir fails on NFS races)
B=/path/to/TemporalGenomicsEvaluation
PLAINMAP=/path/to/software/PlainMap/plainmap.sh
MAPDAMAGE=/path/to/software/miniconda3/bin/mapDamage
PY=/path/to/software/miniconda3/bin/python
TVCOUNT=$B/tools/TvCount/TvCount.py
export PATH=$PATH:/path/to/software/miniconda3/envs/r_phylo_env/bin
T=${SLURM_CPUS_PER_TASK:-10}
ARMS=(original frag frag_deam frag_deam_tvex)
i=$SLURM_ARRAY_TASK_ID
arm=${ARMS[$(( (i - 1) % 4 ))]}

IFS=$'\t' read -r species voucher crr ref <<< "$(sed -n $(( (i - 1) / 4 + 1 ))p $1)"
sp=${species// /_}
prefix=${voucher}_${crr}
d=$B/scar_sim/$sp/$prefix
fq=$(ls $d/$arm.f*q.gz)
name=${prefix}_$arm
out=$d/map_$arm
echo "task $i: $name host=$(hostname) $(date)"
mkdir -p $out
echo $fq > $out/manifest.txt

bash $PLAINMAP --manifest $out/manifest.txt --prefix $name --ref $ref --outdir $out \
    --library-type modern --threads $T
bam=$out/$name.final.RG.bam
[[ -f $bam ]] || { echo "ERROR: no BAM"; exit 1; }

# damage QC on 1M reads, comparable with the pilot
$MAPDAMAGE -i $bam -r $ref -d $out/mapdamage --merge-reference-sequences -n 1000000 --downsample-seed 1 --no-stats
$PY $TVCOUNT --bam $bam --reference $ref --sample-label $name \
    --summary-tsv-out $out/tvcount.summary.tsv --pos-tsv-out $out/tvcount.positions.tsv \
    --plot-pdf-out $out/tvcount.pdf --plot-log-y --max-reads 1000000 --threads $T
echo "done $(date)"
