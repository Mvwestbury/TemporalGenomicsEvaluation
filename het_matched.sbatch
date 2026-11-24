#!/usr/bin/env bash
#SBATCH --job-name=het_match
#SBATCH --partition=all
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=/path/to/TemporalGenomicsEvaluation/logs/het_match_%A_%a.out
# Per-individual ANGSD heterozygosity on depth-matched / rescaled BAMs (settings as het_ind.sbatch).
# usage: sbatch --array=1-(2N) het_matched.sbatch <manifest.tsv>   (N data rows; i<=N all sites, i>N -noTrans 1)
# Depth filter: half to twice the target (post-subsampling) depth.
set -euo pipefail
source /path/to/TemporalGenomicsEvaluation/scripts/safe_mkdir.sh   # race-safe mkdir -p (uutils mkdir fails on NFS races)
B=/path/to/TemporalGenomicsEvaluation
A=/path/to/software/angsd-0.935
T=${SLURM_CPUS_PER_TASK:-4}
M=$1
N=$(( $(grep -vc '^#' $M) - 1 ))
i=$SLURM_ARRAY_TASK_ID
if (( i <= N )); then variant=all; extra=""; row=$i; else variant=tv; extra="-noTrans 1"; row=$(( i - N )); fi
IFS=$'\t' read -r species sample period input bam ref native frac target partner seed <<< "$(grep -v '^#' $M | sed -n $((row + 1))p)"
sp=${species// /_}
rf=$B/references/regions/$(basename $ref .fasta).autosomes.txt
o=$B/het_matched/$input/$variant/$sp; mkdir -p $o
out=$o/$sample
echo "task $i: $sample $period $input $variant native=$native frac=$frac target=$target partner=$partner host=$(hostname)"
[[ -s $bam ]] || { echo "SKIP: $bam missing"; exit 0; }

if awk -v f=$frac 'BEGIN{exit !(f < 0.999)}'; then
    ds=$o/$sample.ds.bam
    samtools view -b -@ $T --subsample $frac --subsample-seed $seed -o $ds $bam
    samtools index $ds
else
    ds=$bam
fi
mind=$(awk -v m=$target 'BEGIN{x=int(m/2 + 0.5); print (x < 1 ? 1 : x)}')
maxd=$(awk -v m=$target 'BEGIN{print int(2*m + 0.5)}')
$A/angsd -i $ds -ref $ref -anc $ref -rf $rf \
    -doSaf 1 -GL 2 -doCounts 1 \
    -minMapQ 30 -minQ 30 -uniqueOnly 1 -remove_bads 1 \
    -setMinDepth $mind -setMaxDepth $maxd $extra \
    -P $T -out $out
$A/misc/realSFS $out.saf.idx -fold 1 -P $T > $out.sfs
H=$(awk '{s=0; for(j=1;j<=NF;j++) s+=$j; printf "%.7g\t%.0f", $2/s, s}' $out.sfs)
printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$species" "$sample" "$period" "$input" "$variant" "$target" "$partner" "$mind" "$maxd" "$H" > $out.het
[[ $ds == $bam ]] || rm -f $ds $ds.bai
cat $out.het
