#!/usr/bin/env bash
#SBATCH --job-name=scar_het
#SBATCH --partition=all
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=/path/to/TemporalGenomicsEvaluation/logs/scar_het_%A_%a.out
# Per-individual heterozygosity for the SCAR arms, same ANGSD settings as het_ind.sbatch.
# All 4 arms of a sample are downsampled to the same autosomal depth: the lowest arm depth
# of that sample, capped at TARGET (default 6x, roughly the historical Fringilla depth).
# usage: sbatch --array=1-(8*N) scar_het.sbatch <samples.tsv> [TARGET] [GL]
#   task -> sample ((i-1)/8)+1, arm ((i-1)%4), variant all (first 4) / tv (last 4)
set -euo pipefail
source /path/to/TemporalGenomicsEvaluation/scripts/safe_mkdir.sh   # race-safe mkdir -p (uutils mkdir fails on NFS races)
B=/path/to/TemporalGenomicsEvaluation
A=/path/to/software/angsd-0.935
T=${SLURM_CPUS_PER_TASK:-4}
TARGET=${2:-6}
GL=${3:-2}   # ANGSD genotype-likelihood model; outputs for GL != 2 go to het_gl$GL/
ARMS=(original frag frag_deam frag_deam_tvex)
i=$SLURM_ARRAY_TASK_ID
k=$(( (i - 1) % 8 ))
arm=${ARMS[$(( k % 4 ))]}
if (( k < 4 )); then variant=all; extra=""; else variant=tv; extra="-noTrans 1"; fi

IFS=$'\t' read -r species voucher crr ref <<< "$(sed -n $(( (i - 1) / 8 + 1 ))p $1)"
sp=${species// /_}
prefix=${voucher}_${crr}
d=$B/scar_sim/$sp/$prefix
rf=$B/references/regions/$(basename $ref .fasta).autosomes.txt
echo "task $i: $prefix $arm $variant host=$(hostname) $(date)"

depth() { awk 'NR==FNR{a[$1]=1; next} ($1":") in a{l+=$3; d+=$3*$7} END{printf "%.4f", d/l}' $rf $1; }
target=$TARGET
for a in "${ARMS[@]}"; do
    x=$(depth $d/map_$a/${prefix}_$a.plainmap.coverage.tsv)
    echo "  $a autosomal depth $x"
    target=$(awk -v t=$target -v x=$x 'BEGIN{print (x < t ? x : t)}')
done
bam=$d/map_$arm/${prefix}_$arm.final.RG.bam
mean=$(depth $d/map_$arm/${prefix}_$arm.plainmap.coverage.tsv)
frac=$(awk -v t=$target -v m=$mean 'BEGIN{printf "%.6f", t/m}')
echo "target depth $target, this arm $mean -> keep fraction $frac"

hd=het; [[ $GL == 2 ]] || hd=het_gl$GL
o=$d/$hd/$variant
mkdir -p $o
ds=$o/${prefix}_$arm.ds.bam
if awk -v f=$frac 'BEGIN{exit !(f < 0.999)}'; then
    samtools view -b -@ $T -s $(awk -v f=$frac 'BEGIN{printf "%.6f", 42 + f}') -o $ds $bam
else
    ln -sf $bam $ds
fi
samtools index $ds

mind=$(awk -v m=$target 'BEGIN{x=int(m/2 + 0.5); print (x < 1 ? 1 : x)}')
maxd=$(awk -v m=$target 'BEGIN{print int(2*m + 0.5)}')
out=$o/${prefix}_$arm
$A/angsd -i $ds -ref $ref -anc $ref -rf $rf \
    -doSaf 1 -GL $GL -doCounts 1 \
    -minMapQ 30 -minQ 30 -uniqueOnly 1 -remove_bads 1 \
    -setMinDepth $mind -setMaxDepth $maxd $extra \
    -P $T -out $out
$A/misc/realSFS $out.saf.idx -fold 1 -P $T > $out.sfs
H=$(awk '{s=0; for(j=1;j<=NF;j++) s+=$j; printf "%.7g\t%.0f", $2/s, s}' $out.sfs)
printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$species" "$prefix" "$arm" "$variant" "$target" "$mind" "$maxd" "$H" > $out.het
rm -f $ds $ds.bai
cat $out.het
echo "done $(date)"
