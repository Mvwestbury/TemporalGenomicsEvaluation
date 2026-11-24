#!/usr/bin/env bash
#SBATCH --job-name=rb_recal
#SBATCH --partition=all
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/path/to/TemporalGenomicsEvaluation/logs/rb_recal_%A_%a.out
# ANGSD heterozygosity on ATLAS-recalibrated empirical robin BAMs (all sites + transversions).
# = rb_het_ind.sbatch (same BAM, native depth, depth window mean/2..2*mean, same ANGSD settings) with base
#   qualities from ATLAS: filterBAM with the bird's native-level calibration_RGInfo.json (rb_atlas.sbatch),
#   as in scar_recal_angsd.sbatch. Bases with ORIGINAL quality < 30 are set to Q0 (scripts/mask_recal_quals.py)
#   and ANGSD runs with -minQ 1, so exactly the bases of the -minQ 30 runs are used. PMD not incorporated.
#   sbatch --array=1-10 rb_recal_angsd.sbatch      (task = row of robin/robin_samples.tsv)
set -euo pipefail
source /path/to/TemporalGenomicsEvaluation/scripts/safe_mkdir.sh   # race-safe mkdir -p (uutils mkdir fails on NFS races)
B=/path/to/TemporalGenomicsEvaluation
A=/path/to/software/angsd-0.935
PY=/path/to/software/miniconda3/bin/python
source $B/atlas/config.sh
T=${SLURM_CPUS_PER_TASK:-4}

IFS=$'\t' read -r species voucher crr period ref <<< "$(awk -F'\t' -v n=$((SLURM_ARRAY_TASK_ID + 1)) 'NR==n' $B/robin/robin_samples.tsv)"
sp=${species// /_}; prefix=${voucher}_${crr}
bam=$B/mapped/modern/$sp/$prefix.final.RG.bam
cov=$B/mapped/modern/$sp/$prefix.plainmap.coverage.tsv
rf=$B/references/regions/$(basename $ref .fasta).autosomes.txt
rg=$B/atlas/results/$sp/$prefix/native/calibration_RGInfo.json
[[ -s $rg ]] || { echo "ERROR: no ATLAS calibration $rg"; exit 1; }
echo "task $SLURM_ARRAY_TASK_ID: $prefix $period host=$(hostname) $(date)"

mean=$(awk 'NR==FNR{a[$1]=1; next} ($1":") in a{l+=$3; d+=$3*$7} END{printf "%.3f", d/l}' $rf $cov)
mind=$(awk -v m=$mean 'BEGIN{x=int(m/2 + 0.5); print (x < 1 ? 1 : x)}')
maxd=$(awk -v m=$mean 'BEGIN{print int(2*m + 0.5)}')
echo "mean autosomal depth $mean -> setMinDepth $mind setMaxDepth $maxd"

w=$B/het_atlasrecal/work/$prefix; rm -rf $w; mkdir -p $w; trap "rm -rf $w" EXIT; cd $w
"$ATLAS" filterBAM --bam $bam --RGInfo $rg --out recal > filterBAM.log 2>&1
grep -q 'Will write recalibrated quality scores' filterBAM.log
$PY $B/scripts/mask_recal_quals.py $bam recal_filtered.bam masked.bam --minq 30
samtools index masked.bam

for variant in all tv; do
    extra=""; [[ $variant == tv ]] && extra="-noTrans 1"
    o=$B/het_atlasrecal/$variant/$sp; mkdir -p $o; out=$o/$prefix
    $A/angsd -i masked.bam -ref $ref -anc $ref -rf $rf \
        -doSaf 1 -GL 2 -doCounts 1 \
        -minMapQ 30 -minQ 1 -uniqueOnly 1 -remove_bads 1 \
        -setMinDepth $mind -setMaxDepth $maxd $extra \
        -P $T -out $out
    $A/misc/realSFS $out.saf.idx -fold 1 -P $T > $out.sfs
    awk -v s=$prefix -v sp=$sp -v p=$period -v v=$variant -v m=$mean -v lo=$mind -v hi=$maxd \
        '{n=$1+$2; print sp"\t"s"\t"p"\t"v"\t"m"\t"lo"\t"hi"\t"n"\t"$2/n}' $out.sfs > $out.het
    cat $out.het
done
echo "done $(date)"
