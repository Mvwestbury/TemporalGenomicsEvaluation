#!/usr/bin/env bash
#SBATCH --job-name=scar_recal
#SBATCH --partition=all
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/path/to/TemporalGenomicsEvaluation/logs/scar_recal_%A_%a.out
# ANGSD heterozygosity on ATLAS-recalibrated SCAR arms.
# Same reads, depth target, sites and ANGSD settings as scar_het.sbatch; the only change is that base
# qualities come from ATLAS (filterBAM with the arm's own calibration_RGInfo.json from scar_atlas.sbatch).
# Bases with ORIGINAL quality < 30 are set to Q0 (scripts/mask_recal_quals.py) and ANGSD runs with -minQ 1,
# so exactly the bases of the -minQ 30 runs are used. PMD is not incorporated (--incorporatePMD is broken
# in ATLAS 2.0.4: recalibrateWithPMD skips parsing and writes SEQ/QUAL '*').
# usage: sbatch --array=1-12 scar_recal_angsd.sbatch scar_sim/fringilla_samples.tsv
#   task -> sample ((i-1)/4)+1, arm ((i-1)%4); runs all sites and transversions
set -euo pipefail
source /path/to/TemporalGenomicsEvaluation/scripts/safe_mkdir.sh   # race-safe mkdir -p (uutils mkdir fails on NFS races)
B=/path/to/TemporalGenomicsEvaluation
A=/path/to/software/angsd-0.935
PY=/path/to/software/miniconda3/bin/python
source $B/atlas/config.sh
T=${SLURM_CPUS_PER_TASK:-4}
TARGET=6
ARMS=(original frag frag_deam frag_deam_tvex)
i=$SLURM_ARRAY_TASK_ID
arm=${ARMS[$(( (i - 1) % 4 ))]}
IFS=$'\t' read -r species voucher crr ref <<< "$(sed -n $(( (i - 1) / 4 + 1 ))p $B/$1)"
sp=${species// /_}; prefix=${voucher}_${crr}; d=$B/scar_sim/$sp/$prefix
rf=$B/references/regions/$(basename $ref .fasta).autosomes.txt
rg=$d/atlas/$arm/calibration_RGInfo.json
[[ -s $rg ]] || { echo "no ATLAS calibration $rg"; exit 1; }
echo "task $i: $prefix $arm host=$(hostname) $(date)"

depth() { awk 'NR==FNR{a[$1]=1; next} ($1":") in a{l+=$3; d+=$3*$7} END{printf "%.4f", d/l}' $rf $1; }
target=$TARGET
for a in "${ARMS[@]}"; do
    x=$(depth $d/map_$a/${prefix}_$a.plainmap.coverage.tsv)
    target=$(awk -v t=$target -v x=$x 'BEGIN{print (x < t ? x : t)}')
done
bam=$d/map_$arm/${prefix}_$arm.final.RG.bam
mean=$(depth $d/map_$arm/${prefix}_$arm.plainmap.coverage.tsv)
frac=$(awk -v t=$target -v m=$mean 'BEGIN{printf "%.6f", t/m}')
echo "target depth $target, this arm $mean -> keep fraction $frac"

w=$d/het_atlasrecal/work_$arm; rm -rf $w; mkdir -p $w; trap "rm -rf $w" EXIT; cd $w
if awk -v f=$frac 'BEGIN{exit !(f < 0.999)}'; then
    samtools view -b -@ $T -s $(awk -v f=$frac 'BEGIN{printf "%.6f", 42 + f}') -o ds.bam $bam   # same seed as scar_het
else
    samtools view -b -@ $T -o ds.bam $bam
fi
samtools index ds.bam
"$ATLAS" filterBAM --bam ds.bam --RGInfo $rg --out recal > filterBAM.log 2>&1
grep -q 'Will write recalibrated quality scores' filterBAM.log
$PY $B/scripts/mask_recal_quals.py ds.bam recal_filtered.bam masked.bam --minq 30
samtools index masked.bam

mind=$(awk -v m=$target 'BEGIN{x=int(m/2 + 0.5); print (x < 1 ? 1 : x)}')
maxd=$(awk -v m=$target 'BEGIN{print int(2*m + 0.5)}')
for variant in all tv; do
    extra=""; [[ $variant == tv ]] && extra="-noTrans 1"
    o=$d/het_atlasrecal/$variant; mkdir -p $o; out=$o/${prefix}_$arm
    $A/angsd -i masked.bam -ref $ref -anc $ref -rf $rf \
        -doSaf 1 -GL 2 -doCounts 1 \
        -minMapQ 30 -minQ 1 -uniqueOnly 1 -remove_bads 1 \
        -setMinDepth $mind -setMaxDepth $maxd $extra \
        -P $T -out $out
    $A/misc/realSFS $out.saf.idx -fold 1 -P $T > $out.sfs
    H=$(awk '{s=0; for(j=1;j<=NF;j++) s+=$j; printf "%.7g\t%.0f", $2/s, s}' $out.sfs)
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$species" "$prefix" "$arm" "$variant" "$target" "$mind" "$maxd" "$H" > $out.het
    cat $out.het
done
echo "done $(date)"
