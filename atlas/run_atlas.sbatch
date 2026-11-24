#!/bin/bash
#SBATCH --job-name=atlas_het
#SBATCH --partition=all
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --output=/path/to/TemporalGenomicsEvaluation/logs/atlas_het_%A_%a.out
# ATLAS heterozygosity for one (sample, depth level). Array index i: sample row (i-1)/5+1 of
# metadata/full_map_samples.tsv, level (i-1)%5 of LEVELS (native 15 10 7.5 5).
# Depth = mean over the evaluation windows (samtools depth -q 20 -Q 30 -s -G 3844), as for the bears.
# Downsampling only where native depth >= target; one samtools subsample per level, distinct seeds.
# Input = PlainMap BAM without mapDamage rescaling (ATLAS fits its own PMD + recalibration per BAM).
set -euo pipefail
source /path/to/TemporalGenomicsEvaluation/scripts/safe_mkdir.sh   # race-safe mkdir -p (uutils mkdir fails on NFS races)
source /path/to/TemporalGenomicsEvaluation/atlas/config.sh
i=$SLURM_ARRAY_TASK_ID
row=$(( (i - 1) / 5 + 1 )); level=${LEVELS[$(( (i - 1) % 5 ))]}
IFS=$'\t' read -r species voucher crr period ref <<< "$(awk -F'\t' -v n=$((row + 1)) 'NR==n' $B/metadata/full_map_samples.tsv)"
sp=${species// /_}; prefix=${voucher}_${crr}; r=$(basename $ref .fasta)
input=$B/mapped/modern/$sp/$prefix.final.RG.bam
C=$ROOT/calib/$r
echo "task $i: $prefix $period level=$level ref=$r host=$(hostname) $(date)"
[[ -f $B/mapped/modern/$sp/$prefix.plainmap.coverage.tsv ]] || { echo "SKIP: not mapped"; exit 0; }
[[ -e $C/CALIBRATION_SUCCESS ]] || { echo "ERROR: no calibration for $r"; exit 1; }

sd=$ROOT/samples/$sp/$prefix; mkdir -p $sd
total=$(awk '{n+=$3-$2} END{print n}' $C/evaluation.bed)
measure() { samtools depth -b $C/evaluation.bed -q 20 -Q 30 -s -G 3844 "$1" | awk -v t=$total '{s+=$3} END{printf "%.6f\n", s/t}'; }
( flock 9; [[ -s $sd/native_depth.txt ]] || measure $input > $sd/native_depth.txt ) 9> $sd/.lock
native=$(<$sd/native_depth.txt)
echo "native window depth $native"

work=$ROOT/results/$sp/$prefix/$level; mkdir -p $work; cd $work
rm -f ANALYSIS_FINISHED FIT_REVIEW_REQUIRED
if [[ $level == native ]]; then
    bam=$input; depth=$native
else
    awk -v d=$native -v t=$level 'BEGIN{exit !(d >= t)}' || { echo "SKIP: native $native < $level"; rm -rf $work; exit 0; }
    frac=$(awk -v d=$native -v t=$level 'BEGIN{printf "%.6f", t/d}')
    seed=$(( 1000 + row * 10 + (i - 1) % 5 ))
    bam=$work/ds.bam
    samtools view -@ 1 -b -F 3844 --subsample $frac --subsample-seed $seed -o $bam $input
    samtools index $bam
    depth=$(measure $bam)
fi
printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$species" "$prefix" "$period" "$level" "$native" "$depth" > depth.tsv
echo "analysis depth $depth"

common=(--bam $bam --fasta $ref --filterMQ 30, --filterBaseQual 20,93
        --filterSoftClips ,0.1 --filterMappingMismatches ,0.05)
"$ATLAS" estimateErrors "${common[@]}" --regions $C/calibration.bed \
    --ploidy 1 --recalModel 'intercept;quality;position;context;fragmentLength;mappingQuality;' \
    --iterations 200 --out calibration > calibration.console.log 2>&1
[[ -s calibration_RGInfo.json ]]
grep -q 'EM has converged' calibration.console.log || echo "Calibration did not report convergence" > FIT_REVIEW_REQUIRED
"$ATLAS" summaryStats "${common[@]}" --RGInfo calibration_RGInfo.json \
    --regions $C/evaluation.bed --mask $C/calibration.bed \
    --estimators HKY85 --genomeWide --out heterozygosity > heterozygosity.console.log 2>&1
"$ATLAS" summaryStats "${common[@]}" --RGInfo calibration_RGInfo.json \
    --regions $C/evaluation.bed --mask $C/calibration.bed \
    --window $C/evaluation.bed --estimators HKY85 --out windows > windows.console.log 2>&1
gzip -t heterozygosity_stats.txt.gz windows_stats.txt.gz
[[ $level == native ]] || rm -f $work/ds.bam $work/ds.bam.bai
touch ANALYSIS_FINISHED
echo "done $(date)"
