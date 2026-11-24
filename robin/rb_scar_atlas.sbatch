#!/bin/bash
#SBATCH --job-name=rb_scar_atlas
#SBATCH --partition=all
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --output=/path/to/TemporalGenomicsEvaluation/logs/rb_scar_atlas_%A_%a.out
# Robin copy of scar_atlas.sbatch (sample table robin/robin_scar_samples.tsv).
# ATLAS on the SCAR arms, same settings as atlas/run_atlas.sbatch. All 4 arms of a bird are
# downsampled to one common evaluation-window depth: min(5x, lowest arm depth). Task i as scar_hearty.sbatch.
set -euo pipefail
source /path/to/TemporalGenomicsEvaluation/scripts/safe_mkdir.sh   # race-safe mkdir -p (uutils mkdir fails on NFS races)
source /path/to/TemporalGenomicsEvaluation/atlas/config.sh
ARMS=(original frag frag_deam frag_deam_tvex)
i=$SLURM_ARRAY_TASK_ID; arm=${ARMS[$(( (i - 1) % 4 ))]}
IFS=$'\t' read -r species voucher crr ref <<< "$(sed -n $(( (i - 1) / 4 + 1 ))p $B/robin/robin_scar_samples.tsv)"
prefix=${voucher}_${crr}; d=$B/scar_sim/${species// /_}/$prefix; r=$(basename $ref .fasta); C=$ROOT/calib/$r
total=$(awk '{n+=$3-$2} END{print n}' $C/evaluation.bed)
measure() { samtools depth -b $C/evaluation.bed -q 20 -Q 30 -s -G 3844 "$1" | awk -v t=$total '{s+=$3} END{printf "%.6f\n", s/t}'; }
mkdir -p $d/atlas
for a in "${ARMS[@]}"; do   # window depth of every arm (cached; arrays of the same bird share it)
  ( flock 9; [[ -s $d/atlas/$a.depth ]] || measure $d/map_$a/${prefix}_$a.final.RG.bam > $d/atlas/$a.depth ) 9> $d/atlas/.lock_$a
done
target=$(cat $d/atlas/*.depth | awk 'BEGIN{m=5} {if ($1 < m) m = $1} END{printf "%.4f", m}')
native=$(<$d/atlas/$arm.depth); input=$d/map_$arm/${prefix}_$arm.final.RG.bam
work=$d/atlas/$arm; mkdir -p $work; cd $work; rm -f ANALYSIS_FINISHED FIT_REVIEW_REQUIRED
echo "task $i: $prefix $arm native=$native target=$target $(date)"
frac=$(awk -v d=$native -v t=$target 'BEGIN{printf "%.6f", t/d}')
if awk -v f=$frac 'BEGIN{exit !(f < 0.999)}'; then
  samtools view -@ 1 -b -F 3844 --subsample $frac --subsample-seed $((2000 + i)) -o ds.bam $input; samtools index ds.bam; bam=ds.bam
else bam=$input; fi
depth=$(measure $bam); printf "%s\t%s\t%s\t%s\n" $prefix $arm $native $depth > depth.tsv; echo "analysis depth $depth"
common=(--bam $bam --fasta $ref --filterMQ 30, --filterBaseQual 20,93 --filterSoftClips ,0.1 --filterMappingMismatches ,0.05)
"$ATLAS" estimateErrors "${common[@]}" --regions $C/calibration.bed --ploidy 1 \
    --recalModel 'intercept;quality;position;context;fragmentLength;mappingQuality;' \
    --iterations 200 --out calibration > calibration.console.log 2>&1
grep -q 'EM has converged' calibration.console.log || echo "Calibration did not report convergence" > FIT_REVIEW_REQUIRED
"$ATLAS" summaryStats "${common[@]}" --RGInfo calibration_RGInfo.json --regions $C/evaluation.bed \
    --mask $C/calibration.bed --estimators HKY85 --genomeWide --out heterozygosity > heterozygosity.console.log 2>&1
gzip -t heterozygosity_stats.txt.gz; rm -f ds.bam ds.bam.bai; touch ANALYSIS_FINISHED; echo "done $(date)"
