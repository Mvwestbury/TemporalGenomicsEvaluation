#!/bin/bash
#SBATCH --job-name=scar_atlas_tq60
#SBATCH --partition=all
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --output=/path/to/TemporalGenomicsEvaluation/logs/scar_atlas_tq60_%A_%a.out
# Top-bin Q60 test of ATLAS on the Liu-style SCAR Fringilla arms (replaces the all-Q60 test, which also gave every sequencing error Q60):
# same cached common target depth, same subsample seed 2000+i, i.e. the same reads), except that each base quality
# of the analysed BAM with Q>=37 (F-I) is set to Q60, lower qualities kept. Output atlas/<arm>_tq60/
# usage: sbatch --array=1-12 scar_atlas_tq60.sbatch
set -euo pipefail
source /path/to/TemporalGenomicsEvaluation/scripts/safe_mkdir.sh   # race-safe mkdir -p (uutils mkdir fails on NFS races)
source /path/to/TemporalGenomicsEvaluation/atlas/config.sh
ARMS=(original frag frag_deam frag_deam_tvex)
i=$SLURM_ARRAY_TASK_ID; arm=${ARMS[$(( (i - 1) % 4 ))]}
IFS=$'\t' read -r species voucher crr ref <<< "$(sed -n $(( (i - 1) / 4 + 1 ))p $B/scar_sim/fringilla_samples.tsv)"
prefix=${voucher}_${crr}; d=$B/scar_sim/${species// /_}/$prefix; r=$(basename $ref .fasta); C=$ROOT/calib/$r
total=$(awk '{n+=$3-$2} END{print n}' $C/evaluation.bed)
measure() { samtools depth -b $C/evaluation.bed -q 20 -Q 30 -s -G 3844 "$1" | awk -v t=$total '{s+=$3} END{printf "%.6f\n", s/t}'; }
for a in "${ARMS[@]}"; do [[ -s $d/atlas/$a.depth ]] || { echo "ERROR: no cached depth for $a"; exit 1; }; done
target=$(cat $d/atlas/{original,frag,frag_deam,frag_deam_tvex}.depth | awk 'BEGIN{m=5} {if ($1 < m) m = $1} END{printf "%.4f", m}')
native=$(<$d/atlas/$arm.depth); input=$d/map_$arm/${prefix}_$arm.final.RG.bam
work=$d/atlas/${arm}_tq60; mkdir -p $work; cd $work; rm -f ANALYSIS_FINISHED FIT_REVIEW_REQUIRED
echo "task $i: $prefix $arm Q60 native=$native target=$target host=$(hostname) $(date)"
frac=$(awk -v d=$native -v t=$target 'BEGIN{printf "%.6f", t/d}')
sub=(); awk -v f=$frac 'BEGIN{exit !(f < 0.999)}' && sub=(--subsample $frac --subsample-seed $((2000 + i)))
# same reads as scar_atlas.sbatch; QUAL Q37-40 (F-I) -> Q60 (]), lower kept
samtools view -@ 1 -h -F 3844 "${sub[@]}" $input \
  | perl -F'\t' -lane 'BEGIN{$,="\t"} $F[10] =~ tr/F-I/]/ unless /^\@/; print @F' \
  | samtools view -@ 1 -b -o ds_tq60.bam -
samtools index ds_tq60.bam; bam=ds_tq60.bam
depth=$(measure $bam); printf "%s\t%s_tq60\t%s\t%s\n" $prefix $arm $native $depth > depth.tsv; echo "analysis depth $depth"
common=(--bam $bam --fasta $ref --filterMQ 30, --filterBaseQual 20,93 --filterSoftClips ,0.1 --filterMappingMismatches ,0.05)
"$ATLAS" estimateErrors "${common[@]}" --regions $C/calibration.bed --ploidy 1 \
    --recalModel 'intercept;quality;position;context;fragmentLength;mappingQuality;' \
    --iterations 200 --out calibration > calibration.console.log 2>&1
grep -q 'EM has converged' calibration.console.log || echo "Calibration did not report convergence" > FIT_REVIEW_REQUIRED
"$ATLAS" summaryStats "${common[@]}" --RGInfo calibration_RGInfo.json --regions $C/evaluation.bed \
    --mask $C/calibration.bed --estimators HKY85 --genomeWide --out heterozygosity > heterozygosity.console.log 2>&1
gzip -t heterozygosity_stats.txt.gz; rm -f ds_tq60.bam ds_tq60.bam.bai; touch ANALYSIS_FINISHED; echo "done $(date)"
