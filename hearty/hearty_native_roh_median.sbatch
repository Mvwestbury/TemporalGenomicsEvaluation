#!/usr/bin/env bash
#SBATCH --job-name=roh_med
#SBATCH --partition=all
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=/path/to/TemporalGenomicsEvaluation/logs/roh_med_%A_%a.out
# ROH on the native-depth HEARTY windows (hearty_native_roh.sbatch) with a relative threshold:
# window is ROH-like if het <= 0.3 x the sample's median window heterozygosity. Other settings default.
set -euo pipefail
source /path/to/TemporalGenomicsEvaluation/scripts/safe_mkdir.sh   # race-safe mkdir -p (uutils mkdir fails on NFS races)
B=/path/to/TemporalGenomicsEvaluation
H=$B/hearty/HEARTY_d8305fe
IFS=$'\t' read -r sp sample period input bam ref mean maxd <<< "$(sed -n $((SLURM_ARRAY_TASK_ID + 1))p $B/hearty/manifest.tsv)"
mind=$(python3 -c "print(max(5, round($mean/2)))")
o=$B/hearty/raw/$sp/native
for t in 0.25; do  # HEARTY writes windows only for -t thresholds (0.05 = minorfreq only)
  w=$(ls $o/${sample}_het$t.mincov${mind}_maxcov${maxd}*.windows.txt.gz | grep -v _trv)
  Rscript $H/ROH_Window_Summary.R --mode single --input $w --output-pdf $o/roh/${sample}_het$t.roh_median0.3.pdf \
      --roh-threshold 0.3 --roh-threshold-mode median
done
