#!/usr/bin/env bash
#SBATCH --job-name=hearty_nat
#SBATCH --partition=all
#SBATCH --cpus-per-task=8
#SBATCH --mem=12G
#SBATCH --output=/path/to/TemporalGenomicsEvaluation/logs/hearty_nat_%A_%a.out
# HEARTY windows + ROH for modern birds at NATIVE depth (no downsampling), plain (non-rescaled) BAMs.
# Reuses the basecall table; min depth = max(5, half mean autosomal depth), max = 2x mean; het threshold 0.25
# (0.05 also written by HEARTY). ROH: ROH_Window_Summary.R defaults (500 kb bins, abs. threshold 1e-4).
set -euo pipefail
source /path/to/TemporalGenomicsEvaluation/scripts/safe_mkdir.sh   # race-safe mkdir -p (uutils mkdir fails on NFS races)
B=/path/to/TemporalGenomicsEvaluation
H=$B/hearty/HEARTY_d8305fe
IFS=$'\t' read -r sp sample period input bam ref mean maxd <<< "$(sed -n $((SLURM_ARRAY_TASK_ID + 1))p $B/hearty/manifest.tsv)"
[[ $period == modern && $input == raw ]] || { echo "not a modern raw row"; exit 1; }
mind=$(python3 -c "print(max(5, round($mean/2)))")
o=$B/hearty/raw/$sp/native; mkdir -p $o || true
echo "task $SLURM_ARRAY_TASK_ID: $sample mean=$mean min=$mind max=$maxd $(date)"
python3 $H/HEARTY.py --basecall-file $B/hearty/raw/$sp/$sample.basecall.txt.gz --out-prefix $sample --outdir $o \
    --min-depth $mind --max-depth $maxd --threshold 0.25 \
    --downstream-workers 6 --progress-interval 600 --force
for t in 0.25; do  # HEARTY writes windows only for -t thresholds (0.05 = minorfreq only)
  w=$(ls $o/${sample}_het$t.mincov${mind}_maxcov${maxd}*.windows.txt.gz | grep -v _trv)
  Rscript $H/ROH_Window_Summary.R --mode single --input $w --output-pdf $o/roh/${sample}_het$t.roh.pdf
done
echo "done $(date)"
