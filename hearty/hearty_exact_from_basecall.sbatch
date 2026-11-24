#!/usr/bin/env bash
#SBATCH --job-name=hearty_ex
#SBATCH --partition=all
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=/path/to/TemporalGenomicsEvaluation/logs/hearty_ex_%A_%a.out
# Exact-mode HEARTY summaries (pinned HEARTY d8305fe) from an existing basecall table; row TASK+1 of hearty/manifest.tsv.
set -euo pipefail
source /path/to/TemporalGenomicsEvaluation/scripts/safe_mkdir.sh   # race-safe mkdir -p (uutils mkdir fails on NFS races)
B=/path/to/TemporalGenomicsEvaluation
HEARTY=$B/hearty/HEARTY_d8305fe/HEARTY.py
T=${SLURM_CPUS_PER_TASK:-8}
IFS=$'\t' read -r sp sample period input bam ref mean maxd <<< "$(sed -n $((SLURM_ARRAY_TASK_ID + 1))p $B/hearty/manifest.tsv)"
o=$B/hearty/$input/$sp
bc=$o/$sample.basecall.txt.gz
[[ -s $bc ]] || { [[ -s $o/ds5_r1_s1/$sample.basecall.txt.gz && ! -L $o/ds5_r1_s1/$sample.basecall.txt.gz ]] && mv $o/ds5_r1_s1/$sample.basecall.txt.gz $bc; }
[[ -s $bc ]]
rm -f $o/ds5_r1_s1/$sample.pos.gz $o/ds5_r1_s1/$sample.counts.gz
echo "task $SLURM_ARRAY_TASK_ID: $sample $period $input max=$maxd $(date)"
python3 $HEARTY --basecall-file $bc --out-prefix $sample --outdir $o/ds5_exact \
    --min-depth 5 --max-depth $maxd --threshold 0.25 --downsample-depth 5 --downsample-mode exact \
    --downstream-workers $((T - 2)) --progress-interval 600 --force
echo "done $(date)"
