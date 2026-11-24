#!/usr/bin/env bash
#SBATCH --job-name=rb_hearty_dir
#SBATCH --partition=all
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --output=/path/to/TemporalGenomicsEvaluation/logs/rb_hearty_dir_%A_%a.out
# Robin copy of hearty/directional.sbatch (manifest hearty/manifest_robin.tsv: rows 1-10 raw, 11-20 rescaled).
# Reference-directional singleton/doubleton spectrum (HEARTY draw: depth 5, 1 rep, seed 1), row TASK+1 of manifest.
set -euo pipefail
source /path/to/TemporalGenomicsEvaluation/scripts/safe_mkdir.sh   # race-safe mkdir -p (uutils mkdir fails on NFS races)
B=/path/to/TemporalGenomicsEvaluation
IFS=$'\t' read -r sp sample period input bam ref mean maxd <<< "$(sed -n $((SLURM_ARRAY_TASK_ID + 1))p $B/hearty/manifest_robin.tsv)"
o=$B/hearty/$input/$sp
echo "task $SLURM_ARRAY_TASK_ID: $sample $input $(date)"
python3 $B/hearty/directional_singletons.py $o/$sample.basecall.txt.gz $o/ds5_r1_s1/$sample.directional.tsv 5 $maxd 5 1 7
head -1 $o/ds5_r1_s1/$sample.directional.tsv; echo "done $(date)"
