#!/usr/bin/env bash
#SBATCH --job-name=rb_damage
#SBATCH --partition=all
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=/path/to/TemporalGenomicsEvaluation/logs/rb_damage_%A_%a.out
# Damage levels of every robin (for the composite figure; same metrics as pilot/pilot_damage_table.tsv):
# mapDamage + TvCount on 1M reads of the bwa mem BAM, with the settings of rb_profile.sbatch.
#   sbatch --array=1-10 rb_damage.sbatch      (task = row of robin/robin_samples.tsv) -> robin/damage/<prefix>/
set -euo pipefail
source /path/to/TemporalGenomicsEvaluation/scripts/safe_mkdir.sh   # race-safe mkdir -p (uutils mkdir fails on NFS races)
B=/path/to/TemporalGenomicsEvaluation
MAPDAMAGE=/path/to/software/miniconda3/bin/mapDamage
PY=/path/to/software/miniconda3/bin/python
TVCOUNT=$B/tools/TvCount/TvCount.py
export PATH=$PATH:/path/to/software/miniconda3/envs/r_phylo_env/bin
T=${SLURM_CPUS_PER_TASK:-8}
IFS=$'\t' read -r species voucher crr period ref <<< "$(awk -F'\t' -v n=$((SLURM_ARRAY_TASK_ID + 1)) 'NR==n' $B/robin/robin_samples.tsv)"
s=${voucher}_${crr}
bam=$B/mapped/modern/${species// /_}/$s.final.RG.bam
q=$B/robin/damage/$s; mkdir -p $q
echo "task $SLURM_ARRAY_TASK_ID: $s $period $(date)"
$MAPDAMAGE -i $bam -r $ref -d $q/mapdamage --merge-reference-sequences -n 1000000 --downsample-seed 1 --no-stats
$PY $TVCOUNT --bam $bam --reference $ref --sample-label $s \
    --summary-tsv-out $q/tvcount.summary.tsv --pos-tsv-out $q/tvcount.positions.tsv \
    --plot-pdf-out $q/tvcount.pdf --plot-log-y --max-reads 1000000 --threads $T
echo "done $(date)"
