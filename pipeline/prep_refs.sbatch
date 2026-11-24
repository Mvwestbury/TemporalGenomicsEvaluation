#!/usr/bin/env bash
#SBATCH --job-name=bwa_index
#SBATCH --partition=all
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=/path/to/TemporalGenomicsEvaluation/logs/bwa_index_%j.out
# bwa index (and faidx if missing) every reference in references/full, in parallel.
set -euo pipefail
source /path/to/TemporalGenomicsEvaluation/scripts/safe_mkdir.sh   # race-safe mkdir -p (uutils mkdir fails on NFS races)
cd /path/to/TemporalGenomicsEvaluation/references/full
for fa in *.fasta; do
    [[ -f $fa.fai ]] || samtools faidx $fa
    [[ -f $fa.bwt ]] || { bwa index $fa 2> $fa.bwa_index.log && echo "indexed $fa"; } &
done
wait
echo "done $(date)"
