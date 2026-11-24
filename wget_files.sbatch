#!/usr/bin/env bash
#SBATCH --job-name=gsa_wget
#SBATCH --partition=all
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --output=/path/to/TemporalGenomicsEvaluation/logs/gsa_wget_%j.out
# Usage: sbatch wget_files.sbatch FILE_LIST [N_PARALLEL]
# Keep N_PARALLEL modest: GSA blocks hosts that open too many connections.
cd /path/to/TemporalGenomicsEvaluation
echo "host $(hostname)  start $(date)"
./wget_files.sh "$1" "${2:-3}"
