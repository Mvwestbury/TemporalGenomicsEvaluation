#!/usr/bin/env bash
#SBATCH --job-name=gsa_fetch
#SBATCH --partition=all
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --output=/path/to/TemporalGenomicsEvaluation/logs/gsa_fetch_%j.out
# Usage: sbatch fetch_chunked.sbatch FILE_LIST full|head [N_PARALLEL] [CRR_FILTER]
# Keep N_PARALLEL modest: GSA blocks hosts that open too many connections.
cd /path/to/TemporalGenomicsEvaluation
echo "host $(hostname)  start $(date)"
./pipeline/fetch_chunked.sh "$1" "$2" "${3:-6}" "${4:-}"
