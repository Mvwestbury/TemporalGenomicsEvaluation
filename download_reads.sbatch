#!/usr/bin/env bash
#SBATCH --job-name=gsa_dl
#SBATCH --partition=all
#SBATCH --cpus-per-task=6
#SBATCH --mem=4G
#SBATCH --output=/path/to/TemporalGenomicsEvaluation/logs/gsa_dl_%j.out
# Usage: sbatch download_reads.sbatch metadata/four_species_runs.tsv [n_parallel]
cd /path/to/TemporalGenomicsEvaluation
echo "host $(hostname)  start $(date)"
curl -s -o /dev/null -w "GSA reachable: HTTP %{http_code}\n" --max-time 60 https://download.cncb.ac.cn/gsa5/
./download_reads.sh "$1" "${2:-6}"
echo "end $(date)"
tail -1 logs/download_*.log
