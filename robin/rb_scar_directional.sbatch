#!/usr/bin/env bash
#SBATCH --job-name=rb_scar_dir
#SBATCH --partition=all
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --output=/path/to/TemporalGenomicsEvaluation/logs/rb_scar_dir_%A_%a.out
# Directional singletons/doubletons on the robin SCAR arms: hearty/directional_singletons.py with the settings of
# hearty/directional.sbatch (depth 5, seed 1, 7 workers; min depth 5, max = round(2 x mean) as rb_scar_hearty)
# on the basecall table written by rb_scar_hearty.sbatch. Task i: sample (i-1)/4+1 of robin/robin_scar_samples.tsv, arm (i-1)%4.
set -euo pipefail
source /path/to/TemporalGenomicsEvaluation/scripts/safe_mkdir.sh   # race-safe mkdir -p (uutils mkdir fails on NFS races)
B=/path/to/TemporalGenomicsEvaluation
ARMS=(original frag frag_deam frag_deam_tvex)
i=$SLURM_ARRAY_TASK_ID; arm=${ARMS[$(( (i - 1) % 4 ))]}
IFS=$'\t' read -r species voucher crr ref <<< "$(sed -n $(( (i - 1) / 4 + 1 ))p $B/robin/robin_scar_samples.tsv)"
prefix=${voucher}_${crr}; d=$B/scar_sim/${species// /_}/$prefix
rf=$B/references/regions/$(basename $ref .fasta).autosomes.txt
mean=$(awk 'NR==FNR{a[$1]=1; next} ($1":") in a{l+=$3; d+=$3*$7} END{printf "%.4f", d/l}' $rf $d/map_$arm/${prefix}_$arm.plainmap.coverage.tsv)
maxd=$(awk -v m=$mean 'BEGIN{print int(2*m + 0.5)}')
bc=$d/hearty/$arm/${prefix}_$arm.basecall.txt.gz
[[ -s $bc ]] || { echo "ERROR: no basecall table $bc"; exit 1; }
echo "task $i: $prefix $arm max=$maxd $(date)"
python3 $B/hearty/directional_singletons.py $bc $d/hearty/$arm/${prefix}_$arm.directional.tsv 5 $maxd 5 1 7
head -1 $d/hearty/$arm/${prefix}_$arm.directional.tsv; echo "done $(date)"
