#!/usr/bin/env bash
#SBATCH --job-name=rb_scar_hearty
#SBATCH --partition=all
#SBATCH --cpus-per-task=8
#SBATCH --mem=12G
#SBATCH --output=/path/to/TemporalGenomicsEvaluation/logs/rb_scar_hearty_%A_%a.out
# Robin copy of scar_hearty.sbatch (sample table robin/robin_scar_samples.tsv).
# HEARTY on the SCAR arms (positive control with known truth = the unmodified arm).
# Full arm BAMs (HEARTY's fixed 5-read draw makes depth irrelevant), same settings as hearty.sbatch,
# exact mode, pinned HEARTY d8305fe. Task i: sample (i-1)/4+1 of scar_sim/fringilla_samples.tsv, arm (i-1)%4.
set -euo pipefail
source /path/to/TemporalGenomicsEvaluation/scripts/safe_mkdir.sh   # race-safe mkdir -p (uutils mkdir fails on NFS races)
B=/path/to/TemporalGenomicsEvaluation
H=$B/hearty/HEARTY_d8305fe/HEARTY.py
ARMS=(original frag frag_deam frag_deam_tvex)
i=$SLURM_ARRAY_TASK_ID; arm=${ARMS[$(( (i - 1) % 4 ))]}
IFS=$'\t' read -r species voucher crr ref <<< "$(sed -n $(( (i - 1) / 4 + 1 ))p $B/robin/robin_scar_samples.tsv)"
prefix=${voucher}_${crr}; d=$B/scar_sim/${species// /_}/$prefix
rf=$B/references/regions/$(basename $ref .fasta).autosomes.txt
mean=$(awk 'NR==FNR{a[$1]=1; next} ($1":") in a{l+=$3; d+=$3*$7} END{printf "%.4f", d/l}' $rf $d/map_$arm/${prefix}_$arm.plainmap.coverage.tsv)
maxd=$(awk -v m=$mean 'BEGIN{print int(2*m + 0.5)}')
echo "task $i: $prefix $arm mean=$mean max=$maxd $(date)"
python3 $H --bam $d/map_$arm/${prefix}_$arm.final.RG.bam --out-prefix ${prefix}_$arm --outdir $d/hearty/$arm \
    --ref $ref --rf $rf --angsd /path/to/software/angsd-0.935/angsd --angsd-params "-minq 30 -minmapq 30" \
    --min-depth 5 --max-depth $maxd --threshold 0.25 --downsample-depth 5 \
    --cores 8 --downstream-workers 6 --progress-interval 600 --reuse-existing
echo "done $(date)"
