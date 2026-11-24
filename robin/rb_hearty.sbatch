#!/usr/bin/env bash
#SBATCH --job-name=rb_hearty
#SBATCH --partition=all
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=/path/to/TemporalGenomicsEvaluation/logs/rb_hearty_%A_%a.out
# Robin copy of hearty/hearty.sbatch (manifest hearty/manifest_robin.tsv: rows 1-10 raw, 11-20 rescaled).
# HEARTY per BAM (row TASK+1 of hearty/manifest_robin.tsv; raw and mapDamage-rescaled inputs).
# ANGSD counts: autosomes, -minq 30 -minmapq 30 (as the ANGSD het runs), reference bases annotated.
# Downstream: min depth 5, max 2x mean autosomal depth, each site downsampled to 5 reads,
# het threshold 0.25 (= >=2 minor reads of 5), from the same basecall table:
#   ds5_r1_s{1,2,3} : one draw, seeds 1-3
#   ds5_r10_s1      : 10 draws averaged (HEARTY default), for comparison
set -euo pipefail
source /path/to/TemporalGenomicsEvaluation/scripts/safe_mkdir.sh   # race-safe mkdir -p (uutils mkdir fails on NFS races)
B=/path/to/TemporalGenomicsEvaluation
HEARTY=/path/to/HEARTY/HEARTY.py
ANGSD=/path/to/software/angsd-0.935/angsd
T=${SLURM_CPUS_PER_TASK:-8}
IFS=$'\t' read -r sp sample period input bam ref mean maxd <<< "$(sed -n $((SLURM_ARRAY_TASK_ID + 1))p $B/hearty/manifest_robin.tsv)"
rf=$B/references/regions/$(basename $ref .fasta).autosomes.txt
o=$B/hearty/$input/$sp; mkdir -p $o
echo "task $SLURM_ARRAY_TASK_ID: $sample $period $input mean=$mean max=$maxd host=$(hostname) $(date)"

# 1. ANGSD counts -> basecall table (+ ds5_r1_s1 summaries)
python3 $HEARTY --bam $bam --out-prefix $sample --outdir $o/ds5_r1_s1 \
    --ref $ref --rf $rf --angsd $ANGSD --angsd-params "-minq 30 -minmapq 30" \
    --min-depth 5 --max-depth $maxd --threshold 0.25 \
    --downsample-depth 5 --downsample-reps 1 --downsample-seed 1 \
    --cores $T --downstream-workers $((T - 2)) --progress-interval 600 --reuse-existing
bc=$o/ds5_r1_s1/$sample.basecall.txt.gz
[[ -s $bc ]]
mv $bc $o/$sample.basecall.txt.gz && ln -sf ../$sample.basecall.txt.gz $bc
rm -f $o/ds5_r1_s1/$sample.pos.gz $o/ds5_r1_s1/$sample.counts.gz

# 2. further summaries from the same table
run() {  # dir reps seed
    python3 $HEARTY --basecall-file $o/$sample.basecall.txt.gz --out-prefix $sample --outdir $o/$1 \
        --min-depth 5 --max-depth $maxd --threshold 0.25 \
        --downsample-depth 5 --downsample-reps $2 --downsample-seed $3 \
        --downstream-workers $((T - 2)) --progress-interval 600 --force
}
run ds5_r1_s2 1 2
run ds5_r1_s3 1 3
run ds5_r10_s1 10 1
echo "done $(date)"
