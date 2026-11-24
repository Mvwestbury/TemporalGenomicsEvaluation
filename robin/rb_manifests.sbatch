#!/usr/bin/env bash
#SBATCH --job-name=rb_manifests
#SBATCH --partition=all
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=/path/to/TemporalGenomicsEvaluation/logs/rb_manifests_%j.out
# After rb_map: the robin manifests, built as for the Liu birds.
#   het_matched/manifest_robin{,_raw,_rescaled}.tsv : rb_depth_match_manifest.py (rank-paired historical/modern,
#                                                     modern subsampled to the partner's depth), split by input
#   hearty/manifest_robin.tsv : rows 1-10 raw PlainMap BAMs, rows 11-20 mapDamage-rescaled BAMs (same sample order),
#                               mean = PlainMap autosomal depth (%.3f), max depth = round(2 x mean)
set -euo pipefail
source /path/to/TemporalGenomicsEvaluation/scripts/safe_mkdir.sh   # race-safe mkdir -p (uutils mkdir fails on NFS races)
B=/path/to/TemporalGenomicsEvaluation
PY=/path/to/software/miniconda3/bin/python
M=$B/het_matched/manifest_robin
$PY $B/robin/rb_depth_match_manifest.py > $M.tsv
for inp in raw rescaled; do awk -F'\t' -v x=$inp 'NR==1 || $4==x' $M.tsv > ${M}_$inp.tsv; done
wc -l $M*.tsv

H=$B/hearty/manifest_robin.tsv
printf 'species\tsample\tperiod\tinput\tbam\tref\tmean_depth\tmax_depth\n' > $H
for inp in raw rescaled; do
    tail -n +2 $B/robin/robin_samples.tsv | while IFS=$'\t' read -r species voucher crr period ref; do
        sp=${species// /_}; prefix=${voucher}_${crr}
        rf=$B/references/regions/$(basename $ref .fasta).autosomes.txt
        cov=$B/mapped/modern/$sp/$prefix.plainmap.coverage.tsv
        mean=$(awk 'NR==FNR{a[$1]=1; next} ($1":") in a{l+=$3; d+=$3*$7} END{printf "%.3f", d/l}' $rf $cov)
        maxd=$(awk -v m=$mean 'BEGIN{print int(2*m + 0.5)}')
        bam=$B/mapped/modern/$sp/$prefix.final.RG.bam
        [[ $inp == raw ]] || bam=$B/mapped/rescaled/$sp/$prefix.rescaled.bam
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' $sp $prefix $period $inp $bam $ref $mean $maxd >> $H
    done
done
cat $H | cut -f2-4,7,8
[[ $(wc -l < $H) -eq 21 ]]
