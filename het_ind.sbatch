#!/usr/bin/env bash
#SBATCH --job-name=het_ind
#SBATCH --partition=all
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=/path/to/TemporalGenomicsEvaluation/logs/het_ind_%A_%a.out
# Per-individual autosomal heterozygosity with ANGSD (single-sample SAF -> folded SFS, H = sfs[2]/sum).
# Array index 1-40: all sites, 41-80: transversions only (-noTrans 1); the sample is
# row ((i-1) % 40) + 1 of metadata/full_map_samples.tsv. Depth filter: mean/2 to 2*mean,
# with the mean autosomal depth taken from PlainMap's coverage table.
set -euo pipefail
source /path/to/TemporalGenomicsEvaluation/scripts/safe_mkdir.sh   # race-safe mkdir -p (uutils mkdir fails on NFS races)
B=/path/to/TemporalGenomicsEvaluation
A=/path/to/software/angsd-0.935
T=${SLURM_CPUS_PER_TASK:-4}
i=$SLURM_ARRAY_TASK_ID
if (( i <= 40 )); then variant=all; extra=""; else variant=tv; extra="-noTrans 1"; fi
row=$(( (i - 1) % 40 + 1 ))

IFS=$'\t' read -r species voucher crr period ref <<< "$(awk -F'\t' -v n=$((row + 1)) 'NR==n' $B/metadata/full_map_samples.tsv)"
sp=${species// /_}
prefix=${voucher}_${crr}
bam=$B/mapped/modern/$sp/$prefix.final.RG.bam
cov=$B/mapped/modern/$sp/$prefix.plainmap.coverage.tsv
rf=$B/references/regions/$(basename $ref .fasta).autosomes.txt
out=$B/het/$variant/$sp/$prefix
echo "task $i: $species $voucher $crr $period variant=$variant host=$(hostname) $(date)"
[[ -f $bam && -f $cov ]] || { echo "SKIP: $bam not mapped yet"; exit 0; }
mkdir -p $(dirname $out)

mean=$(awk 'NR==FNR{a[$1]=1; next} ($1":") in a{l+=$3; d+=$3*$7} END{printf "%.3f", d/l}' $rf $cov)
mind=$(awk -v m=$mean 'BEGIN{x=int(m/2 + 0.5); print (x < 1 ? 1 : x)}')
maxd=$(awk -v m=$mean 'BEGIN{print int(2*m + 0.5)}')
echo "mean autosomal depth $mean -> setMinDepth $mind setMaxDepth $maxd"

$A/angsd -i $bam -ref $ref -anc $ref -rf $rf \
    -doSaf 1 -GL 2 -doCounts 1 \
    -minMapQ 30 -minQ 30 -uniqueOnly 1 -remove_bads 1 \
    -setMinDepth $mind -setMaxDepth $maxd $extra \
    -P $T -out $out
$A/misc/realSFS $out.saf.idx -fold 1 -P $T > $out.sfs
awk -v s=$prefix -v sp=$sp -v p=$period -v v=$variant -v m=$mean -v lo=$mind -v hi=$maxd \
    '{n=$1+$2; print sp"\t"s"\t"p"\t"v"\t"m"\t"lo"\t"hi"\t"n"\t"$2/n}' $out.sfs > $out.het
cat $out.het
echo "done $(date)"
