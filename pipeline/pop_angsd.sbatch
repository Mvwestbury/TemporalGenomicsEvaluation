#!/usr/bin/env bash
#SBATCH --job-name=pop_angsd
#SBATCH --partition=all
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --output=/path/to/TemporalGenomicsEvaluation/logs/pop_angsd_%A_%a.out
# Replicates the authors' population-level ANGSD analysis (authors_code/4.diversity ...sh.txt):
# one SAF per species and period with their filters, then Watterson's theta and Tajima's D in
# 50 kb windows (folded SFS) and their "heterozygosity" = singleton class of the unfolded SFS
# (a[2]/sum(a)), polarised by the reference, -GL 1 as in their code.
# Deviation: restricted to the autosome list (chromosome-scale autosomal scaffolds >= 1.2 Mb; -rf),
# as in all our other analyses (their code ran the whole reference, ~600-700k tiny contigs for the babblers).
# Array index 1-8: species (alphabetical) x period (historical, modern).
#   pop_angsd.sbatch            -> per-population steps
#   pop_angsd.sbatch fst        -> index 1-4: Fst historical vs modern for that species
set -euo pipefail
source /path/to/TemporalGenomicsEvaluation/scripts/safe_mkdir.sh   # race-safe mkdir -p (uutils mkdir fails on NFS races)
B=/path/to/TemporalGenomicsEvaluation
A=/path/to/software/angsd-0.935
T=${SLURM_CPUS_PER_TASK:-8}
species=$(awk -F'\t' 'NR>1{print $1}' $B/metadata/full_map_samples.tsv | sort -u | sed -n "$(( (SLURM_ARRAY_TASK_ID - 1) / 2 + 1 ))p")
sp=${species// /_}
ref=$(awk -F'\t' -v s="$species" '$1==s{print $5; exit}' $B/metadata/full_map_samples.tsv)
mkdir -p $B/pop/$sp

if [[ ${1:-} == fst ]]; then
    sp=$(awk -F'\t' 'NR>1{print $1}' $B/metadata/full_map_samples.tsv | sort -u | sed -n "${SLURM_ARRAY_TASK_ID}p"); sp=${sp// /_}
    o=$B/pop/$sp/$sp
    $A/misc/realSFS $o.modern.saf.idx $o.historical.saf.idx -P $T > $o.fst.ml
    $A/misc/realSFS fst index $o.modern.saf.idx $o.historical.saf.idx -sfs $o.fst.ml -fstout $o.fst -P $T
    $A/misc/realSFS fst stats $o.fst.fst.idx > $o.fst.txt
    cat $o.fst.txt; exit 0
fi

period=$( (( SLURM_ARRAY_TASK_ID % 2 )) && echo historical || echo modern )
o=$B/pop/$sp/$sp.$period
awk -F'\t' -v s="$species" -v p=$period -v B=$B '$1==s && $4==p{sp=$1; gsub(/ /,"_",sp); print B"/mapped/modern/"sp"/"$2"_"$3".final.RG.bam"}' \
    $B/metadata/full_map_samples.tsv > $o.bamlist
n=$(wc -l < $o.bamlist)
echo "$species $period: $n BAMs, ref $ref, host $(hostname) $(date)"
while read -r b; do [[ -f $b ]] || { echo "ERROR: missing $b"; exit 1; }; done < $o.bamlist

$A/angsd -b $o.bamlist -GL 1 -doSaf 1 -P $T -remove_bads 1 -only_proper_pairs 1 \
    -minMapQ 30 -minQ 20 -minInd $n -anc $ref -rf $B/references/regions/$(basename $ref .fasta).autosomes.txt -out $o
# theta and Tajima's D (folded)
$A/misc/realSFS $o.saf.idx -P $T -fold 1 > $o.sfs
$A/misc/realSFS saf2theta $o.saf.idx -sfs $o.sfs -outname $o -fold 1
$A/misc/thetaStat print $o.thetas.idx > $o.thetas.txt
$A/misc/thetaStat do_stat $o.thetas.idx -win 50000 -step 20000 -outnames $o.Win_50K
# their "heterozygosity": unfolded SFS, singleton class
$A/misc/realSFS $o.saf.idx -P $T > $o.ml
awk -v s=$sp -v p=$period '{n=0; for(i=1;i<=NF;i++) n+=$i; print s"\t"p"\t"n"\t"$2/n}' $o.ml > $o.H
cat $o.H
echo "done $(date)"
