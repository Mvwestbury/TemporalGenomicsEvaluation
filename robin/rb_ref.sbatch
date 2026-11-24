#!/usr/bin/env bash
#SBATCH --job-name=rb_ref
#SBATCH --partition=all
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=/path/to/TemporalGenomicsEvaluation/logs/rb_ref_%j.out
# Black robin reference for the pipeline used for the Liu et al. birds (as prep_refs.sbatch):
#   NCBI GCA_025920805.1 (Ptraversi_NRM_v1, scaffold-level) -> references/full/Ptraversi_NRM_v1.fasta
#   (headers cut to the accession), samtools faidx, bwa index.
# Autosomes = the 99 autosomal scaffolds >= 1 Mb (sex-linked scaffolds removed) from the GT_damage robin
# reference folder -> references/regions/Ptraversi_NRM_v1.autosomes.txt ("name:" per line, as the Liu refs).
set -euo pipefail
source /path/to/TemporalGenomicsEvaluation/scripts/safe_mkdir.sh   # race-safe mkdir -p (uutils mkdir fails on NFS races)
B=/path/to/TemporalGenomicsEvaluation
G=/path/to/GT_damage/robin/reference
fa=$B/references/full/Ptraversi_NRM_v1.fasta
rf=$B/references/regions/Ptraversi_NRM_v1.autosomes.txt
mkdir -p $B/robin/reference_info
cp $G/autosomes_scaffolds.txt $G/Ptraversi_sex_linked_scaffolds.txt $B/robin/reference_info/

if [[ ! -s $fa.fai ]]; then
    zcat $G/GCA_025920805.1_Ptraversi_NRM_v1_genomic.fna.gz | awk '/^>/{print $1; next} {print}' > $fa.tmp
    mv $fa.tmp $fa
    samtools faidx $fa
fi
awk '{print $1":"}' $G/autosomes_scaffolds.txt > $rf
# every autosome must exist in the reference, and none may be sex-linked
awk 'NR==FNR{a[$1]=1; next} !($1 in a){print "missing: "$1; bad=1} END{exit bad}' $fa.fai $G/autosomes_scaffolds.txt
awk 'NR==FNR{s[$1]=1; next} ($1 in s){print "sex-linked: "$1; bad=1} END{exit bad}' $G/Ptraversi_sex_linked_scaffolds.txt $G/autosomes_scaffolds.txt
awk 'NR==FNR{a[$1]=1; next} ($1 in a){n++; l+=$2} END{printf "%d autosomal scaffolds, %.1f Mb\n", n, l/1e6}' $G/autosomes_scaffolds.txt $fa.fai
awk '{l+=$2} END{printf "%d sequences, %.1f Mb total\n", NR, l/1e6}' $fa.fai

[[ -f $fa.bwt ]] || bwa index $fa 2> $fa.bwa_index.log
echo "done $(date)"
