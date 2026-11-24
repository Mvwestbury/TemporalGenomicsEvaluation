#!/usr/bin/env bash
#SBATCH --job-name=bcf_call
#SBATCH --partition=all
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --output=/path/to/TemporalGenomicsEvaluation/logs/bcf_call_%A_%a.out
# Joint calling of the 5 modern birds of a species on one autosome (row i of chroms.tsv).
# bcftools defaults (incl. BAQ) except -q 30 -Q 30; multiallelic caller; variant sites only.
# Per-chromosome temp BAMs with only that @SQ line: the babbler assemblies have contig names
# ('+', ',') that are illegal in VCF headers.
set -euo pipefail
source /path/to/TemporalGenomicsEvaluation/scripts/safe_mkdir.sh   # race-safe mkdir -p (uutils mkdir fails on NFS races)
B=/path/to/TemporalGenomicsEvaluation; o=$B/bcftools_roh
IFS=$'\t' read -r sp ref chr <<< "$(sed -n ${SLURM_ARRAY_TASK_ID}p $o/chroms.tsv)"
mkdir -p $o/$sp/chr || true
t=$o/$sp/tmp_$chr; rm -rf $t; mkdir -p $t; trap "rm -rf $t" EXIT
k=0
for bam in $(awk -F'\t' -v s=$sp '$1==s && $3=="modern" && $4=="raw"{print $5}' $B/hearty/manifest.tsv); do
  samtools view -h $bam $chr | awk -v c=$chr '!/^@SQ/ || $2=="SN:"c' | samtools view -b -o $t/$k.bam 2>/dev/null  # mates on other chromosomes: harmless warnings
  samtools index $t/$k.bam; echo $t/$k.bam >> $t/bams.txt; k=$((k+1))
done
bcftools mpileup -f $ref -r $chr -b $t/bams.txt -q 30 -Q 30 -a FORMAT/AD,FORMAT/DP -Ou \
  | bcftools call -mv -Ob -o $o/$sp/chr/$chr.bcf.tmp
mv $o/$sp/chr/$chr.bcf.tmp $o/$sp/chr/$chr.bcf
