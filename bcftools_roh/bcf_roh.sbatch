#!/usr/bin/env bash
#SBATCH --job-name=bcf_roh
#SBATCH --partition=all
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/path/to/TemporalGenomicsEvaluation/logs/bcf_roh_%A_%a.out
# Per species (task 1-4): concat autosomes; biallelic SNPs with QUAL>=30; AF from the 5 birds' genotypes
# (+fill-tags). Then per bird: keep sites where that bird's DP is within [max(5, mean/2), 2*mean] and
# genotyped, and run bcftools roh -G30 with (1) --AF-dflt 0.4 and (2) --AF-tag AF (5-bird estimate).
set -euo pipefail
source /path/to/TemporalGenomicsEvaluation/scripts/safe_mkdir.sh   # race-safe mkdir -p (uutils mkdir fails on NFS races)
B=/path/to/TemporalGenomicsEvaluation; o=$B/bcftools_roh
sp=$(cut -f1 $o/chroms.tsv | sort -u | sed -n ${SLURM_ARRAY_TASK_ID}p); w=$o/$sp
n=$(awk -v s=$sp '$1==s' $o/chroms.tsv | wc -l); have=$(ls $w/chr/*.bcf | wc -l)
[[ $n == $have ]] || { echo "only $have of $n chromosomes"; exit 1; }
awk -v s=$sp -v w=$w '$1==s{print w"/chr/"$3".bcf"}' $o/chroms.tsv > $w/parts.txt
bcftools concat -f $w/parts.txt -Ou | bcftools view -m2 -M2 -v snps -i 'QUAL>=30' -Ou \
  | bcftools +fill-tags -Ob -o $w/snps.bcf -- -t AF,AN
bcftools index -f $w/snps.bcf
echo "SNPs: $(bcftools index -n $w/snps.bcf)"
mkdir -p $w/roh
while IFS=$'\t' read -r bam mean; do
  s=$(samtools view -H $bam | awk -F'\t' '/^@RG/{for(i=1;i<=NF;i++) if($i~/^SM:/){print substr($i,4); exit}}')
  lo=$(awk -v m=$mean 'BEGIN{x=int(m/2+0.5); print (x<5?5:x)}'); hi=$(awk -v m=$mean 'BEGIN{print int(2*m+0.5)}')
  bcftools view -s "$s" -i "FMT/DP>=$lo && FMT/DP<=$hi && GT!=\"mis\"" -Ob -o $w/roh/$s.bcf $w/snps.bcf
  echo "$s DP $lo-$hi: $(bcftools view -H $w/roh/$s.bcf | wc -l) sites, het $(bcftools view -H -g het $w/roh/$s.bcf | wc -l)"
  bcftools roh -G30 --AF-dflt 0.4 -O r -o $w/roh/$s.AFdflt.txt $w/roh/$s.bcf
  bcftools roh -G30 --AF-tag AF  -O r -o $w/roh/$s.AFest.txt  $w/roh/$s.bcf
done < <(awk -F'\t' -v s=$sp '$1==s && $3=="modern" && $4=="raw"{print $5"\t"$7}' $B/hearty/manifest.tsv)
echo "done $(date)"
