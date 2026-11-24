#!/usr/bin/env bash
#SBATCH --job-name=rb_rescale
#SBATCH --partition=all
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --output=/path/to/TemporalGenomicsEvaluation/logs/rb_rescale_%A_%a.out
# Robin copy of rescale_md.sbatch (only the sample table differs).
# mapDamage quality rescaling of the full PlainMap BAMs, as in Liu et al. (mapDamage --rescale, defaults).
# Sample = row TASK+1 of robin/robin_samples.tsv.
#  1. damage model estimated on 2M randomly sampled reads (-n 2000000, seed 1)
#  2. --rescale-only on every read, run in parallel over T groups of whole contigs (balanced by length)
#  3. merge + index -> mapped/rescaled/<sp>/<voucher>_<CRR>.rescaled.bam
set -euo pipefail
source /path/to/TemporalGenomicsEvaluation/scripts/safe_mkdir.sh   # race-safe mkdir -p (uutils mkdir fails on NFS races)
B=/path/to/TemporalGenomicsEvaluation
MD=/path/to/software/miniconda3/envs/mapdamage/bin/mapDamage
export PATH=/path/to/software/miniconda3/envs/mapdamage/bin:$PATH   # Rscript with inline/gam/RcppGSL for the Bayesian step
T=${SLURM_CPUS_PER_TASK:-8}
i=$SLURM_ARRAY_TASK_ID

IFS=$'\t' read -r species voucher crr period ref <<< "$(awk -F'\t' -v n=$((i + 1)) 'NR==n' $B/robin/robin_samples.tsv)"
sp=${species// /_}
prefix=${voucher}_${crr}
bam=$B/mapped/modern/$sp/$prefix.final.RG.bam
out=$B/mapped/rescaled/$sp/$prefix.rescaled.bam
echo "task $i: $species $voucher $crr $period host=$(hostname) $(date)"
[[ -f $bam && -f $B/mapped/modern/$sp/$prefix.plainmap.coverage.tsv ]] || { echo "SKIP: $bam not mapped yet"; exit 0; }
[[ -s $out.bai ]] && { echo "DONE already: $out"; exit 0; }

md=$B/mapped/rescaled/$sp/${prefix}_mapdamage
work=$B/mapped/rescaled/$sp/${prefix}_work
mkdir -p $md $work
rm -f $work/*

# 1. damage model on a subsample
$MD -i $bam -r $ref -d $md --merge-reference-sequences -n 2000000 --downsample-seed 1
[[ -s $md/Stats_out_MCMC_correct_prob.csv ]] || { echo "ERROR: no MCMC stats"; exit 1; }

# 2. split mapped contigs into T length-balanced groups (greedy), rescale each in parallel
samtools idxstats $bam | awk '$3>0 && $1!="*"{print $1"\t"$2}' | sort -k2,2nr |
  awk -v T=$T -v W=$work '{m=1; for(g=2;g<=T;g++) if(s[g]<s[m]) m=g; s[m]+=$2; print $1"\t0\t"$2 > W"/g"m".bed"}'
ls $work/g*.bed | xargs -P $T -I{} bash -c '
  g=$(basename {} .bed)
  samtools view -b -M -L {} -o '$work'/$g.in.bam '$bam' &&
  samtools index '$work'/$g.in.bam &&
  '$MD' --rescale-only -i '$work'/$g.in.bam -r '$ref' -d '$md' --rescale-out '$work'/$g.rescaled.bam > '$work'/$g.log 2>&1 &&
  echo "$g ok" || { echo "$g FAILED"; exit 255; }'

# 3. merge, index, sanity-check read counts
samtools merge -f -@ $T -o $out.tmp.bam $work/g*.rescaled.bam
mv $out.tmp.bam $out
samtools index -@ $T $out
nin=$(samtools idxstats $bam | awk '$1!="*"{s+=$3} END{print s}')
nout=$(samtools idxstats $out | awk '$1!="*"{s+=$3} END{print s}')
echo "mapped reads in $nin out $nout"
[[ $nin -eq $nout ]] || { echo "ERROR: read count mismatch"; exit 1; }
rm -rf $work
echo "finished $(date)"
