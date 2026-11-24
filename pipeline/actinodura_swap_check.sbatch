#!/usr/bin/env bash
#SBATCH --job-name=acy_swap
#SBATCH --partition=all
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=/path/to/TemporalGenomicsEvaluation/logs/acy_swap_%j.out
# Were the Actinodura assemblies swapped? Maps 250k R1 reads of the A. strigula (GLGS5100) and
# A. egertoni (GLGS5798) assembly individuals to GLG_GLGS6485 and MT275; mismatches per aligned bp (MAPQ>=20).
set -u   # no pipefail: head cuts the read streams short on purpose
source /path/to/TemporalGenomicsEvaluation/scripts/safe_mkdir.sh   # race-safe mkdir -p (uutils mkdir fails on NFS races)
B=/path/to/TemporalGenomicsEvaluation
W=$B/refcheck; mkdir -p $W; cd $W
N=250000
for x in GLGS5100:CRR1355784 GLGS5798:CRR1356269; do
  id=${x%:*}; crr=${x#*:}
  [[ -s $id.fq.gz ]] || curl -s -r 0-99999999 --retry 20 https://download.cncb.ac.cn/gsa6/CRA020079/$crr/${crr}_r1.fq.gz |
      (gzip -dc 2>/dev/null || true) | head -n $((N * 4)) | gzip > $id.fq.gz
done
printf 'reads\treference\tmapped_q20\tmm_per_bp\tfrac_reads_lt1pct\tfrac_reads_gt6pct\n'
for s in GLGS5100 GLGS5798; do
  for ref in GLG_GLGS6485 MT275; do
    bwa mem -t 8 $B/references/full/$ref.fasta $s.fq.gz 2>/dev/null | samtools view -F 2308 -q 20 - |
      awk -v s=$s -v r=$ref '{for(i=12;i<=NF;i++) if($i~/^NM:i:/) nm=substr($i,6); al=0; c=$6;
          while(match(c,/[0-9]+[MID]/)){x=substr(c,RSTART,RLENGTH); if(x~/[MI]$/) al+=x+0; c=substr(c,RSTART+RLENGTH)}
          n++; s1+=nm; b+=al; d=nm/al; lo+=(d<0.01); hi+=(d>0.06)}
          END{printf "%s\t%s\t%d\t%.4f\t%.2f\t%.2f\n",s,r,n,s1/b,lo/n,hi/n}'
  done
done
