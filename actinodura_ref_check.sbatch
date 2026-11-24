#!/usr/bin/env bash
#SBATCH --job-name=acy_refcheck
#SBATCH --partition=all
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=/path/to/TemporalGenomicsEvaluation/logs/acy_refcheck_%j.out
# Is the Actinodura reference (GLG_GLGS6485) the same taxon as the Actinodura samples?
# Maps 250k R1 reads of (a) GLGS6485 itself and (b) an Actinodura modern control to the
# Actinodura, Trochalopteron and Yuhina references and reports mismatches per aligned bp (MAPQ>=20).
set -u   # no pipefail: head cuts the read streams short on purpose
source /path/to/TemporalGenomicsEvaluation/scripts/safe_mkdir.sh   # race-safe mkdir -p (uutils mkdir fails on NFS races)
B=/path/to/TemporalGenomicsEvaluation
W=$B/refcheck; mkdir -p $W; cd $W
N=250000
[[ -s GLGS6485.fq.gz ]] || { curl -s -r 0-99999999 --retry 20 https://download.cncb.ac.cn/gsa6/CRA020079/CRR1356155/CRR1356155_r1.fq.gz |
    (gzip -dc 2>/dev/null || true) | head -n $((N * 4)) | gzip > GLGS6485.fq.gz; }
[[ -s KIZ_23559.fq.gz ]] || { gzip -dc $B/pilot_heads/Actinodura_cyanouroptera/CRR2269019/CRR2269019_r1.fq.gz | head -n $((N * 4)) | gzip > KIZ_23559.fq.gz; }
echo "reads: GLGS6485 $(( $(gzip -dc GLGS6485.fq.gz | wc -l) / 4 ))  KIZ_23559 $(( $(gzip -dc KIZ_23559.fq.gz | wc -l) / 4 ))"
printf 'reads\treference\tmapped_q20\tmm_per_bp\tfrac_reads_lt1pct\tfrac_reads_gt6pct\n'
for s in GLGS6485 KIZ_23559; do
  for ref in GLG_GLGS6485 MT275 GLG_GLGS1218; do
    bwa mem -t 8 $B/references/full/$ref.fasta $s.fq.gz 2>/dev/null | samtools view -F 2308 -q 20 - |
      awk -v s=$s -v r=$ref '{for(i=12;i<=NF;i++) if($i~/^NM:i:/) nm=substr($i,6); al=0; c=$6;
          while(match(c,/[0-9]+[MID]/)){x=substr(c,RSTART,RLENGTH); if(x~/[MI]$/) al+=x+0; c=substr(c,RSTART+RLENGTH)}
          n++; s1+=nm; b+=al; d=nm/al; lo+=(d<0.01); hi+=(d>0.06)}
          END{printf "%s\t%s\t%d\t%.4f\t%.2f\t%.2f\n",s,r,n,s1/b,lo/n,hi/n}'
  done
done
