#!/usr/bin/env bash
#SBATCH --job-name=acy_newref
#SBATCH --partition=all
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/path/to/TemporalGenomicsEvaluation/logs/acy_newref_%j.out
# Candidate Actinodura references: download the A. strigula (GLGS5100) and A. egertoni (GLGS5798)
# assemblies, index them, and report mismatches per aligned bp (MAPQ>=20) for A. cyanouroptera reads.
set -u
source /path/to/TemporalGenomicsEvaluation/scripts/safe_mkdir.sh   # race-safe mkdir -p (uutils mkdir fails on NFS races)
B=/path/to/TemporalGenomicsEvaluation
R=$B/references/full
cd $R
for id in GLG_GLGS5100 GLG_GLGS5798; do
  if [[ ! -f $id.fasta ]]; then
    url=$(grep "10.57760sciencedb.17191/$id.zip" $B/metadata/sciencedb17191_links.txt)
    curl -s --retry 20 -o $id.zip "$url"
    /path/to/software/miniconda3/bin/python -I -c "import zipfile,sys; z=zipfile.ZipFile(sys.argv[1]); n=[x for x in z.namelist() if x.endswith(('.fa','.fasta'))][0]; open(sys.argv[2],'wb').write(z.read(n))" $id.zip $id.fasta && rm $id.zip
  fi
  [[ -f $id.fasta.fai ]] || samtools faidx $id.fasta
  [[ -f $id.fasta.bwt ]] || bwa index $id.fasta 2> $id.fasta.bwa_index.log &
done
wait
cd $B/refcheck
gzip -dc $B/pilot_heads/Actinodura_cyanouroptera/CRR2269020/CRR2269020_r1.fq.gz | head -n 1000000 | gzip > KIZ_003002.fq.gz
printf 'reads\treference\tmapped_q20\tmm_per_bp\tfrac_reads_lt1pct\tfrac_reads_gt6pct\n'
for s in KIZ_23559 KIZ_003002 GLGS5100 GLGS5798; do
  for ref in GLG_GLGS5100 GLG_GLGS5798; do
    bwa mem -t 8 $R/$ref.fasta $s.fq.gz 2>/dev/null | samtools view -F 2308 -q 20 - |
      awk -v s=$s -v r=$ref '{for(i=12;i<=NF;i++) if($i~/^NM:i:/) nm=substr($i,6); al=0; c=$6;
          while(match(c,/[0-9]+[MID]/)){x=substr(c,RSTART,RLENGTH); if(x~/[MI]$/) al+=x+0; c=substr(c,RSTART+RLENGTH)}
          n++; s1+=nm; b+=al; d=nm/al; lo+=(d<0.01); hi+=(d>0.06)}
          END{if(n) printf "%s\t%s\t%d\t%.4f\t%.2f\t%.2f\n",s,r,n,s1/b,lo/n,hi/n}'
  done
done
/path/to/software/miniconda3/bin/python $B/scripts/assembly_stats.py $R/GLG_GLGS5100.fasta $R/GLG_GLGS5798.fasta
