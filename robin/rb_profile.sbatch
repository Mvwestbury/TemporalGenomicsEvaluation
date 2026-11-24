#!/usr/bin/env bash
#SBATCH --job-name=rb_profile
#SBATCH --partition=all
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=/path/to/TemporalGenomicsEvaluation/logs/rb_profile_%j.out
# SCAR inputs from the historical robin BRH08 (bwa mem BAM from rb_map), built as for the Liu Fringilla arms:
#   damage  : mapDamage + TvCount on 1M reads (the scar_map.sbatch QC settings), then scripts/scar_profile.py
#             -> scar_sim/profiles/Robin_BRH08_deam.txt     (5' C->T, 3' G->A excess over the read interior)
#                scar_sim/profiles/Robin_BRH08_deam_tvex.txt  (+ G->T / C->A excess from TvCount)
#   lengths : MAPQ>=30 primary reads, 0.5% sample (samtools -s 1.005), SE read length / PE |TLEN| (first mate)
#             -> scar_sim/profiles/BRH08_ERR6335811.fraglen.txt
set -euo pipefail
source /path/to/TemporalGenomicsEvaluation/scripts/safe_mkdir.sh   # race-safe mkdir -p (uutils mkdir fails on NFS races)
B=/path/to/TemporalGenomicsEvaluation
MAPDAMAGE=/path/to/software/miniconda3/bin/mapDamage
PY=/path/to/software/miniconda3/bin/python
TVCOUNT=$B/tools/TvCount/TvCount.py
export PATH=$PATH:/path/to/software/miniconda3/envs/r_phylo_env/bin
T=${SLURM_CPUS_PER_TASK:-8}
s=BRH08_ERR6335811
ref=$B/references/full/Ptraversi_NRM_v1.fasta
bam=$B/mapped/modern/Petroica_traversi/$s.final.RG.bam
P=$B/scar_sim/profiles
q=$P/robin_$s; mkdir -p $q
[[ -s $bam.bai || -s ${bam%.bam}.bai ]] || samtools index $bam

$MAPDAMAGE -i $bam -r $ref -d $q/mapdamage --merge-reference-sequences -n 1000000 --downsample-seed 1 --no-stats
$PY $TVCOUNT --bam $bam --reference $ref --sample-label $s \
    --summary-tsv-out $q/tvcount.summary.tsv --pos-tsv-out $q/tvcount.positions.tsv \
    --plot-pdf-out $q/tvcount.pdf --plot-log-y --max-reads 1000000 --threads $T
$PY -I $B/scripts/scar_profile.py --misinc $q/mapdamage/misincorporation.txt --tvcount $q/tvcount.positions.tsv > $P/Robin_BRH08_deam_tvex.txt
$PY -I $B/scripts/scar_profile.py --misinc $q/mapdamage/misincorporation.txt > $P/Robin_BRH08_deam.txt
wc -l $P/Robin_BRH08_deam*.txt
grep -E $'^(5p\t(1|2|3)\tC\tT|3p\t(1|2|3)\tG\tA|5p\t(1|2|26)\tG\tT)\t' $P/Robin_BRH08_deam_tvex.txt

samtools view -s 1.005 -F 0xF04 -q 30 $bam | awk 'and($2,1)==0 {print length($10); next} and($2,2) && and($2,64) && $9!=0 {t=$9<0?-$9:$9; print t}' > $P/$s.fraglen.txt
sort -n $P/$s.fraglen.txt | awk '{a[NR]=$1; m+=$1} END{printf "fraglen n=%d mean=%.1f median=%d\n", NR, m/NR, a[int(NR/2)]}'
echo "done $(date)"
