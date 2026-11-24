#!/bin/bash
#SBATCH --job-name=atlas_gerp
#SBATCH --partition=all
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --array=1-2
#SBATCH --output=/path/to/TemporalGenomicsEvaluation/logs/atlas_gerp_%A_%a.out
# Download Ensembl 115 sauropsid GERP scores + genome for one source species and rank the top 20M sites.
set -euo pipefail
source /path/to/TemporalGenomicsEvaluation/scripts/safe_mkdir.sh   # race-safe mkdir -p (uutils mkdir fails on NFS races)
source /path/to/TemporalGenomicsEvaluation/atlas/config.sh
case $SLURM_ARRAY_TASK_ID in
    1) sp=serinus_canaria; asm=SCA1; Sp=Serinus_canaria ;;
    2) sp=parus_major; asm=Parus_major1.1; Sp=Parus_major ;;
esac
d=$ROOT/resources/$sp
mkdir -p $d; cd $d
gerp_url=https://ftp.ensembl.org/pub/release-115/compara/conservation_scores/27_sauropsids.gerp_conservation_score
fasta_url=https://ftp.ensembl.org/pub/release-115/fasta/$sp/dna
fetch() { [[ -s $2 ]] || { curl -fL --retry 5 --connect-timeout 60 -C - "$1" -o $2.part; mv $2.part $2; }; }
bw=gerp_conservation_scores.$sp.$asm.bw
fa=$Sp.$asm.dna.toplevel.fa
fetch $gerp_url/MD5SUM gerp.MD5SUM
fetch $gerp_url/$bw $bw
grep " $bw\$\|  $bw\$" gerp.MD5SUM | md5sum -c -
fetch $fasta_url/$fa.gz $fa.gz
gzip -t $fa.gz
[[ -s $fa ]] || { gzip -cd $fa.gz > $fa.tmp; mv $fa.tmp $fa; }
samtools faidx $fa
python $SCRIPTS/conserved_sites.py source --input $bw --output gerp_source_top20M.bed --count 20000000
cat gerp_source_top20M.json
touch GERP_RANKING_SUCCESS
