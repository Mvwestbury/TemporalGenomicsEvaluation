#!/usr/bin/env bash
#SBATCH --job-name=scar_prep
#SBATCH --partition=all
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=/path/to/TemporalGenomicsEvaluation/logs/scar_prep_%A_%a.out
# Synthetic "historical" versions of modern samples with SCAR (degradation control).
# usage: sbatch --array=1-N scar_prep.sbatch <samples.tsv> <fraglen.txt> <profile_deam.txt> <profile_deam_tvex.txt>
#   samples.tsv rows: species  voucher  CRR  ref   (no header)
# 1. fastp: adapter trimming + merging of overlapping pairs (-l 30)
# 2. SE input = merged + unmerged/unpaired R1 + reverse-complemented unmerged/unpaired R2
#    (all reads in the orientation of the original molecule's top strand)
# 3. SCAR arms: frag (empirical historical lengths), frag+deam, frag+deam+tvex
#    (damage applied to the same fragmented reads, so the arms differ only by damage)
set -euo pipefail
source /path/to/TemporalGenomicsEvaluation/scripts/safe_mkdir.sh   # race-safe mkdir -p (uutils mkdir fails on NFS races)
B=/path/to/TemporalGenomicsEvaluation
SCAR=/path/to/software/SCAR-private/bin/scar
FASTP=/path/to/software/fastp
T=${SLURM_CPUS_PER_TASK:-16}
SAMPLES=$1; FRAG=$2; PDEAM=$3; PDEAMTV=$4

IFS=$'\t' read -r species voucher crr ref <<< "$(sed -n ${SLURM_ARRAY_TASK_ID}p $SAMPLES)"
sp=${species// /_}
prefix=${voucher}_${crr}
reads=$B/raw_reads/$sp/$crr
d=$B/scar_sim/$sp/$prefix
mkdir -p $d
echo "task $SLURM_ARRAY_TASK_ID: $prefix host=$(hostname) $(date)"
for f in $reads/${crr}_r1.fq.gz $reads/${crr}_r2.fq.gz; do
    [[ -f $f.md5ok ]] || { echo "SKIP: $f not verified"; exit 0; }
done
cp $FRAG $d/fraglen.txt; cp $PDEAM $d/profile_deam.txt; cp $PDEAMTV $d/profile_deam_tvex.txt

# 1. merge + trim (skipped if a complete original.fq.gz exists from an earlier run)
if [[ ! -s $d/original.fq.gz || -e $d/m.fq.gz ]]; then
$FASTP -i $reads/${crr}_r1.fq.gz -I $reads/${crr}_r2.fq.gz \
    --merge --merged_out $d/m.fq.gz -o $d/u1.fq.gz -O $d/u2.fq.gz \
    --unpaired1 $d/s1.fq.gz --unpaired2 $d/s2.fq.gz \
    --detect_adapter_for_pe -l 30 -w $T -j $d/fastp.json -h $d/fastp.html 2> $d/fastp.log

# 2. one SE file; R2-derived reads reverse-complemented and renamed (_r2) to keep names unique
#    (fastp leaves a 0-byte file, not an empty gzip, when an output gets no reads -> skip those)
nz() { for f; do [[ -s $f ]] && zcat $f; done; return 0; }
{ nz $d/m.fq.gz $d/u1.fq.gz $d/s1.fq.gz
  seqtk seq -r <(nz $d/u2.fq.gz $d/s2.fq.gz) | awk 'NR%4==1{sub(/[ \t].*/, ""); $0=$0"_r2"} {print}'
} | bgzip -@ $T > $d/original.fq.gz
rm -f $d/m.fq.gz $d/u1.fq.gz $d/u2.fq.gz $d/s1.fq.gz $d/s2.fq.gz
fi

# 3. SCAR
cd $d
$SCAR -i original.fq.gz -o frag --fd empirical --fdf fraglen.txt --mfl 30 --mutation-rate 0 --seed 11 --threads $T > scar_frag.log 2>&1
$SCAR -i frag.fastq.gz -o frag_deam    -d profile_deam.txt    --seed 12 --threads $T > scar_deam.log 2>&1
$SCAR -i frag.fastq.gz -o frag_deam_tvex -d profile_deam_tvex.txt --seed 12 --threads $T > scar_deam_tvex.log 2>&1

for f in original.fq.gz frag.fastq.gz frag_deam.fastq.gz frag_deam_tvex.fastq.gz; do
    zcat $f | awk -v f=$f 'NR%4==2{n++; b+=length($0)} END{printf "%s reads=%d bases=%d mean_len=%.1f\n", f, n, b, b/n}'
done
echo "done $(date)"
