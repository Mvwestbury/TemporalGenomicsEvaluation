#!/bin/bash
#SBATCH --job-name=rb_calib
#SBATCH --partition=all
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --array=1
#SBATCH --output=/path/to/TemporalGenomicsEvaluation/logs/rb_calib_%A_%a.out
# Robin copy of atlas/calib.sbatch (GERP source: parus_major, the config.sh default). Per target reference: lift the GERP top-20M source sites (ATLAS liftOver, 100-bp flanks, bwa mem, MAPQ>=30),
# keep the top 10M on autosomes as calibration sites, and write fixed 1-Mb evaluation windows
# (one centred window per 10-Mb autosome block), as in the black bear pipeline.
set -euo pipefail
source /path/to/TemporalGenomicsEvaluation/scripts/safe_mkdir.sh   # race-safe mkdir -p (uutils mkdir fails on NFS races)
source /path/to/TemporalGenomicsEvaluation/atlas/config.sh
refs=(Ptraversi_NRM_v1)
r=${refs[$((SLURM_ARRAY_TASK_ID - 1))]}
REF=$B/references/full/$r.fasta
src=$(gerp_source $REF)
S=$ROOT/resources/$src
[[ -e $S/GERP_RANKING_SUCCESS ]]
d=$ROOT/calib/$r
mkdir -p $d; cd $d
sed 's/:$//' $B/references/regions/$r.autosomes.txt > autosomes.txt

fa=$(ls $S/*.dna.toplevel.fa)
[[ -s source_top20M.fastq ]] || "$ATLAS" liftOver --mode Bed2Fastq --bed $S/gerp_source_top20M.bed \
    --fasta $fa --flank 100 --out source_top20M
bwa mem -M -t 16 $REF source_top20M.fastq |
    samtools view -u -q 30 -F 3844 - |
    samtools sort -@ 4 -m 2G -o liftover.bam -
samtools index liftover.bam
"$ATLAS" liftOver --mode Bam2Bed --bam liftover.bam --filterMQ 30, --out lifted
python $SCRIPTS/conserved_sites.py lifted --input lifted.bed --output calibration.bed \
    --count $CALIB_SITES --autosomes autosomes.txt
cat calibration.json

python - $REF.fai autosomes.txt evaluation.bed <<'EOF'
import sys
fai, auto, out = sys.argv[1:]
allowed = set(open(auto).read().split())
n = 0
with open(out, "w") as o:
    for line in open(fai):
        chrom, length = line.split()[:2]
        length = int(length)
        if chrom not in allowed:
            continue
        for left in range(0, length, 10_000_000):
            right = min(left + 10_000_000, length)
            if right - left < 2_000_000:
                continue
            start = (left + right - 1_000_000) // 2
            o.write(f"{chrom}\t{start}\t{start + 1_000_000}\n"); n += 1
print(f"{n} evaluation windows")
EOF
rm -f source_top20M.fastq
touch CALIBRATION_SUCCESS
