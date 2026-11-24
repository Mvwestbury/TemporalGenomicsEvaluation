#!/usr/bin/env bash
#SBATCH --job-name=mm_diag
#SBATCH --partition=all
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=/path/to/TemporalGenomicsEvaluation/logs/mm_diag_%j.out
source /path/to/TemporalGenomicsEvaluation/scripts/safe_mkdir.sh   # race-safe mkdir -p (uutils mkdir fails on NFS races)
# Mismatch diagnostics for the modern controls of all four species plus two historical Fringilla.
B=/path/to/TemporalGenomicsEvaluation
P=$B/pilot/modern
mkdir -p $B/pilot/mismatch_diag
/path/to/software/miniconda3/bin/python $B/scripts/mismatch_diagnostics.py $B/pilot/mismatch_diag/diag \
    $P/Fringilla_montifringilla/KIZ_GLGS5757_CRR2269155.final.RG.bam $P/Fringilla_montifringilla/KIZ_GLGS5758_CRR2269156.final.RG.bam \
    $P/Fringilla_montifringilla/KIZ_013083_CRR2269150.final.RG.bam $P/Fringilla_montifringilla/KIZ_013086_CRR2269149.final.RG.bam \
    $P/Yuhina_occipitalis/KIZ_23021_CRR2269428.final.RG.bam $P/Yuhina_occipitalis/KIZ_23059_CRR2269425.final.RG.bam \
    $P/Trochalopteron_erythrocephalum/KIZ_GLGS5059_CRR2269388.final.RG.bam $P/Trochalopteron_erythrocephalum/KIZ_GLGS5728_CRR1356041.final.RG.bam
