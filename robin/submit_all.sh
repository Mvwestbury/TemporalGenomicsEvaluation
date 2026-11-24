#!/usr/bin/env bash
# Submit the full black robin chain (pipeline used for the Liu et al. birds, Liu-consistent). Job IDs -> robin/jobs.tsv
#   empirical : ref + reads -> map (10) -> manifests, mapDamage rescale, ANGSD native (all/tv), ATLAS (5 depth levels),
#               ANGSD depth-matched (raw, rescaled), HEARTY (raw, rescaled) -> HEARTY exact + directional singletons
#   SCAR      : BRH08 profile (after its map task) -> scar_prep (5 modern) -> scar_map (20 arms) ->
#               ANGSD (all/tv, 6x cap), HEARTY, ATLAS (5x cap) -> ATLAS-recalibrated ANGSD, directional singletons
set -euo pipefail
source /path/to/TemporalGenomicsEvaluation/scripts/safe_mkdir.sh   # race-safe mkdir -p (uutils mkdir fails on NFS races)
B=/path/to/TemporalGenomicsEvaluation
R=$B/robin; P=$B/scar_sim/profiles; S=$R/robin_scar_samples.tsv
cd $B
sb() { local name=$1; shift; local id; id=$(sbatch --parsable "$@"); printf '%s\t%s\n' "$name" "$id" >> $R/jobs.tsv; echo $id; }
: > $R/jobs.tsv

ref=$(sb ref $R/rb_ref.sbatch)
reads=$(sb reads $R/rb_reads.sbatch)
calib=$(sb calib --dependency=afterok:$ref $R/rb_atlas_calib.sbatch)
map=$(sb map --array=1-10 --dependency=afterok:$ref:$reads $R/rb_map.sbatch)

# empirical analyses
man=$(sb manifests --dependency=afterok:$map $R/rb_manifests.sbatch)
resc=$(sb rescale --array=1-10 --dependency=afterok:$map $R/rb_rescale_md.sbatch)
sb het_ind --array=1-20 --dependency=afterok:$map $R/rb_het_ind.sbatch > /dev/null
sb atlas --array=1-50%20 --dependency=afterok:$map:$calib $R/rb_atlas.sbatch > /dev/null
sb het_matched_raw --array=1-10 --dependency=afterok:$man $B/pipeline/het_matched.sbatch $B/het_matched/manifest_robin_raw.tsv > /dev/null
sb het_matched_rescaled --array=1-20 --dependency=afterok:$man:$resc $B/pipeline/het_matched.sbatch $B/het_matched/manifest_robin_rescaled.tsv > /dev/null
hr=$(sb hearty_raw --array=1-10 --dependency=afterok:$man $R/rb_hearty.sbatch)
hs=$(sb hearty_rescaled --array=11-20 --dependency=afterok:$man:$resc $R/rb_hearty.sbatch)
sb hearty_exact --array=1-20 --dependency=afterok:$hr:$hs $R/rb_hearty_exact_from_basecall.sbatch > /dev/null
sb directional --array=1-20 --dependency=afterok:$hr:$hs $R/rb_directional.sbatch > /dev/null

# SCAR simulations: BRH08 = row 2 of robin_samples.tsv -> map task 2
prof=$(sb profile --dependency=afterok:${map}_2 $R/rb_profile.sbatch)
prep=$(sb scar_prep --array=1-5 --dependency=afterok:$prof:$reads $B/pipeline/scar_prep.sbatch $S \
    $P/BRH08_ERR6335811.fraglen.txt $P/Robin_BRH08_deam.txt $P/Robin_BRH08_deam_tvex.txt)
smap=$(sb scar_map --array=1-20%10 --dependency=afterok:$prep $B/pipeline/scar_map.sbatch $S)
sb scar_het --array=1-40 --dependency=afterok:$smap $B/pipeline/scar_het.sbatch $S 6 > /dev/null
shea=$(sb scar_hearty --array=1-20 --dependency=afterok:$smap $R/rb_scar_hearty.sbatch)
satl=$(sb scar_atlas --array=1-20 --dependency=afterok:$smap:$calib $R/rb_scar_atlas.sbatch)
sb scar_recal --array=1-20 --dependency=afterok:$satl $B/pipeline/scar_recal_angsd.sbatch robin/robin_scar_samples.tsv > /dev/null
sb scar_directional --array=1-20 --dependency=afterok:$shea $R/rb_scar_directional.sbatch > /dev/null
column -t $R/jobs.tsv
