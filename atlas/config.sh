#!/bin/bash
# ATLAS heterozygosity for the Liu et al. birds, adapted from Black_bears/Analysis/ATLAS_heterozygosity.
PROJECT=/path/to
B=$PROJECT/TemporalGenomicsEvaluation
ROOT=$B/atlas
ENV=$PROJECT/Software/atlas_env
ATLAS=$PROJECT/Software/atlas/build/atlas          # ATLAS 2.0.4, commit a37654b9
SCRIPTS=$ROOT/scripts
LEVELS=(native 15 10 7.5 5)                       # downsampling series (x), only where native depth allows
CALIB_SITES=10000000
export PATH="$ENV/bin:/usr/bin:/bin"
export LD_LIBRARY_PATH="$ENV/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 LC_ALL=C

# GERP source genome (Ensembl release 115, 27 sauropsids) for each target reference
gerp_source() {
    case $(basename "$1" .fasta) in
        bFriMon1) echo serinus_canaria ;;        # Fringillidae, same family as Fringilla
        *)        echo parus_major ;;            # babblers (Sylvioidea): great tit is the closest passerine
    esac
}
