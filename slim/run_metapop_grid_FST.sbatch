#!/usr/bin/env bash
# 400 independent array tasks: 20 scenarios x 20 replicates.
# At most 20 SLiM processes running concurrently (1 CPU, 1 GiB each).
# Submit from the folder containing the SLiM model:
#   mkdir -p logs tmp_out
#   sbatch --export=ALL,INPUT_DIR=/absolute/path/to/genomic/inputs run_metapop_grid_independent.sbatch

#SBATCH --job-name=metapop_grid_fst
#SBATCH --array=0-799%20
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:40:00
#SBATCH --output=logs/metapop_%A_%a.out
#SBATCH --error=logs/metapop_%A_%a.err

# time will depend on burn-in time and ancestral population size

set -euo pipefail

# Site-specific: if 'slim' is only available as a module, uncomment and adjust.
# module load SLiM

WORK_DIR="${SLURM_SUBMIT_DIR}"
SLIM_SCRIPT="${WORK_DIR}/metapop_5demes_WF.slim"
INPUT_DIR="${WORK_DIR}"
OUT_DIR="${WORK_DIR}/tmp_out_mjolnir_fst"
SLIM_BIN="PATH/slim"
mkdir -p "$OUT_DIR"

if ! command -v "$SLIM_BIN" >/dev/null 2>&1; then
  echo "ERROR: SLiM binary missing ($SLIM_BIN); load the SLiM module or set SLIM_BIN." >&2
  exit 1
fi
if [[ ! -f "$SLIM_SCRIPT" ]]; then
  echo "ERROR: Missing model: $SLIM_SCRIPT" >&2
  exit 1
fi
for required in \
    coords_SLiM_GCF_000247815.1_FicAlb1.5_chr12.txt \
    out_rho4Ne_win2e+05_col.hun.snpPair_runALL_chr12.txt; do
  if [[ ! -f "$INPUT_DIR/$required" ]]; then
    echo "ERROR: Missing genomic input: $INPUT_DIR/$required" >&2
    exit 1
  fi
done

BOTTLE_SIZES=(10 20 50 100 1000)
BOTTLE_MIGRATIONS=(0.0 0.001 0.01 0.1)
BOTTLE_LENGTHS=(10 20)
REPS_PER_SCENARIO=20

TASK_ID="${SLURM_ARRAY_TASK_ID:?Run as SLURM job array}"

N_SIZE="${#BOTTLE_SIZES[@]}"
N_MIG="${#BOTTLE_MIGRATIONS[@]}"
N_LENGTH="${#BOTTLE_LENGTHS[@]}"

N_SCENARIOS=$(( N_SIZE * N_MIG * N_LENGTH ))

if (( TASK_ID < 0 || TASK_ID >= N_SCENARIOS * REPS_PER_SCENARIO )); then
  echo "ERROR: Out-of-range array task ID $TASK_ID" >&2
  exit 1
fi

SCENARIO=$(( TASK_ID / REPS_PER_SCENARIO ))
REP=$(( TASK_ID % REPS_PER_SCENARIO + 1 ))

# Decode duration, population size, and migration rate
TIDX=$(( SCENARIO / (N_SIZE * N_MIG) ))
WITHIN_T=$(( SCENARIO % (N_SIZE * N_MIG) ))
NIDX=$(( WITHIN_T / N_MIG ))
MIDX=$(( WITHIN_T % N_MIG ))

BTL_GENS="${BOTTLE_LENGTHS[$TIDX]}"
N_BTL="${BOTTLE_SIZES[$NIDX]}"
M_BTL="${BOTTLE_MIGRATIONS[$MIDX]}"

SEED=$(( 100000 + TASK_ID + 1 ))

N_ANC="${N_ANC:-1000}"
M_ANC="${M_ANC:-0.1}"
BURN_LEN="${BURN_LEN:-50000}"  # pilot setting; lengthen if equilibrium is needed
RUN_TYPE="${RUN_TYPE:-1}"
SAMPLE_N="${SAMPLE_N:-5}"
ROH_EST="${ROH_EST:-1}"
ROH_CUTOFF="${ROH_CUTOFF:-500000}"
OUTPUT_EVERY="${OUTPUT_EVERY:-1}"

printf 'task=%s scenario=%s replicate=%s N_btl=%s M_btl=%s seed=%s start=%s\n' \
  "$TASK_ID" "$SCENARIO" "$REP" "$N_BTL" "$M_BTL" "$SEED" "$(date -Is)"

# Skip a previously completed replicate when resubmitting an array subset.
ANC_TAG="${N_ANC}-${N_ANC}-${N_ANC}-${N_ANC}-${N_ANC}"
BTL_TAG="${N_BTL}-${N_BTL}-${N_BTL}-${N_BTL}-${N_BTL}"
MA_TAG="${M_ANC//./p}"
MB_TAG="${M_BTL//./p}"
PREFIX="${OUT_DIR}/metapop_Na${ANC_TAG}_Nb${BTL_TAG}_Ma${MA_TAG}_Mb${MB_TAG}_T${BTL_GENS}_B${BURN_LEN}_R${RUN_TYPE}_S${SEED}"
if [[ -s "${PREFIX}_comparison.tsv" && -s "${PREFIX}_fst.tsv" ]]; then
  echo "SKIP completed task=$TASK_ID seed=$SEED"
  exit 0
fi

"$SLIM_BIN" -s "$SEED" \
  -d "basePath=\"${INPUT_DIR}\"" \
  -d "outDir=\"${OUT_DIR}\"" \
  -d "burnLen=${BURN_LEN}" \
  -d "btl_gens=${BTL_GENS}" \
  -d "N_anc1=${N_ANC}" \
  -d "N_anc2=${N_ANC}" \
  -d "N_anc3=${N_ANC}" \
  -d "N_anc4=${N_ANC}" \
  -d "N_anc5=${N_ANC}" \
  -d "N_btl1=${N_BTL}" \
  -d "N_btl2=${N_BTL}" \
  -d "N_btl3=${N_BTL}" \
  -d "N_btl4=${N_BTL}" \
  -d "N_btl5=${N_BTL}" \
  -d "M_anc=${M_ANC}" \
  -d "M_btl=${M_BTL}" \
  -d "runType=${RUN_TYPE}" \
  -d "sampleN=${SAMPLE_N}" \
  -d "ROHest=${ROH_EST}" \
  -d "rohCutoff=${ROH_CUTOFF}" \
  -d "outputEvery=${OUTPUT_EVERY}" \
  "$SLIM_SCRIPT"

printf 'DONE task=%s scenario=%s replicate=%s at %s\n' "$TASK_ID" "$SCENARIO" "$REP" "$(date -Is)"
