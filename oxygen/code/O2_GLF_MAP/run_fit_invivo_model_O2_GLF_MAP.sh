#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FIT_SCRIPT="${SCRIPT_DIR}/fit_invivo_model_O2_GLF_MAP.R"

RUN_PREFIX="fit_invivo_model_O2_GLF_MAP_$(date +%Y%m%d_%H%M%S)"
OUT_ROOT="${SCRIPT_DIR}/../../results"
DATA_DIR="${SCRIPT_DIR}/../../../data/InVivoData_Gemcitabine"
SEEDS_FILE="${SCRIPT_DIR}/seeds.csv"
SEEDS_CSV=""
N_CORES=8
USE_DEOPTIM=TRUE
DEOPTIM_PARALLEL=TRUE
FIT_TREATMENT=FALSE
DOSE_ZERO_ONLY=TRUE
PAIRED_ONLY=TRUE
TRUNCATE_AT_TREATMENT=FALSE
PLOIDY_AT_HARVEST=TRUE
ITERMAX=220
NP=256
N_STARTS=80
OPTIM_MAXIT=15000
SIGMA_BURDEN=0.35
SIGMA_PLOIDY=0.15
USE_SOFT_PRIOR=TRUE
LAMBDA_PRIOR=0.03
TAU_O2=3
O2_CURVE_TYPE=glogistic
O2_CAP_PCT=5
O2_ANCHOR_N=1e6
PARAMETER_TABLE="${SCRIPT_DIR}/parameter_table.csv"

EXTRA_ARGS=()
for arg in "$@"; do
  case "$arg" in
    --run_prefix=*) RUN_PREFIX="${arg#*=}" ;;
    --out_root=*) OUT_ROOT="${arg#*=}" ;;
    --data_dir=*) DATA_DIR="${arg#*=}" ;;
    --seeds_file=*) SEEDS_FILE="${arg#*=}" ;;
    --seeds_csv=*) SEEDS_CSV="${arg#*=}" ;;
    --n_cores=*) N_CORES="${arg#*=}" ;;
    --use_deoptim=*) USE_DEOPTIM="${arg#*=}" ;;
    --deoptim_parallel=*) DEOPTIM_PARALLEL="${arg#*=}" ;;
    --fit_treatment=*) FIT_TREATMENT="${arg#*=}" ;;
    --dose_zero_only=*) DOSE_ZERO_ONLY="${arg#*=}" ;;
    --paired_only=*) PAIRED_ONLY="${arg#*=}" ;;
    --truncate_at_treatment=*) TRUNCATE_AT_TREATMENT="${arg#*=}" ;;
    --ploidy_at_harvest=*) PLOIDY_AT_HARVEST="${arg#*=}" ;;
    --itermax=*) ITERMAX="${arg#*=}" ;;
    --np=*|--NP=*) NP="${arg#*=}" ;;
    --n_starts=*) N_STARTS="${arg#*=}" ;;
    --optim_maxit=*) OPTIM_MAXIT="${arg#*=}" ;;
    --sigma_burden=*) SIGMA_BURDEN="${arg#*=}" ;;
    --sigma_ploidy=*) SIGMA_PLOIDY="${arg#*=}" ;;
    --use_soft_prior=*) USE_SOFT_PRIOR="${arg#*=}" ;;
    --lambda_prior=*) LAMBDA_PRIOR="${arg#*=}" ;;
    --tau_O2=*) TAU_O2="${arg#*=}" ;;
    --o2_curve_type=*) O2_CURVE_TYPE="${arg#*=}" ;;
    --o2_cap_pct=*) O2_CAP_PCT="${arg#*=}" ;;
    --o2_anchor_N=*) O2_ANCHOR_N="${arg#*=}" ;;
    --parameter_table=*) PARAMETER_TABLE="${arg#*=}" ;;
    *) EXTRA_ARGS+=("$arg") ;;
  esac
done

read_seeds_from_file() {
  local f="$1"
  Rscript --vanilla - "$f" <<'RS'
args <- commandArgs(trailingOnly = TRUE)
f <- args[[1]]
if (!file.exists(f)) { cat(""); quit(save = "no", status = 0) }
ln <- readLines(f, warn = FALSE)
ln <- gsub("\\r", "", ln)
ln <- ln[nzchar(trimws(ln))]
if (length(ln) == 0L) { cat(""); quit(save = "no", status = 0) }
raw <- paste(ln, collapse = ",")
parts <- trimws(unlist(strsplit(raw, "[,]", fixed = FALSE)))
parts <- parts[nzchar(parts)]
parts <- parts[!tolower(parts) %in% c("seed", "seeds")]
nums <- suppressWarnings(as.integer(parts))
nums <- nums[is.finite(nums)]
if (length(nums) == 0L) { cat(""); quit(save = "no", status = 0) }
cat(paste(nums, collapse = ","))
RS
}

SEEDS_FROM_FILE=""
if [[ -n "${SEEDS_FILE}" ]]; then
  SEEDS_FROM_FILE="$(read_seeds_from_file "${SEEDS_FILE}")"
fi
if [[ -n "${SEEDS_FROM_FILE}" ]]; then
  SEEDS_USE="${SEEDS_FROM_FILE}"
  SEED_SOURCE="file:${SEEDS_FILE}"
else
  SEEDS_USE="${SEEDS_CSV}"
  SEED_SOURCE="arg:--seeds_csv"
fi
if [[ -z "${SEEDS_USE}" ]]; then
  echo "ERROR: no seeds found. Provide ${SEEDS_FILE} or --seeds_csv." >&2
  exit 1
fi

mkdir -p "${OUT_ROOT}"

echo "Running O2_GLF_MAP"
echo "  Fit script: ${FIT_SCRIPT}"
echo "  Seeds: ${SEEDS_USE} (${SEED_SOURCE})"
echo "  Parameter table: ${PARAMETER_TABLE}"

IFS=',' read -r -a seed_arr <<< "${SEEDS_USE}"
for seed in "${seed_arr[@]}"; do
  seed="$(echo "$seed" | tr -d '[:space:]')"
  [[ -z "$seed" ]] && continue
  run_dir="${OUT_ROOT}/${RUN_PREFIX}_seed${seed}"
  mkdir -p "${run_dir}"
  cmd=(
    Rscript "${FIT_SCRIPT}"
    "--seed=${seed}"
    "--out_dir=${run_dir}"
    "--data_dir=${DATA_DIR}"
    "--n_cores=${N_CORES}"
    "--use_deoptim=${USE_DEOPTIM}"
    "--deoptim_parallel=${DEOPTIM_PARALLEL}"
    "--fit_treatment=${FIT_TREATMENT}"
    "--dose_zero_only=${DOSE_ZERO_ONLY}"
    "--paired_only=${PAIRED_ONLY}"
    "--truncate_at_treatment=${TRUNCATE_AT_TREATMENT}"
    "--ploidy_at_harvest=${PLOIDY_AT_HARVEST}"
    "--itermax=${ITERMAX}"
    "--NP=${NP}"
    "--n_starts=${N_STARTS}"
    "--optim_maxit=${OPTIM_MAXIT}"
    "--sigma_burden=${SIGMA_BURDEN}"
    "--sigma_ploidy=${SIGMA_PLOIDY}"
    "--use_soft_prior=${USE_SOFT_PRIOR}"
    "--lambda_prior=${LAMBDA_PRIOR}"
    "--tau_O2=${TAU_O2}"
    "--o2_curve_type=${O2_CURVE_TYPE}"
    "--o2_cap_pct=${O2_CAP_PCT}"
    "--o2_anchor_N=${O2_ANCHOR_N}"
    "--parameter_table=${PARAMETER_TABLE}"
  )
  if [[ ${#EXTRA_ARGS[@]} -gt 0 ]]; then
    cmd+=("${EXTRA_ARGS[@]}")
  fi
  echo "[$(date '+%F %T')] seed=${seed}: start"
  echo "Command: ${cmd[*]}"
  "${cmd[@]}"
  echo "[$(date '+%F %T')] seed=${seed}: done"
done

echo "All done. Output root: ${OUT_ROOT}"
