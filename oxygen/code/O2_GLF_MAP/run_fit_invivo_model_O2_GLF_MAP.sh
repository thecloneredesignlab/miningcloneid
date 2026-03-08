#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OXYGEN_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
FIT_SCRIPT="${SCRIPT_DIR}/fit_invivo_model_O2_GLF_MAP.R"
VIZ_SCRIPT="${SCRIPT_DIR}/viz_invivo_model_O2_GLF_MAP_results.R"

DEFAULT_CONFIG="${OXYGEN_DIR}/config/O2_GLF.yaml"
CONFIG_FILE="${DEFAULT_CONFIG}"

RUN_PREFIX="fit_invivo_model_O2_GLF_MAP_$(date +%Y%m%d_%H%M%S)"
OUT_ROOT="${OXYGEN_DIR}/results"
DATA_DIR="${SCRIPT_DIR}/../../../data/InVivoData_Gemcitabine"
SEEDS_FILE="${OXYGEN_DIR}/data/O2_GLF/seeds.csv"
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
SIGMA_PLOIDY=0.08
BURDEN_EXCLUDE_DAY0=TRUE
USE_SOFT_PRIOR=TRUE
LAMBDA_PRIOR=0.03
TAU_O2=3
O2_CURVE_TYPE=glogistic
O2_CAP_PCT=5
O2_ANCHOR_N=1e6
PARAMETER_TABLE="${OXYGEN_DIR}/data/O2_GLF/parameter_table.csv"
APPEND_RUN_PREFIX_TIMESTAMP=FALSE
RUN_PREFIX_TIMESTAMP_FORMAT="%Y%m%d_%H%M%S"
AUTO_VIZ=TRUE
VIZ_REPORT_DT=1
VIZ_TOP_N=6

trim() {
  local s="$1"
  s="${s#"${s%%[![:space:]]*}"}"
  s="${s%"${s##*[![:space:]]}"}"
  printf '%s' "$s"
}

resolve_path() {
  local p="$1"
  local base="$2"
  if [[ -z "$p" ]]; then
    printf '%s' "$p"
    return
  fi
  if [[ "$p" == /* ]]; then
    printf '%s' "$p"
    return
  fi
  if [[ "$p" == ~* ]]; then
    eval printf '%s' "$p"
    return
  fi
  printf '%s' "${base}/$p"
}

load_config() {
  local cfg_path="$1"
  local cfg_dir="$2"
  [[ -f "$cfg_path" ]] || { echo "ERROR: config not found: $cfg_path" >&2; exit 1; }

  while IFS= read -r raw || [[ -n "$raw" ]]; do
    local line="${raw%%#*}"
    line="$(trim "$line")"
    [[ -z "$line" ]] && continue
    [[ "$line" != *:* ]] && continue

    local key="${line%%:*}"
    local val="${line#*:}"
    key="$(trim "$key")"
    val="$(trim "$val")"

    if [[ "$val" == "\""*"\"" && "$val" == *"\"" ]]; then
      val="${val:1:${#val}-2}"
    elif [[ "$val" == "'"*"'" && "$val" == *"'" ]]; then
      val="${val:1:${#val}-2}"
    fi

    case "$key" in
      run_prefix) RUN_PREFIX="$val" ;;
      out_root) OUT_ROOT="$(resolve_path "$val" "$cfg_dir")" ;;
      data_dir) DATA_DIR="$(resolve_path "$val" "$cfg_dir")" ;;
      seeds_file) SEEDS_FILE="$(resolve_path "$val" "$cfg_dir")" ;;
      seeds_csv) SEEDS_CSV="$val" ;;
      n_cores) N_CORES="$val" ;;
      use_deoptim) USE_DEOPTIM="$val" ;;
      deoptim_parallel) DEOPTIM_PARALLEL="$val" ;;
      fit_treatment) FIT_TREATMENT="$val" ;;
      dose_zero_only) DOSE_ZERO_ONLY="$val" ;;
      paired_only) PAIRED_ONLY="$val" ;;
      truncate_at_treatment) TRUNCATE_AT_TREATMENT="$val" ;;
      ploidy_at_harvest) PLOIDY_AT_HARVEST="$val" ;;
      itermax) ITERMAX="$val" ;;
      np|NP) NP="$val" ;;
      n_starts) N_STARTS="$val" ;;
      optim_maxit) OPTIM_MAXIT="$val" ;;
      sigma_burden) SIGMA_BURDEN="$val" ;;
      sigma_ploidy) SIGMA_PLOIDY="$val" ;;
      burden_exclude_day0) BURDEN_EXCLUDE_DAY0="$val" ;;
      use_soft_prior) USE_SOFT_PRIOR="$val" ;;
      lambda_prior) LAMBDA_PRIOR="$val" ;;
      tau_O2) TAU_O2="$val" ;;
      o2_curve_type) O2_CURVE_TYPE="$val" ;;
      o2_cap_pct) O2_CAP_PCT="$val" ;;
      o2_anchor_N) O2_ANCHOR_N="$val" ;;
      parameter_table) PARAMETER_TABLE="$(resolve_path "$val" "$cfg_dir")" ;;
      append_run_prefix_timestamp) APPEND_RUN_PREFIX_TIMESTAMP="$val" ;;
      run_prefix_timestamp_format) RUN_PREFIX_TIMESTAMP_FORMAT="$val" ;;
      auto_viz) AUTO_VIZ="$val" ;;
      viz_report_dt) VIZ_REPORT_DT="$val" ;;
      viz_top_n) VIZ_TOP_N="$val" ;;
    esac
  done < "$cfg_path"
}

for arg in "$@"; do
  case "$arg" in
    --config=*) CONFIG_FILE="${arg#*=}" ;;
  esac
done

if [[ "$CONFIG_FILE" != /* ]]; then
  CONFIG_FILE="$(resolve_path "$CONFIG_FILE" "$PWD")"
fi
CONFIG_DIR="$(cd "$(dirname "$CONFIG_FILE")" && pwd)"
load_config "$CONFIG_FILE" "$CONFIG_DIR"

EXTRA_ARGS=()
for arg in "$@"; do
  case "$arg" in
    --config=*) ;;
    --run_prefix=*) RUN_PREFIX="${arg#*=}" ;;
    --out_root=*) OUT_ROOT="$(resolve_path "${arg#*=}" "$PWD")" ;;
    --data_dir=*) DATA_DIR="$(resolve_path "${arg#*=}" "$PWD")" ;;
    --seeds_file=*) SEEDS_FILE="$(resolve_path "${arg#*=}" "$PWD")" ;;
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
    --burden_exclude_day0=*) BURDEN_EXCLUDE_DAY0="${arg#*=}" ;;
    --use_soft_prior=*) USE_SOFT_PRIOR="${arg#*=}" ;;
    --lambda_prior=*) LAMBDA_PRIOR="${arg#*=}" ;;
    --tau_O2=*) TAU_O2="${arg#*=}" ;;
    --o2_curve_type=*) O2_CURVE_TYPE="${arg#*=}" ;;
    --o2_cap_pct=*) O2_CAP_PCT="${arg#*=}" ;;
    --o2_anchor_N=*) O2_ANCHOR_N="${arg#*=}" ;;
    --parameter_table=*) PARAMETER_TABLE="$(resolve_path "${arg#*=}" "$PWD")" ;;
    --append_run_prefix_timestamp=*) APPEND_RUN_PREFIX_TIMESTAMP="${arg#*=}" ;;
    --run_prefix_timestamp_format=*) RUN_PREFIX_TIMESTAMP_FORMAT="${arg#*=}" ;;
    --auto_viz=*) AUTO_VIZ="${arg#*=}" ;;
    --viz_report_dt=*) VIZ_REPORT_DT="${arg#*=}" ;;
    --viz_top_n=*) VIZ_TOP_N="${arg#*=}" ;;
    *) EXTRA_ARGS+=("$arg") ;;
  esac
done

if [[ "${APPEND_RUN_PREFIX_TIMESTAMP}" == "TRUE" || "${APPEND_RUN_PREFIX_TIMESTAMP}" == "true" || "${APPEND_RUN_PREFIX_TIMESTAMP}" == "1" ]]; then
  ts_suffix="$(date +"${RUN_PREFIX_TIMESTAMP_FORMAT}")"
  RUN_PREFIX="${RUN_PREFIX}_${ts_suffix}"
fi

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
RUN_DIR="${OUT_ROOT}/${RUN_PREFIX}"
mkdir -p "${RUN_DIR}"
RUN_LOG="${RUN_DIR}/run_status.log"
touch "${RUN_LOG}"
exec > >(tee -a "${RUN_LOG}") 2>&1

echo "Running O2_GLF_MAP"
echo "  Config: ${CONFIG_FILE}"
echo "  Fit script: ${FIT_SCRIPT}"
echo "  Data dir: ${DATA_DIR}"
echo "  Seeds: ${SEEDS_USE} (${SEED_SOURCE})"
echo "  Parameter table: ${PARAMETER_TABLE}"
echo "  Run dir: ${RUN_DIR}"
echo "  Run log: ${RUN_LOG}"
echo "  Run prefix timestamp suffix: ${APPEND_RUN_PREFIX_TIMESTAMP} (format=${RUN_PREFIX_TIMESTAMP_FORMAT})"
echo "  Auto viz: ${AUTO_VIZ} (report_dt=${VIZ_REPORT_DT}, top_n=${VIZ_TOP_N})"
echo "  burden_exclude_day0: ${BURDEN_EXCLUDE_DAY0}"

IFS=',' read -r -a seed_arr <<< "${SEEDS_USE}"
for seed in "${seed_arr[@]}"; do
  seed="$(echo "$seed" | tr -d '[:space:]')"
  [[ -z "$seed" ]] && continue
  run_dir="${RUN_DIR}/seed${seed}"
  mkdir -p "${run_dir}"
  fit_log="${run_dir}/fit_status.log"
  viz_log="${run_dir}/viz_status.log"
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
    "--burden_exclude_day0=${BURDEN_EXCLUDE_DAY0}"
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
  echo "[$(date '+%F %T')] seed=${seed}: fit_log=${fit_log}"
  echo "Command: ${cmd[*]}"
  "${cmd[@]}" 2>&1 | tee "${fit_log}"
  echo "[$(date '+%F %T')] seed=${seed}: done"

  if [[ "${AUTO_VIZ}" == "TRUE" || "${AUTO_VIZ}" == "true" || "${AUTO_VIZ}" == "1" ]]; then
    viz_cmd=(
      Rscript "${VIZ_SCRIPT}"
      "--fit_dir=${run_dir}"
      "--data_dir=${DATA_DIR}"
      "--report_dt=${VIZ_REPORT_DT}"
      "--top_n=${VIZ_TOP_N}"
      "--n_cores=1"
    )
    echo "[$(date '+%F %T')] seed=${seed}: viz start"
    echo "[$(date '+%F %T')] seed=${seed}: viz_log=${viz_log}"
    echo "Viz command: ${viz_cmd[*]}"
    "${viz_cmd[@]}" 2>&1 | tee "${viz_log}"
    echo "[$(date '+%F %T')] seed=${seed}: viz done"
  fi
done

echo "All done. Run directory: ${RUN_DIR}"
