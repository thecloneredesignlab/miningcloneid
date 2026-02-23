#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FIT_SCRIPT="${SCRIPT_DIR}/fit_invivo_model_buffering_align_with_Richard_tmb_hierarchical_chain.R"

usage() {
  cat <<'EOF'
Usage:
  bash run_fit_invivo_model_buffering_align_with_Richard_tmb_hierarchical_chain.sh [--key=value ...]

Examples:
  bash run_fit_invivo_model_buffering_align_with_Richard_tmb_hierarchical_chain.sh \
    --seeds_csv=1,2,3 \
    --n_cores=24 \
    --run_prefix=fit_richard_tmb_hier_chain

Supported options:
  out_root, run_prefix, data_dir, seeds_csv, n_cores, max_scenarios
  fit_treatment, dose_zero_only, paired_only, truncate_at_treatment, ploidy_at_harvest
  w_burden_chain, w_ploidy_chain
  loss_rescale, loss_scale_burden, loss_scale_ploidy, loss_scale_eps
  c_vol_2N_mm3, burden_log_eps, huber_k_burden_log
  n_alt_iter, n_starts_local, n_starts_global, maxit_local, maxit_global
  use_deoptim_local, use_deoptim_global
  deoptim_itermax_local, deoptim_np_local
  deoptim_itermax_global, deoptim_np_global, deoptim_trace
  lambda_shrink, tau_floor, tmb_tau_min, tmb_log_tau_prior_sd, tmb_maxit, tmb_rebuild
  select_rule
  N_MIN, N_MAX, N_UNIT, dt, O2, K
EOF
}

parse_cli_args() {
  for arg in "$@"; do
    case "${arg}" in
      -h|--help)
        usage
        exit 0
        ;;
      --*=*)
        local key="${arg%%=*}"
        local val="${arg#*=}"
        key="${key#--}"
        local env_key
        env_key="$(echo "${key}" | tr '[:lower:]-' '[:upper:]_')"
        case "${env_key}" in
          OUT_ROOT|RUN_PREFIX|DATA_DIR|SEEDS_CSV|N_CORES|MAX_SCENARIOS|FIT_TREATMENT|DOSE_ZERO_ONLY|PAIRED_ONLY|TRUNCATE_AT_TREATMENT|PLOIDY_AT_HARVEST|W_BURDEN_CHAIN|W_PLOIDY_CHAIN|LOSS_RESCALE|LOSS_SCALE_BURDEN|LOSS_SCALE_PLOIDY|LOSS_SCALE_EPS|C_VOL_2N_MM3|BURDEN_LOG_EPS|HUBER_K_BURDEN_LOG|N_ALT_ITER|N_STARTS_LOCAL|N_STARTS_GLOBAL|MAXIT_LOCAL|MAXIT_GLOBAL|USE_DEOPTIM_LOCAL|USE_DEOPTIM_GLOBAL|DEOPTIM_ITERMAX_LOCAL|DEOPTIM_NP_LOCAL|DEOPTIM_ITERMAX_GLOBAL|DEOPTIM_NP_GLOBAL|DEOPTIM_TRACE|LAMBDA_SHRINK|TAU_FLOOR|TMB_TAU_MIN|TMB_LOG_TAU_PRIOR_SD|TMB_MAXIT|TMB_REBUILD|SELECT_RULE|N_MIN|N_MAX|N_UNIT|DT|O2|K)
            export "${env_key}=${val}"
            ;;
          *)
            echo "ERROR: unknown option --${key}" >&2
            usage
            exit 1
            ;;
        esac
        ;;
      *)
        echo "ERROR: unsupported argument '${arg}'. Use --key=value format." >&2
        usage
        exit 1
        ;;
    esac
  done
}

parse_cli_args "$@"

if [[ ! -f "${FIT_SCRIPT}" ]]; then
  echo "ERROR: fit script not found: ${FIT_SCRIPT}" >&2
  exit 1
fi

OUT_ROOT="${OUT_ROOT:-${SCRIPT_DIR}/../results}"
RUN_PREFIX="${RUN_PREFIX:-fit_invivo_model_buffering_align_with_Richard_tmb_hierarchical_chain}"
DATA_DIR="${DATA_DIR:-}"
SEEDS_CSV="${SEEDS_CSV:-1}"
N_CORES="${N_CORES:-}"
MAX_SCENARIOS="${MAX_SCENARIOS:-}"

FIT_TREATMENT="${FIT_TREATMENT:-FALSE}"
DOSE_ZERO_ONLY="${DOSE_ZERO_ONLY:-TRUE}"
PAIRED_ONLY="${PAIRED_ONLY:-TRUE}"
TRUNCATE_AT_TREATMENT="${TRUNCATE_AT_TREATMENT:-FALSE}"
PLOIDY_AT_HARVEST="${PLOIDY_AT_HARVEST:-TRUE}"

W_BURDEN_CHAIN="${W_BURDEN_CHAIN:-1,0.8,0.6,0.4,0.2,0.175,0.15,0.1,0.05,0}"
W_PLOIDY_CHAIN="${W_PLOIDY_CHAIN:-0,0.2,0.4,0.6,0.8,0.825,0.85,0.9,0.95,1}"

LOSS_RESCALE="${LOSS_RESCALE:-TRUE}"
LOSS_SCALE_BURDEN="${LOSS_SCALE_BURDEN:-}"
LOSS_SCALE_PLOIDY="${LOSS_SCALE_PLOIDY:-}"
LOSS_SCALE_EPS="${LOSS_SCALE_EPS:-1e-8}"
C_VOL_2N_MM3="${C_VOL_2N_MM3:-}"
BURDEN_LOG_EPS="${BURDEN_LOG_EPS:-}"
HUBER_K_BURDEN_LOG="${HUBER_K_BURDEN_LOG:-}"

N_ALT_ITER="${N_ALT_ITER:-3}"
N_STARTS_LOCAL="${N_STARTS_LOCAL:-6}"
N_STARTS_GLOBAL="${N_STARTS_GLOBAL:-12}"
MAXIT_LOCAL="${MAXIT_LOCAL:-2500}"
MAXIT_GLOBAL="${MAXIT_GLOBAL:-6000}"
USE_DEOPTIM_LOCAL="${USE_DEOPTIM_LOCAL:-TRUE}"
USE_DEOPTIM_GLOBAL="${USE_DEOPTIM_GLOBAL:-TRUE}"
DEOPTIM_ITERMAX_LOCAL="${DEOPTIM_ITERMAX_LOCAL:-120}"
DEOPTIM_NP_LOCAL="${DEOPTIM_NP_LOCAL:-80}"
DEOPTIM_ITERMAX_GLOBAL="${DEOPTIM_ITERMAX_GLOBAL:-180}"
DEOPTIM_NP_GLOBAL="${DEOPTIM_NP_GLOBAL:-160}"
DEOPTIM_TRACE="${DEOPTIM_TRACE:-TRUE}"
LAMBDA_SHRINK="${LAMBDA_SHRINK:-1.0}"
TAU_FLOOR="${TAU_FLOOR:-0.05}"
TMB_TAU_MIN="${TMB_TAU_MIN:-1e-3}"
TMB_LOG_TAU_PRIOR_SD="${TMB_LOG_TAU_PRIOR_SD:-2.0}"
TMB_MAXIT="${TMB_MAXIT:-200}"
TMB_REBUILD="${TMB_REBUILD:-FALSE}"
SELECT_RULE="${SELECT_RULE:-min_objective_data}"

# Optional model/data domain overrides
N_MIN="${N_MIN:-}"
N_MAX="${N_MAX:-}"
N_UNIT="${N_UNIT:-}"
DT="${DT:-}"
O2="${O2:-}"
K="${K:-}"

if [[ -z "${DATA_DIR}" ]]; then
  DATA_DIR="$(cd "${SCRIPT_DIR}/../../data/InVivoData_Gemcitabine" && pwd)"
fi

mkdir -p "${OUT_ROOT}"

trim() {
  local s="${1:-}"
  s="${s#"${s%%[![:space:]]*}"}"
  s="${s%"${s##*[![:space:]]}"}"
  printf '%s' "${s}"
}

IFS=',' read -r -a SEEDS <<< "${SEEDS_CSV}"

echo "Running TMB hierarchical chain fit"
echo "  Fit script: ${FIT_SCRIPT}"
echo "  Data dir:   ${DATA_DIR}"
echo "  Seeds:      ${SEEDS_CSV}"
echo "  Weights:    (${W_BURDEN_CHAIN}) vs (${W_PLOIDY_CHAIN})"
echo "  paired_only=${PAIRED_ONLY}, fit_treatment=${FIT_TREATMENT}, n_cores=${N_CORES:-auto}"
echo "  Solvers:    local=${USE_DEOPTIM_LOCAL}, global=${USE_DEOPTIM_GLOBAL} (DEoptim)"
echo "  Burden obs model args: c_vol_2N_mm3=${C_VOL_2N_MM3:-default}, burden_log_eps=${BURDEN_LOG_EPS:-default}, huber_k_burden_log=${HUBER_K_BURDEN_LOG:-default}"
echo

for seed_raw in "${SEEDS[@]}"; do
  seed="$(trim "${seed_raw}")"
  [[ -z "${seed}" ]] && continue

  out_dir="${OUT_ROOT}/${RUN_PREFIX}_seed${seed}"
  mkdir -p "${out_dir}"
  log_file="${out_dir}/run_realtime.log"
  cmd=(
    Rscript "${FIT_SCRIPT}"
    "--seed=${seed}"
    "--out_dir=${out_dir}"
    "--data_dir=${DATA_DIR}"
    "--fit_treatment=${FIT_TREATMENT}"
    "--dose_zero_only=${DOSE_ZERO_ONLY}"
    "--paired_only=${PAIRED_ONLY}"
    "--truncate_at_treatment=${TRUNCATE_AT_TREATMENT}"
    "--ploidy_at_harvest=${PLOIDY_AT_HARVEST}"
    "--w_burden=${W_BURDEN_CHAIN}"
    "--w_ploidy=${W_PLOIDY_CHAIN}"
    "--loss_rescale=${LOSS_RESCALE}"
    "--loss_scale_eps=${LOSS_SCALE_EPS}"
    "--n_alt_iter=${N_ALT_ITER}"
    "--n_starts_local=${N_STARTS_LOCAL}"
    "--n_starts_global=${N_STARTS_GLOBAL}"
    "--maxit_local=${MAXIT_LOCAL}"
    "--maxit_global=${MAXIT_GLOBAL}"
    "--use_deoptim_local=${USE_DEOPTIM_LOCAL}"
    "--use_deoptim_global=${USE_DEOPTIM_GLOBAL}"
    "--deoptim_itermax_local=${DEOPTIM_ITERMAX_LOCAL}"
    "--deoptim_np_local=${DEOPTIM_NP_LOCAL}"
    "--deoptim_itermax_global=${DEOPTIM_ITERMAX_GLOBAL}"
    "--deoptim_np_global=${DEOPTIM_NP_GLOBAL}"
    "--deoptim_trace=${DEOPTIM_TRACE}"
    "--lambda_shrink=${LAMBDA_SHRINK}"
    "--tau_floor=${TAU_FLOOR}"
    "--tmb_tau_min=${TMB_TAU_MIN}"
    "--tmb_log_tau_prior_sd=${TMB_LOG_TAU_PRIOR_SD}"
    "--tmb_maxit=${TMB_MAXIT}"
    "--tmb_rebuild=${TMB_REBUILD}"
    "--select_rule=${SELECT_RULE}"
  )

  if [[ -n "${N_CORES}" ]]; then
    cmd+=("--n_cores=${N_CORES}")
  fi
  if [[ -n "${MAX_SCENARIOS}" ]]; then
    cmd+=("--max_scenarios=${MAX_SCENARIOS}")
  fi
  if [[ -n "${LOSS_SCALE_BURDEN}" ]]; then
    cmd+=("--loss_scale_burden=${LOSS_SCALE_BURDEN}")
  fi
  if [[ -n "${LOSS_SCALE_PLOIDY}" ]]; then
    cmd+=("--loss_scale_ploidy=${LOSS_SCALE_PLOIDY}")
  fi
  if [[ -n "${C_VOL_2N_MM3}" ]]; then
    cmd+=("--c_vol_2N_mm3=${C_VOL_2N_MM3}")
  fi
  if [[ -n "${BURDEN_LOG_EPS}" ]]; then
    cmd+=("--burden_log_eps=${BURDEN_LOG_EPS}")
  fi
  if [[ -n "${HUBER_K_BURDEN_LOG}" ]]; then
    cmd+=("--huber_k_burden_log=${HUBER_K_BURDEN_LOG}")
  fi
  if [[ -n "${N_MIN}" ]]; then
    cmd+=("--N_MIN=${N_MIN}")
  fi
  if [[ -n "${N_MAX}" ]]; then
    cmd+=("--N_MAX=${N_MAX}")
  fi
  if [[ -n "${N_UNIT}" ]]; then
    cmd+=("--N_UNIT=${N_UNIT}")
  fi
  if [[ -n "${DT}" ]]; then
    cmd+=("--dt=${DT}")
  fi
  if [[ -n "${O2}" ]]; then
    cmd+=("--O2=${O2}")
  fi
  if [[ -n "${K}" ]]; then
    cmd+=("--K=${K}")
  fi

  echo "[$(date '+%F %T')] seed=${seed}: start"
  echo "[$(date '+%F %T')] seed=${seed}: realtime_log=${log_file}"
  if command -v stdbuf >/dev/null 2>&1; then
    run_cmd=(stdbuf -oL -eL "${cmd[@]}")
  else
    run_cmd=("${cmd[@]}")
  fi
  {
    echo "[$(date '+%F %T')] seed=${seed}: start"
    if command -v stdbuf >/dev/null 2>&1; then
      echo "[$(date '+%F %T')] seed=${seed}: line_buffering=TRUE (stdbuf)"
    else
      echo "[$(date '+%F %T')] seed=${seed}: line_buffering=FALSE (stdbuf not found)"
    fi
    printf 'Command: '
    printf '%q ' "${cmd[@]}"
    echo
    "${run_cmd[@]}"
  } 2>&1 | tee -a "${log_file}"
  echo "[$(date '+%F %T')] seed=${seed}: done -> ${out_dir}" | tee -a "${log_file}"
  echo

done

echo "All done."
