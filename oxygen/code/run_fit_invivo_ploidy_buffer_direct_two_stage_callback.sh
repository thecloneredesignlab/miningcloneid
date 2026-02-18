#!/usr/bin/env bash
set -euo pipefail

# Direct launcher for fit_invivo_ploidy_buffer.R using stage-wise warm starts:
# (1,0)->(0.8,0.2)->(0.6,0.4)->(0.4,0.6)->(0.2,0.8)->(0,1)->callback(1,1)

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FIT_SCRIPT="${SCRIPT_DIR}/fit_invivo_ploidy_buffer.R"

usage() {
  cat <<'EOF'
Usage:
  bash run_fit_invivo_ploidy_buffer_direct_two_stage_callback.sh [--key=value ...]

Examples:
  bash run_fit_invivo_ploidy_buffer_direct_two_stage_callback.sh \
    --seeds_csv=1,2,3 \
    --n_cores=34 \
    --out_root=/path/to/results \
    --run_prefix=fit_stagewise_chain

  # Custom weight chain + callback
  bash run_fit_invivo_ploidy_buffer_direct_two_stage_callback.sh \
    --w_burden_chain=1,0.8,0.6,0.4,0.2,0 \
    --w_ploidy_chain=0,0.2,0.4,0.6,0.8,1 \
    --callback_w_burden=1 \
    --callback_w_ploidy=1

Supported --key=value options:
  out_root, run_prefix, seeds_csv, k, n_cores, max_scenarios
  pass_itermax, callback_itermax, np
  pass_n_starts, callback_n_starts
  pass_optim_maxit, callback_optim_maxit
  use_deoptim, deoptim_parallel
  dose_zero_only, truncate_at_treatment, ploidy_at_harvest
  loss_rescale, loss_scale_burden, loss_scale_ploidy, loss_scale_eps
  w_burden_chain, w_ploidy_chain, callback_w_burden, callback_w_ploidy
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
          OUT_ROOT|RUN_PREFIX|SEEDS_CSV|K|N_CORES|MAX_SCENARIOS|PASS_ITERMAX|CALLBACK_ITERMAX|NP|PASS_N_STARTS|CALLBACK_N_STARTS|PASS_OPTIM_MAXIT|CALLBACK_OPTIM_MAXIT|USE_DEOPTIM|DEOPTIM_PARALLEL|DOSE_ZERO_ONLY|TRUNCATE_AT_TREATMENT|PLOIDY_AT_HARVEST|LOSS_RESCALE|LOSS_SCALE_BURDEN|LOSS_SCALE_PLOIDY|LOSS_SCALE_EPS|W_BURDEN_CHAIN|W_PLOIDY_CHAIN|CALLBACK_W_BURDEN|CALLBACK_W_PLOIDY|OMP_NUM_THREADS|OPENBLAS_NUM_THREADS|MKL_NUM_THREADS|VECLIB_MAXIMUM_THREADS)
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

# ----------------------------
# Defaults (override via env)
# ----------------------------
OUT_ROOT="${OUT_ROOT:-${SCRIPT_DIR}/../results}"
RUN_PREFIX="${RUN_PREFIX:-fit_invivo_stagewise_chain}"

SEEDS_CSV="${SEEDS_CSV:-1,2,3}"
K="${K:-1e12}"
N_CORES="${N_CORES:-}"
MAX_SCENARIOS="${MAX_SCENARIOS:-}"

PASS_ITERMAX="${PASS_ITERMAX:-120}"
CALLBACK_ITERMAX="${CALLBACK_ITERMAX:-220}"
NP="${NP:-140}"
PASS_N_STARTS="${PASS_N_STARTS:-40}"
CALLBACK_N_STARTS="${CALLBACK_N_STARTS:-80}"
PASS_OPTIM_MAXIT="${PASS_OPTIM_MAXIT:-8000}"
CALLBACK_OPTIM_MAXIT="${CALLBACK_OPTIM_MAXIT:-14000}"

USE_DEOPTIM="${USE_DEOPTIM:-FALSE}"
DEOPTIM_PARALLEL="${DEOPTIM_PARALLEL:-FALSE}"

DOSE_ZERO_ONLY="${DOSE_ZERO_ONLY:-TRUE}"
TRUNCATE_AT_TREATMENT="${TRUNCATE_AT_TREATMENT:-FALSE}"
PLOIDY_AT_HARVEST="${PLOIDY_AT_HARVEST:-TRUE}"
LOSS_RESCALE="${LOSS_RESCALE:-TRUE}"
LOSS_SCALE_BURDEN="${LOSS_SCALE_BURDEN:-}"
LOSS_SCALE_PLOIDY="${LOSS_SCALE_PLOIDY:-}"
LOSS_SCALE_EPS="${LOSS_SCALE_EPS:-1e-8}"

# Progressive chain weights (same length, comma-separated)
W_BURDEN_CHAIN="${W_BURDEN_CHAIN:-1,0.8,0.6,0.4,0.2,0}"
W_PLOIDY_CHAIN="${W_PLOIDY_CHAIN:-0,0.2,0.4,0.6,0.8,1}"
CALLBACK_W_BURDEN="${CALLBACK_W_BURDEN:-1}"
CALLBACK_W_PLOIDY="${CALLBACK_W_PLOIDY:-1}"

# Avoid CPU oversubscription on HPC by default.
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export VECLIB_MAXIMUM_THREADS="${VECLIB_MAXIMUM_THREADS:-1}"

trim() {
  local s="${1:-}"
  s="${s#"${s%%[![:space:]]*}"}"
  s="${s%"${s##*[![:space:]]}"}"
  printf '%s' "${s}"
}

sanitize_tag() {
  local s="${1:-NA}"
  s="${s// /}"
  s="${s//,/}"
  s="${s//./p}"
  s="${s//-/m}"
  printf '%s' "${s}"
}

run_fit_cmd() {
  local label="$1"
  local out_dir="$2"
  local wb="$3"
  local wp="$4"
  local init_tsv="${5:-}"
  local itermax="$6"
  local n_starts="$7"
  local optim_maxit="$8"
  local seed="$9"
  shift 9

  local cmd=(
    Rscript "${FIT_SCRIPT}"
    "--two_stage=FALSE"
    "--w_burden=${wb}"
    "--w_ploidy=${wp}"
    "--K=${K}"
    "--use_deoptim=${USE_DEOPTIM}"
    "--deoptim_parallel=${DEOPTIM_PARALLEL}"
    "--truncate_at_treatment=${TRUNCATE_AT_TREATMENT}"
    "--dose_zero_only=${DOSE_ZERO_ONLY}"
    "--ploidy_at_harvest=${PLOIDY_AT_HARVEST}"
    "--loss_rescale=${LOSS_RESCALE}"
    "--loss_scale_eps=${LOSS_SCALE_EPS}"
    "--itermax=${itermax}"
    "--NP=${NP}"
    "--n_starts=${n_starts}"
    "--optim_maxit=${optim_maxit}"
    "--seed=${seed}"
    "--out_dir=${out_dir}"
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
  if [[ -n "${init_tsv}" ]]; then
    cmd+=("--init_params_tsv=${init_tsv}")
  fi

  mkdir -p "${out_dir}"
  local log_file="${out_dir}.log"
  {
    echo "[$(date '+%F %T')] ${label}"
    echo "Command:"
    printf '  %q' "${cmd[@]}"
    echo
    echo
    "${cmd[@]}"
  } 2>&1 | tee "${log_file}"
}

append_callback_metrics() {
  local callback_dir="$1"
  local seed="$2"
  local metrics_file="$3"
  Rscript - "${callback_dir}" "${seed}" >> "${metrics_file}" <<'RS'
args <- commandArgs(trailingOnly = TRUE)
run_dir <- args[[1]]
seed <- args[[2]]

sum_path <- file.path(run_dir, "fit_summary.tsv")
if (!file.exists(sum_path)) {
  cat(seed, "\tNA\tNA\tNA\tNA\tNA\n", sep = "")
  quit(save = "no")
}
sum_df <- read.delim(sum_path, stringsAsFactors = FALSE, check.names = FALSE)
vals <- setNames(sum_df$value, sum_df$metric)

obj <- as.numeric(vals[["objective"]])
obj_b <- as.numeric(vals[["objective_burden"]])
obj_p <- as.numeric(vals[["objective_ploidy"]])

rmse_4n <- NA_real_
bf <- file.path(run_dir, "burden_fit.tsv")
if (file.exists(bf)) {
  b <- read.delim(bf, stringsAsFactors = FALSE, check.names = FALSE)
  idx <- which(b$cohort == "4N" & is.finite(b$pred_norm) & is.finite(b$obs_norm))
  if (length(idx) > 0) {
    rmse_4n <- sqrt(mean((b$pred_norm[idx] - b$obs_norm[idx])^2))
  }
}

nll_4n <- NA_real_
pf <- file.path(run_dir, "terminal_ploidy_fit.tsv")
if (file.exists(pf)) {
  p <- read.delim(pf, stringsAsFactors = FALSE, check.names = FALSE)
  p <- p[p$cohort == "4N", , drop = FALSE]
  if (nrow(p) > 0) {
    key <- interaction(p$harvest, p$dose, drop = TRUE)
    split_idx <- split(seq_len(nrow(p)), key)
    nll_each <- vapply(split_idx, function(ix) {
      x <- p[ix, , drop = FALSE]
      obs <- as.numeric(x$obs_count)
      if (sum(obs) <= 0) return(NA_real_)
      pred <- pmax(as.numeric(x$pred_fraction), 1e-12)
      -sum(obs * log(pred)) / sum(obs)
    }, numeric(1))
    nll_4n <- mean(nll_each, na.rm = TRUE)
    if (!is.finite(nll_4n)) nll_4n <- NA_real_
  }
}

cat(
  seed, "\t",
  format(obj, digits = 10), "\t",
  format(obj_b, digits = 10), "\t",
  format(obj_p, digits = 10), "\t",
  format(rmse_4n, digits = 10), "\t",
  format(nll_4n, digits = 10), "\n",
  sep = ""
)
RS
}

mkdir -p "${OUT_ROOT}"
IFS=',' read -r -a SEEDS <<< "${SEEDS_CSV}"
IFS=',' read -r -a WB_CHAIN <<< "${W_BURDEN_CHAIN}"
IFS=',' read -r -a WP_CHAIN <<< "${W_PLOIDY_CHAIN}"

if [[ "${#WB_CHAIN[@]}" -eq 0 || "${#WP_CHAIN[@]}" -eq 0 ]]; then
  echo "ERROR: W_BURDEN_CHAIN / W_PLOIDY_CHAIN cannot be empty." >&2
  exit 1
fi
if [[ "${#WB_CHAIN[@]}" -ne "${#WP_CHAIN[@]}" ]]; then
  echo "ERROR: W_BURDEN_CHAIN and W_PLOIDY_CHAIN must have the same length." >&2
  exit 1
fi

METRICS_TSV="${OUT_ROOT}/${RUN_PREFIX}_callback_metrics.tsv"
echo -e "seed\tobjective\tobjective_burden\tobjective_ploidy\trmse_4N_burden\tmean_nll_4N_ploidy" > "${METRICS_TSV}"

echo "Running stage-wise warm-start chain + callback"
echo "  Seeds: ${SEEDS_CSV}"
echo "  Chain: (${W_BURDEN_CHAIN}) vs (${W_PLOIDY_CHAIN})"
echo "  Callback: w_burden=${CALLBACK_W_BURDEN}, w_ploidy=${CALLBACK_W_PLOIDY}"
echo "  Loss rescale: ${LOSS_RESCALE} (scale_b=${LOSS_SCALE_BURDEN:-auto}, scale_p=${LOSS_SCALE_PLOIDY:-auto}, eps=${LOSS_SCALE_EPS})"
echo "  Out root: ${OUT_ROOT}"
echo

for seed_raw in "${SEEDS[@]}"; do
  seed="$(trim "${seed_raw}")"
  [[ -z "${seed}" ]] && continue

  run_root="${OUT_ROOT}/${RUN_PREFIX}_seed${seed}"
  mkdir -p "${run_root}"
  init_tsv=""

  echo "[$(date '+%F %T')] seed=${seed}: stage-wise chain start"
  for i in "${!WB_CHAIN[@]}"; do
    pass_id=$((i + 1))
    wb="$(trim "${WB_CHAIN[$i]}")"
    wp="$(trim "${WP_CHAIN[$i]}")"

    wb_tag="$(sanitize_tag "${wb}")"
    wp_tag="$(sanitize_tag "${wp}")"
    step_dir="${run_root}/step$(printf '%02d' "${pass_id}")_wb${wb_tag}_wp${wp_tag}"
    step_label="seed=${seed} step${pass_id}/${#WB_CHAIN[@]} w=(${wb},${wp})"

    run_fit_cmd "${step_label}" "${step_dir}" "${wb}" "${wp}" "${init_tsv}" "${PASS_ITERMAX}" "${PASS_N_STARTS}" "${PASS_OPTIM_MAXIT}" "${seed}"
    init_tsv="${step_dir}/fit_parameter_stages.tsv"
    if [[ ! -f "${init_tsv}" ]]; then
      echo "ERROR: missing warm-start file after step ${pass_id}: ${init_tsv}" >&2
      exit 1
    fi
  done

  callback_dir="${run_root}/callback_equal"
  callback_label="seed=${seed} callback w=(${CALLBACK_W_BURDEN},${CALLBACK_W_PLOIDY})"
  run_fit_cmd "${callback_label}" "${callback_dir}" "${CALLBACK_W_BURDEN}" "${CALLBACK_W_PLOIDY}" "${init_tsv}" "${CALLBACK_ITERMAX}" "${CALLBACK_N_STARTS}" "${CALLBACK_OPTIM_MAXIT}" "${seed}"
  append_callback_metrics "${callback_dir}" "${seed}" "${METRICS_TSV}"
done

echo "Done."
echo "Metrics summary: ${METRICS_TSV}"
