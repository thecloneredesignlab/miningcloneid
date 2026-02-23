#!/usr/bin/env bash
set -euo pipefail

# Direct launcher for fit_invivo_model_buffering_align_with_Richard.R
# using iterative stage-wise warm starts:
# (1,0)->(0.8,0.2)->...->(0,1), then callback (1,1).

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FIT_SCRIPT="${SCRIPT_DIR}/fit_invivo_model_buffering_align_with_Richard.R"

usage() {
  cat <<'EOF'
Usage:
  bash run_fit_invivo_model_buffering_align_with_Richard_it_stagewise_chain_callback.sh [--key=value ...]

Examples:
  bash run_fit_invivo_model_buffering_align_with_Richard_it_stagewise_chain_callback.sh \
    --seeds_csv=1,2,3 \
    --n_cores=24 \
    --out_root=/path/to/results \
    --run_prefix=fit_invivo_model_buffering_align_with_Richard_it_stagewise_chain

  # Custom chain + callback
  bash run_fit_invivo_model_buffering_align_with_Richard_it_stagewise_chain_callback.sh \
    --w_burden_chain=1,0.8,0.6,0.4,0.2,0.175,0.15,0.1,0.05,0 \
    --w_ploidy_chain=0,0.2,0.4,0.6,0.8,0.825,0.85,0.9,0.95,1 \
    --callback_w_burden=1 \
    --callback_w_ploidy=1

Supported --key=value options:
  out_root, run_prefix, data_dir, seeds_csv, k, n_cores, max_scenarios
  pass_itermax, callback_itermax, np
  pass_n_starts, callback_n_starts
  pass_optim_maxit, callback_optim_maxit
  use_deoptim, deoptim_parallel
  fit_treatment, dose_zero_only, paired_only, truncate_at_treatment, ploidy_at_harvest
  loss_rescale, loss_scale_burden, loss_scale_ploidy, loss_scale_eps
  c_vol_2N_mm3, burden_log_eps, huber_k_burden_log
  w_burden_chain, w_ploidy_chain, callback_w_burden, callback_w_ploidy
  auto_tune_iters
  resume_from_pass, resume_init_tsv_template, resume_skip_existing

Resume options:
  resume_from_pass: 1-based pass index to start from (default: 1).
  resume_init_tsv_template: warm-start file path template; supports {seed} placeholder.
  resume_skip_existing: TRUE/FALSE. If TRUE, existing step/callback outputs are reused.
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
          OUT_ROOT|RUN_PREFIX|DATA_DIR|SEEDS_CSV|K|N_CORES|MAX_SCENARIOS|PASS_ITERMAX|CALLBACK_ITERMAX|NP|PASS_N_STARTS|CALLBACK_N_STARTS|PASS_OPTIM_MAXIT|CALLBACK_OPTIM_MAXIT|USE_DEOPTIM|DEOPTIM_PARALLEL|FIT_TREATMENT|DOSE_ZERO_ONLY|PAIRED_ONLY|TRUNCATE_AT_TREATMENT|PLOIDY_AT_HARVEST|LOSS_RESCALE|LOSS_SCALE_BURDEN|LOSS_SCALE_PLOIDY|LOSS_SCALE_EPS|C_VOL_2N_MM3|BURDEN_LOG_EPS|HUBER_K_BURDEN_LOG|W_BURDEN_CHAIN|W_PLOIDY_CHAIN|CALLBACK_W_BURDEN|CALLBACK_W_PLOIDY|AUTO_TUNE_ITERS|RESUME_FROM_PASS|RESUME_INIT_TSV_TEMPLATE|RESUME_SKIP_EXISTING|OMP_NUM_THREADS|OPENBLAS_NUM_THREADS|MKL_NUM_THREADS|VECLIB_MAXIMUM_THREADS)
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
RUN_PREFIX="${RUN_PREFIX:-fit_invivo_model_buffering_align_with_Richard_it_stagewise_chain}"
DATA_DIR="${DATA_DIR:-}"

SEEDS_CSV="${SEEDS_CSV:-1,2,3}"
K="${K:-1e12}"
N_CORES="${N_CORES:-}"
MAX_SCENARIOS="${MAX_SCENARIOS:-}"

# If empty, values will be auto-estimated from input files.
PASS_ITERMAX="${PASS_ITERMAX:-}"
CALLBACK_ITERMAX="${CALLBACK_ITERMAX:-}"
NP="${NP:-}"
PASS_N_STARTS="${PASS_N_STARTS:-}"
CALLBACK_N_STARTS="${CALLBACK_N_STARTS:-}"
PASS_OPTIM_MAXIT="${PASS_OPTIM_MAXIT:-}"
CALLBACK_OPTIM_MAXIT="${CALLBACK_OPTIM_MAXIT:-}"
AUTO_TUNE_ITERS="${AUTO_TUNE_ITERS:-TRUE}"

USE_DEOPTIM="${USE_DEOPTIM:-FALSE}"
DEOPTIM_PARALLEL="${DEOPTIM_PARALLEL:-FALSE}"
FIT_TREATMENT="${FIT_TREATMENT:-FALSE}"

DOSE_ZERO_ONLY="${DOSE_ZERO_ONLY:-TRUE}"
PAIRED_ONLY="${PAIRED_ONLY:-TRUE}"
TRUNCATE_AT_TREATMENT="${TRUNCATE_AT_TREATMENT:-FALSE}"
PLOIDY_AT_HARVEST="${PLOIDY_AT_HARVEST:-TRUE}"
LOSS_RESCALE="${LOSS_RESCALE:-TRUE}"
LOSS_SCALE_BURDEN="${LOSS_SCALE_BURDEN:-}"
LOSS_SCALE_PLOIDY="${LOSS_SCALE_PLOIDY:-}"
LOSS_SCALE_EPS="${LOSS_SCALE_EPS:-1e-8}"
C_VOL_2N_MM3="${C_VOL_2N_MM3:-}"
BURDEN_LOG_EPS="${BURDEN_LOG_EPS:-}"
HUBER_K_BURDEN_LOG="${HUBER_K_BURDEN_LOG:-}"

# Required chain from request
W_BURDEN_CHAIN="${W_BURDEN_CHAIN:-1,0.8,0.6,0.4,0.2,0.175,0.15,0.1,0.05,0}"
W_PLOIDY_CHAIN="${W_PLOIDY_CHAIN:-0,0.2,0.4,0.6,0.8,0.825,0.85,0.9,0.95,1}"
CALLBACK_W_BURDEN="${CALLBACK_W_BURDEN:-1}"
CALLBACK_W_PLOIDY="${CALLBACK_W_PLOIDY:-1}"

RESUME_FROM_PASS="${RESUME_FROM_PASS:-1}"
RESUME_INIT_TSV_TEMPLATE="${RESUME_INIT_TSV_TEMPLATE:-}"
RESUME_SKIP_EXISTING="${RESUME_SKIP_EXISTING:-FALSE}"

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

is_true() {
  local x="${1:-}"
  x="$(echo "${x}" | tr '[:upper:]' '[:lower:]')"
  [[ "${x}" == "1" || "${x}" == "true" || "${x}" == "t" || "${x}" == "yes" || "${x}" == "y" ]]
}

step_dir_for_pass() {
  local run_root="$1"
  local pass_id="$2"
  local wb="$3"
  local wp="$4"
  local wb_tag
  local wp_tag
  wb_tag="$(sanitize_tag "${wb}")"
  wp_tag="$(sanitize_tag "${wp}")"
  printf '%s/step%02d_wb%s_wp%s' "${run_root}" "${pass_id}" "${wb_tag}" "${wp_tag}"
}

resolve_seed_template_path() {
  local tmpl="$1"
  local seed="$2"
  local out="${tmpl//\{seed\}/${seed}}"
  out="${out//__SEED__/${seed}}"
  printf '%s' "${out}"
}

estimate_runtime_defaults() {
  local fit_script="$1"
  local data_dir="$2"
  Rscript - "${fit_script}" "${data_dir}" <<'RS'
args <- commandArgs(trailingOnly = TRUE)
fit_script <- args[[1]]
data_dir <- args[[2]]

source(fit_script)
paired_only_env <- tolower(Sys.getenv("PAIRED_ONLY", unset = "TRUE")) %in% c("1", "true", "t", "yes", "y")

cfg <- list(
  model_path = file.path(dirname(fit_script), "model_buffering_align_with_Richard.R"),
  N_UNIT = 22L, N_MIN = 22L, N_MAX = 154L,
  DT = 0.5, O2_fixed = 1.0, K = 1e12, crowding = "logistic",
  init_total_size = 1e6, dose_ref = 30, tx_mult_min = 0.05, min_pop = 1e-12,
  huber_k = 0.1, w_burden = 1, w_ploidy = 1, w_burden_schedule = 1, w_ploidy_schedule = 1, n_weight_passes = 1L,
  loss_rescale = TRUE, loss_scale_burden = NA_real_, loss_scale_ploidy = NA_real_, loss_scale_eps = 1e-8, loss_scale_source = "unset",
  optim_trace = TRUE, optim_trace_every = 1L, eps_prob = 1e-12, trace_obj = FALSE,
  fit_full_pmis = FALSE, fit_treatment = FALSE,
  dose_zero_only = TRUE, paired_only = paired_only_env, truncate_at_treatment = FALSE, ploidy_at_harvest = TRUE,
  two_stage = FALSE, stage1_w_burden = 1.0, stage1_w_ploidy = 0.0, stage2_w_burden = 0.0, stage2_w_ploidy = 1.0,
  use_deoptim = FALSE, deoptim_parallel = FALSE, itermax = 40L, NP = 80L, n_cores = 1L, seed = 1L,
  max_scenarios = Inf
)

dt_path <- file.path(data_dir, "dt_Gem_VT_20260209_v5.xlsx")
ploidy_path <- file.path(data_dir, "all_ploidy.tsv")
sc <- prepare_data(dt_path, ploidy_path, cfg)

n_sc <- length(sc)
n_pl <- sum(vapply(sc, function(s) length(s$ploidy_obs_N) > 0, logical(1)))
n_pts <- sum(vapply(sc, function(s) length(s$obs_days), integer(1)))

# Complexity heuristic for Richard model (9 params by default, plus warm-start chain).
score <- n_sc + 0.5 * n_pts + 2.0 * n_pl
if (!is.finite(score)) score <- 80

if (score <= 60) {
  pass_iter <- 120L
  callback_iter <- 220L
  np <- 140L
  pass_starts <- 30L
  callback_starts <- 60L
  pass_maxit <- 7000L
  callback_maxit <- 12000L
} else if (score <= 120) {
  pass_iter <- 150L
  callback_iter <- 280L
  np <- 180L
  pass_starts <- 45L
  callback_starts <- 90L
  pass_maxit <- 9000L
  callback_maxit <- 15000L
} else {
  pass_iter <- 180L
  callback_iter <- 340L
  np <- 220L
  pass_starts <- 60L
  callback_starts <- 120L
  pass_maxit <- 12000L
  callback_maxit <- 20000L
}

cat(sprintf("N_SCENARIOS=%d\n", n_sc))
cat(sprintf("N_PLOIDY_SCENARIOS=%d\n", n_pl))
cat(sprintf("N_BURDEN_POINTS=%d\n", n_pts))
cat(sprintf("PASS_ITERMAX=%d\n", pass_iter))
cat(sprintf("CALLBACK_ITERMAX=%d\n", callback_iter))
cat(sprintf("NP=%d\n", np))
cat(sprintf("PASS_N_STARTS=%d\n", pass_starts))
cat(sprintf("CALLBACK_N_STARTS=%d\n", callback_starts))
cat(sprintf("PASS_OPTIM_MAXIT=%d\n", pass_maxit))
cat(sprintf("CALLBACK_OPTIM_MAXIT=%d\n", callback_maxit))
RS
}

if [[ -z "${DATA_DIR}" ]]; then
  DATA_DIR="$(cd "${SCRIPT_DIR}/../../data/InVivoData_Gemcitabine" && pwd)"
fi

if is_true "${AUTO_TUNE_ITERS}"; then
  echo "[$(date '+%F %T')] Estimating runtime defaults from input files under: ${DATA_DIR}"
  if est_out="$(estimate_runtime_defaults "${FIT_SCRIPT}" "${DATA_DIR}")"; then
    N_SCENARIOS=""
    N_PLOIDY_SCENARIOS=""
    N_BURDEN_POINTS=""
    while IFS='=' read -r key val; do
      case "${key}" in
        N_SCENARIOS) N_SCENARIOS="${val}" ;;
        N_PLOIDY_SCENARIOS) N_PLOIDY_SCENARIOS="${val}" ;;
        N_BURDEN_POINTS) N_BURDEN_POINTS="${val}" ;;
        PASS_ITERMAX) [[ -z "${PASS_ITERMAX}" ]] && PASS_ITERMAX="${val}" ;;
        CALLBACK_ITERMAX) [[ -z "${CALLBACK_ITERMAX}" ]] && CALLBACK_ITERMAX="${val}" ;;
        NP) [[ -z "${NP}" ]] && NP="${val}" ;;
        PASS_N_STARTS) [[ -z "${PASS_N_STARTS}" ]] && PASS_N_STARTS="${val}" ;;
        CALLBACK_N_STARTS) [[ -z "${CALLBACK_N_STARTS}" ]] && CALLBACK_N_STARTS="${val}" ;;
        PASS_OPTIM_MAXIT) [[ -z "${PASS_OPTIM_MAXIT}" ]] && PASS_OPTIM_MAXIT="${val}" ;;
        CALLBACK_OPTIM_MAXIT) [[ -z "${CALLBACK_OPTIM_MAXIT}" ]] && CALLBACK_OPTIM_MAXIT="${val}" ;;
      esac
    done <<< "${est_out}"
    echo "  Data complexity: scenarios=${N_SCENARIOS:-NA}, ploidy_scenarios=${N_PLOIDY_SCENARIOS:-NA}, burden_points=${N_BURDEN_POINTS:-NA}"
  else
    echo "WARNING: runtime estimation failed; falling back to static defaults." >&2
  fi
fi

# Static fallback if still unset.
PASS_ITERMAX="${PASS_ITERMAX:-150}"
CALLBACK_ITERMAX="${CALLBACK_ITERMAX:-280}"
NP="${NP:-180}"
PASS_N_STARTS="${PASS_N_STARTS:-45}"
CALLBACK_N_STARTS="${CALLBACK_N_STARTS:-90}"
PASS_OPTIM_MAXIT="${PASS_OPTIM_MAXIT:-9000}"
CALLBACK_OPTIM_MAXIT="${CALLBACK_OPTIM_MAXIT:-15000}"

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
    "--fit_treatment=${FIT_TREATMENT}"
    "--w_burden=${wb}"
    "--w_ploidy=${wp}"
    "--K=${K}"
    "--use_deoptim=${USE_DEOPTIM}"
    "--deoptim_parallel=${DEOPTIM_PARALLEL}"
    "--truncate_at_treatment=${TRUNCATE_AT_TREATMENT}"
    "--dose_zero_only=${DOSE_ZERO_ONLY}"
    "--paired_only=${PAIRED_ONLY}"
    "--ploidy_at_harvest=${PLOIDY_AT_HARVEST}"
    "--loss_rescale=${LOSS_RESCALE}"
    "--loss_scale_eps=${LOSS_SCALE_EPS}"
    "--itermax=${itermax}"
    "--NP=${NP}"
    "--n_starts=${n_starts}"
    "--optim_maxit=${optim_maxit}"
    "--seed=${seed}"
    "--out_dir=${out_dir}"
    "--data_dir=${DATA_DIR}"
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

RESUME_FROM_PASS="$(trim "${RESUME_FROM_PASS}")"
if ! [[ "${RESUME_FROM_PASS}" =~ ^[0-9]+$ ]]; then
  echo "ERROR: RESUME_FROM_PASS must be an integer >= 1." >&2
  exit 1
fi

N_PASS_TOTAL="${#WB_CHAIN[@]}"
if (( RESUME_FROM_PASS < 1 || RESUME_FROM_PASS > N_PASS_TOTAL + 1 )); then
  echo "ERROR: RESUME_FROM_PASS must be in [1, $((N_PASS_TOTAL + 1))]." >&2
  exit 1
fi

METRICS_TSV="${OUT_ROOT}/${RUN_PREFIX}_callback_metrics.tsv"
echo -e "seed\tobjective\tobjective_burden\tobjective_ploidy\trmse_4N_burden\tmean_nll_4N_ploidy" > "${METRICS_TSV}"

echo "Running Richard-model stage-wise warm-start chain + callback"
echo "  Seeds: ${SEEDS_CSV}"
echo "  Chain: (${W_BURDEN_CHAIN}) vs (${W_PLOIDY_CHAIN})"
echo "  Callback: w_burden=${CALLBACK_W_BURDEN}, w_ploidy=${CALLBACK_W_PLOIDY}"
  echo "  Fit treatment: ${FIT_TREATMENT}"
echo "  paired_only: ${PAIRED_ONLY}"
echo "  Loss rescale: ${LOSS_RESCALE} (scale_b=${LOSS_SCALE_BURDEN:-auto}, scale_p=${LOSS_SCALE_PLOIDY:-auto}, eps=${LOSS_SCALE_EPS})"
echo "  Burden obs model args: c_vol_2N_mm3=${C_VOL_2N_MM3:-default}, burden_log_eps=${BURDEN_LOG_EPS:-default}, huber_k_burden_log=${HUBER_K_BURDEN_LOG:-default}"
echo "  Iter settings: pass_itermax=${PASS_ITERMAX}, callback_itermax=${CALLBACK_ITERMAX}, NP=${NP}, pass_n_starts=${PASS_N_STARTS}, callback_n_starts=${CALLBACK_N_STARTS}, pass_optim_maxit=${PASS_OPTIM_MAXIT}, callback_optim_maxit=${CALLBACK_OPTIM_MAXIT}"
echo "  Resume from pass: ${RESUME_FROM_PASS} (skip_existing=${RESUME_SKIP_EXISTING}, init_template=${RESUME_INIT_TSV_TEMPLATE:-NA})"
echo "  Data dir: ${DATA_DIR}"
echo "  Out root: ${OUT_ROOT}"
echo

for seed_raw in "${SEEDS[@]}"; do
  seed="$(trim "${seed_raw}")"
  [[ -z "${seed}" ]] && continue

  run_root="${OUT_ROOT}/${RUN_PREFIX}_seed${seed}"
  mkdir -p "${run_root}"
  init_tsv=""

  if [[ -n "${RESUME_INIT_TSV_TEMPLATE}" ]]; then
    init_tsv="$(resolve_seed_template_path "${RESUME_INIT_TSV_TEMPLATE}" "${seed}")"
    if [[ ! -f "${init_tsv}" ]]; then
      echo "ERROR: seed=${seed} resume_init_tsv_template resolved to missing file: ${init_tsv}" >&2
      exit 1
    fi
  fi

  if (( RESUME_FROM_PASS > 1 )) && [[ -z "${init_tsv}" ]]; then
    prev_pass=$((RESUME_FROM_PASS - 1))
    prev_wb="$(trim "${WB_CHAIN[$((prev_pass - 1))]}")"
    prev_wp="$(trim "${WP_CHAIN[$((prev_pass - 1))]}")"
    prev_step_dir="$(step_dir_for_pass "${run_root}" "${prev_pass}" "${prev_wb}" "${prev_wp}")"
    prev_step_tsv="${prev_step_dir}/fit_parameter_stages.tsv"
    if [[ ! -f "${prev_step_tsv}" ]]; then
      echo "ERROR: seed=${seed} resume requires previous pass output: ${prev_step_tsv}" >&2
      echo "Hint: provide --resume_init_tsv_template=... if resuming from another run root." >&2
      exit 1
    fi
    init_tsv="${prev_step_tsv}"
  fi

  echo "[$(date '+%F %T')] seed=${seed}: stage-wise chain start (resume_from_pass=${RESUME_FROM_PASS})"
  for ((i=RESUME_FROM_PASS-1; i<${#WB_CHAIN[@]}; i++)); do
    pass_id=$((i + 1))
    wb="$(trim "${WB_CHAIN[$i]}")"
    wp="$(trim "${WP_CHAIN[$i]}")"

    step_dir="$(step_dir_for_pass "${run_root}" "${pass_id}" "${wb}" "${wp}")"
    step_label="seed=${seed} step${pass_id}/${#WB_CHAIN[@]} w=(${wb},${wp})"
    step_tsv="${step_dir}/fit_parameter_stages.tsv"

    if is_true "${RESUME_SKIP_EXISTING}" && [[ -f "${step_tsv}" ]]; then
      echo "[$(date '+%F %T')] seed=${seed}: skip existing step${pass_id} -> ${step_tsv}"
      init_tsv="${step_tsv}"
      continue
    fi

    run_fit_cmd "${step_label}" "${step_dir}" "${wb}" "${wp}" "${init_tsv}" "${PASS_ITERMAX}" "${PASS_N_STARTS}" "${PASS_OPTIM_MAXIT}" "${seed}"
    init_tsv="${step_tsv}"
    if [[ ! -f "${init_tsv}" ]]; then
      echo "ERROR: missing warm-start file after step ${pass_id}: ${init_tsv}" >&2
      exit 1
    fi
  done

  callback_dir="${run_root}/callback_equal"
  callback_label="seed=${seed} callback w=(${CALLBACK_W_BURDEN},${CALLBACK_W_PLOIDY})"
  if is_true "${RESUME_SKIP_EXISTING}" && [[ -f "${callback_dir}/fit_summary.tsv" ]]; then
    echo "[$(date '+%F %T')] seed=${seed}: skip existing callback -> ${callback_dir}/fit_summary.tsv"
  else
    run_fit_cmd "${callback_label}" "${callback_dir}" "${CALLBACK_W_BURDEN}" "${CALLBACK_W_PLOIDY}" "${init_tsv}" "${CALLBACK_ITERMAX}" "${CALLBACK_N_STARTS}" "${CALLBACK_OPTIM_MAXIT}" "${seed}"
  fi
  append_callback_metrics "${callback_dir}" "${seed}" "${METRICS_TSV}"
done

echo "Done."
echo "Metrics summary: ${METRICS_TSV}"
