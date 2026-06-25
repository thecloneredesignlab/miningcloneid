#!/usr/bin/env bash
# Unified local O2 runner for in vivo, in vitro, and joint fitting.

set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  bash run_o2_fit.sh --fitting_mode=invivo [options]
  bash run_o2_fit.sh --fitting_mode=invitro [options]
  bash run_o2_fit.sh --fitting_mode=joint --joint_fitting_mode=JOINT [options]
  bash run_o2_fit.sh --fitting_mode=joint --joint_fitting_mode=DIRECT [options]
  bash run_o2_fit.sh --fitting_mode=joint --joint_fitting_mode=MULTI_WARMUP [options]

Required modes:
  --fitting_mode=invivo|invitro|joint
  --joint_fitting_mode=OFF|JOINT|DIRECT
  If fitting_mode=joint and joint_fitting_mode is omitted, DIRECT is used.

Joint mode behavior:
  OFF    Do not run joint fitting. This is forced when fitting_mode is not joint.
  JOINT  Run or reuse in vivo and in vitro single fits, run extra_results,
         select each best seed, then run joint fitting from the selected
         single-fit anchors. Provided best_seed_dir skips that side completely;
         provided run_dir skips fitting but still selects the best seed after
         extra_results; missing sides are run and selected before joint.
  DIRECT Run only the current joint fitter directly from the config.
         SINGLE is accepted as a legacy alias for DIRECT.
  MULTI_WARMUP
         Run or reuse in vivo and in vitro source runs, build a source ratio
         UMAP/cluster manifest, then run one joint fit per selected warm-up
         pair. User-specified best seed dirs are not accepted in this mode.

Common options:
  --project_root=/path/to/repo
  Defaults to the repository containing this runner script.
  --config_path=/path/to/O2_supply_demand.yaml
  --out_root=/path/to/oxygen/results
  --r_module=R/4.4
  --dry_run=TRUE|FALSE
  --run_extra_results=TRUE|FALSE
  --force_extra_results=TRUE|FALSE
  --append_run_prefix_timestamp=TRUE|FALSE

Single-fit options:
  --invivo_run_prefix=name
  --invitro_run_prefix=name
  --invivo_run_dir=/path/to/existing_invivo_run
  --invitro_run_dir=/path/to/existing_invitro_run
  --seeds_csv=1,2,3
  --invivo_seeds_csv=1,2,3
  --invitro_seeds_csv=1,2,3
  --invivo_total_seeds=1
  --invitro_total_seeds=1
  --select_required_files=best_params.tsv
  --invivo_objective_columns=objective
  --invitro_objective_columns=objective_total,objective

Joint options:
  --joint_run_prefix=name
  --joint_seeds_csv=1,2,3
  --joint_total_seeds=1
  --parameter_table=/path/to/invitro_parameter_table.csv
  --fit_objects_dir=/path/to/fit_objects
  --flow_density_path=/path/to/g0g1_ploidy_density_grid.csv
  --invivo_best_seed_dir=/path/to/invivo/seed50
  --invitro_best_seed_dir=/path/to/invitro/seed350
  --joint_warmup_seed_label=invivo_seed50__invitro_seed350
  --joint_soft_coupling_sigma_default=0.65
  --joint_soft_coupling_welsch_c=0.4
  --joint_soft_coupling_parameters_table=/path/to/joint_soft_coupling_parameters_table.csv
  --joint_warmup_sigmaN=0.0304
  --joint_soft_coupling_delta_params=default|all|none|param1,param2
  --multi_warmup_top_n=10
  --multi_warmup_invivo_top_n=10  (0 disables in vivo source clustering; not both sides)
  --multi_warmup_invitro_top_n=10 (0 disables in vitro source clustering; not both sides)
  --multi_warmup_umap_seed=1
  --multi_warmup_invivo_k=auto
  --multi_warmup_invitro_anchor_ranks=1
  --multi_warmup_include_phase2=TRUE|FALSE
  --multi_warmup_phase2_invitro_anchor_ranks=auto

Local defaults run one seed per mode. Increase *_total_seeds or pass *_seeds_csv
explicitly when running a multi-seed local fit.
EOF
}

truthy() {
  case "${1:-FALSE}" in
    TRUE|true|True|1|yes|YES|y|Y|on|ON) return 0 ;;
    *) return 1 ;;
  esac
}

is_null_value() {
  local val
  val="$(echo "${1:-}" | tr '[:upper:]' '[:lower:]' | tr -d '[:space:]')"
  [[ -z "${val}" || "${val}" == "null" || "${val}" == "none" || "${val}" == "na" ]]
}

normalize_fitting_mode() {
  local val
  val="$(echo "${1:-}" | tr '[:upper:]' '[:lower:]')"
  val="${val// /}"
  val="${val//-/}"
  val="${val//_/}"
  case "${val}" in
    invivo) echo "invivo" ;;
    invitro) echo "invitro" ;;
    joint) echo "joint" ;;
    *) echo "" ;;
  esac
}

require_positive_int() {
  local name="$1"
  local value="$2"
  if ! [[ "${value}" =~ ^[0-9]+$ ]] || (( value <= 0 )); then
    echo "${name} must be a positive integer, got: ${value}" >&2
    exit 2
  fi
}

require_nonnegative_int() {
  local name="$1"
  local value="$2"
  if ! [[ "${value}" =~ ^[0-9]+$ ]]; then
    echo "${name} must be a non-negative integer, got: ${value}" >&2
    exit 2
  fi
}

print_command() {
  local label="$1"
  shift
  printf "%s:" "${label}"
  printf " %q" "$@"
  printf "\n"
}

sanitize_label() {
  local value="${1:-}"
  value="${value// /_}"
  value="$(printf "%s" "${value}" | tr -c 'A-Za-z0-9_.-' '_')"
  value="${value##_}"
  value="${value%%_}"
  if [[ -z "${value}" ]]; then
    value="seed"
  fi
  printf "%s" "${value}"
}

csv_from_total() {
  local total="$1"
  local out=""
  local i
  for ((i = 1; i <= total; i++)); do
    if [[ -n "${out}" ]]; then
      out+=","
    fi
    out+="${i}"
  done
  printf "%s" "${out}"
}

validate_seed_csv() {
  local name="$1"
  local csv="$2"
  local item
  IFS=',' read -r -a items <<< "${csv}"
  if [[ "${#items[@]}" -eq 0 ]]; then
    echo "${name} must not be empty." >&2
    exit 2
  fi
  for item in "${items[@]}"; do
    item="$(echo "${item}" | tr -d '[:space:]')"
    require_positive_int "${name}" "${item}"
  done
}

derive_joint_warmup_seed_label() {
  local invivo_label
  local invitro_label
  invivo_label="$(sanitize_label "$(basename "${INVIVO_BEST_SEED_DIR}")")"
  invitro_label="$(sanitize_label "$(basename "${INVITRO_BEST_SEED_DIR}")")"
  printf "invivo_%s__invitro_%s" "${invivo_label}" "${invitro_label}"
}

label_joint_run_prefix() {
  if truthy "${JOINT_WARMUP_ENABLE}" && ! is_null_value "${JOINT_WARMUP_SEED_LABEL}"; then
    if [[ "${JOINT_RUN_PREFIX}" != *"${JOINT_WARMUP_SEED_LABEL}"* ]]; then
      JOINT_RUN_PREFIX="${JOINT_RUN_PREFIX}__${JOINT_WARMUP_SEED_LABEL}"
    fi
  fi
}

parse_args() {
  for arg in "$@"; do
    case "${arg}" in
      --help|-h)
        usage
        exit 0
        ;;
      --fitting_mode=*) FITTING_MODE="${arg#*=}" ;;
      --joint_fitting_mode=*) JOINT_FITTING_MODE="${arg#*=}" ;;
      --project_root=*) PROJECT_ROOT="${arg#*=}" ;;
      --config_path=*|--config=*) CONFIG_PATH="${arg#*=}" ;;
      --out_root=*) OUT_ROOT="${arg#*=}" ;;
      --r_module=*) R_MODULE="${arg#*=}" ;;
      --dry_run=*) DRY_RUN="${arg#*=}" ;;
      --run_extra_results=*) RUN_EXTRA_RESULTS="${arg#*=}" ;;
      --force_extra_results=*) FORCE_EXTRA_RESULTS="${arg#*=}" ;;
      --append_run_prefix_timestamp=*) APPEND_RUN_PREFIX_TIMESTAMP="${arg#*=}" ;;
      --run_prefix=*) RUN_PREFIX="${arg#*=}" ;;
      --invivo_run_dir=*) INVIVO_RUN_DIR="${arg#*=}" ;;
      --invitro_run_dir=*) INVITRO_RUN_DIR="${arg#*=}" ;;
      --invivo_run_prefix=*) INVIVO_RUN_PREFIX="${arg#*=}" ;;
      --invitro_run_prefix=*) INVITRO_RUN_PREFIX="${arg#*=}" ;;
      --joint_run_prefix=*) JOINT_RUN_PREFIX="${arg#*=}" ;;
      --total_seeds=*) TOTAL_SEEDS="${arg#*=}" ;;
      --invivo_total_seeds=*) INVIVO_TOTAL_SEEDS="${arg#*=}" ;;
      --invitro_total_seeds=*) INVITRO_TOTAL_SEEDS="${arg#*=}" ;;
      --joint_total_seeds=*) JOINT_TOTAL_SEEDS="${arg#*=}" ;;
      --seeds_csv=*) SEEDS_CSV="${arg#*=}" ;;
      --invivo_seeds_csv=*) INVIVO_SEEDS_CSV="${arg#*=}" ;;
      --invitro_seeds_csv=*) INVITRO_SEEDS_CSV="${arg#*=}" ;;
      --joint_seeds_csv=*) JOINT_SEEDS_CSV="${arg#*=}" ;;
      --n_cores=*) N_CORES="${arg#*=}" ;;
      --invivo_n_cores=*) INVIVO_N_CORES="${arg#*=}" ;;
      --invitro_n_cores=*) INVITRO_N_CORES="${arg#*=}" ;;
      --joint_n_cores=*) JOINT_N_CORES="${arg#*=}" ;;
      --parameter_table=*) PARAMETER_TABLE="${arg#*=}" ;;
      --invitro_parameter_table=*) PARAMETER_TABLE="${arg#*=}" ;;
      --parameter_table_invitro=*) PARAMETER_TABLE="${arg#*=}" ;;
      --fit_objects_dir=*) FIT_OBJECTS_DIR="${arg#*=}" ;;
      --flow_density_path=*) FLOW_DENSITY_PATH="${arg#*=}" ;;
      --invivo_best_seed_dir=*|--joint_warmup_invivo_seed_dir=*|--joint_warmup_invivo_best_seed_dir=*) INVIVO_BEST_SEED_DIR="${arg#*=}"; USER_INVIVO_BEST_SEED_DIR="TRUE" ;;
      --invitro_best_seed_dir=*|--joint_warmup_invitro_seed_dir=*|--joint_warmup_invitro_best_seed_dir=*|--joint_warmup_vitro_seed_dir=*) INVITRO_BEST_SEED_DIR="${arg#*=}"; USER_INVITRO_BEST_SEED_DIR="TRUE" ;;
      --joint_warmup_enable=*) JOINT_WARMUP_ENABLE="${arg#*=}" ;;
      --joint_warmup_seed_label=*|--joint_seed_label=*|--seed_label=*) JOINT_WARMUP_SEED_LABEL="${arg#*=}" ;;
      --joint_warmup_sigmaN=*) JOINT_WARMUP_SIGMAN="${arg#*=}" ;;
      --joint_soft_coupling_sigma_default=*) JOINT_SOFT_COUPLING_SIGMA_DEFAULT="${arg#*=}" ;;
      --joint_soft_coupling_welsch_c=*) JOINT_SOFT_COUPLING_WELSCH_C="${arg#*=}" ;;
      --joint_soft_coupling_parameters_table=*|--joint_soft_coupling_parameters_table_path=*) JOINT_SOFT_COUPLING_PARAMETERS_TABLE="${arg#*=}" ;;
      --joint_soft_coupling_delta_params=*) JOINT_SOFT_COUPLING_DELTA_PARAMS="${arg#*=}" ;;
      --multi_warmup_top_n=*) MULTI_WARMUP_TOP_N="${arg#*=}" ;;
      --multi_warmup_invivo_top_n=*|--invivo_top_n=*) MULTI_WARMUP_INVIVO_TOP_N="${arg#*=}" ;;
      --multi_warmup_invitro_top_n=*|--invitro_top_n=*) MULTI_WARMUP_INVITRO_TOP_N="${arg#*=}" ;;
      --multi_warmup_umap_seed=*|--umap_seed=*) MULTI_WARMUP_UMAP_SEED="${arg#*=}" ;;
      --multi_warmup_invivo_k=*) MULTI_WARMUP_INVIVO_K="${arg#*=}" ;;
      --multi_warmup_invitro_k=*) MULTI_WARMUP_INVITRO_K="${arg#*=}" ;;
      --multi_warmup_invitro_anchor_ranks=*) MULTI_WARMUP_INVITRO_ANCHOR_RANKS="${arg#*=}" ;;
      --multi_warmup_include_phase2=*) MULTI_WARMUP_INCLUDE_PHASE2="${arg#*=}" ;;
      --multi_warmup_phase2_invitro_anchor_ranks=*) MULTI_WARMUP_PHASE2_INVITRO_ANCHOR_RANKS="${arg#*=}" ;;
      --select_required_files=*) SELECT_REQUIRED_FILES="${arg#*=}" ;;
      --invivo_objective_columns=*) INVIVO_OBJECTIVE_COLUMNS="${arg#*=}" ;;
      --invitro_objective_columns=*) INVITRO_OBJECTIVE_COLUMNS="${arg#*=}" ;;
      --itermax=*) ITERMAX="${arg#*=}" ;;
      --de_reltol=*) DE_RELTOL="${arg#*=}" ;;
      --de_steptol=*) DE_STEPTOL="${arg#*=}" ;;
      --np=*|--NP=*) NP="${arg#*=}" ;;
      --auto_viz=*) AUTO_VIZ="${arg#*=}" ;;
      *)
        echo "Unknown argument: ${arg}" >&2
        usage >&2
        exit 2
        ;;
    esac
  done
}

load_r_module() {
  if command -v ml >/dev/null 2>&1; then
    ml "${R_MODULE}"
  elif command -v module >/dev/null 2>&1; then
    module load "${R_MODULE}"
  fi
}

ensure_rscript() {
  if truthy "${DRY_RUN}"; then
    return
  fi
  load_r_module
  if ! command -v Rscript >/dev/null 2>&1; then
    echo "Rscript not found after loading ${R_MODULE}." >&2
    exit 1
  fi
}

run_or_print() {
  local label="$1"
  shift
  print_command "${label}" "$@"
  if ! truthy "${DRY_RUN}"; then
    "$@"
  fi
}

resolve_existing_dir() {
  local label="$1"
  local path="$2"
  if [[ "${path}" != /* && -d "${PROJECT_ROOT}/${path}" ]]; then
    path="${PROJECT_ROOT}/${path}"
  fi
  if [[ ! -d "${path}" ]]; then
    echo "Missing ${label}: ${path}" >&2
    exit 1
  fi
  (cd "${path}" && pwd)
}

first_line() {
  local path="$1"
  if [[ ! -f "${path}" ]]; then
    echo "Missing file: ${path}" >&2
    exit 1
  fi
  local line
  line="$(head -n 1 "${path}")"
  line="$(printf "%s" "${line}" | tr -d '\r' | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')"
  if [[ -z "${line}" ]]; then
    echo "Empty first line in ${path}" >&2
    exit 1
  fi
  printf "%s" "${line}"
}

maybe_flow_density_args() {
  FLOW_DENSITY_ARGS=()
  if [[ -n "${FLOW_DENSITY_PATH}" && -f "${FLOW_DENSITY_PATH}" ]]; then
    FLOW_DENSITY_ARGS=("--flow_density_path=${FLOW_DENSITY_PATH}")
  elif [[ -n "${FLOW_DENSITY_PATH}" ]]; then
    echo "Flow density file not found; continuing without flow-density input: ${FLOW_DENSITY_PATH}"
  fi
}

run_extra_results() {
  local label="$1"
  local run_dir="$2"
  local summary_path="${run_dir}/extra_results/seed_summary.tsv"
  local log_path="${run_dir}/extra_results_run.log"
  if ! truthy "${RUN_EXTRA_RESULTS}"; then
    echo "run_extra_results=FALSE; skipping ${label} extra_results."
    return
  fi
  if [[ -f "${summary_path}" ]] && ! truthy "${FORCE_EXTRA_RESULTS}"; then
    echo "extra_results already exists; skipping ${label}: ${summary_path}"
    return
  fi
  local cmd=(Rscript "${EXTRA_RESULTS_SCRIPT}" "--run_dir=${run_dir}")
  print_command "Run ${label} extra_results" "${cmd[@]}"
  if truthy "${DRY_RUN}"; then
    return
  fi
  "${cmd[@]}" > "${log_path}" 2>&1
  if [[ ! -f "${summary_path}" ]]; then
    echo "extra_results did not create ${summary_path}. Log: ${log_path}" >&2
    tail -40 "${log_path}" >&2 || true
    exit 1
  fi
}

select_best_seed() {
  local label="$1"
  local run_dir="$2"
  local objective_columns="$3"
  local log_path="${run_dir}/select_best_${label}.log"
  local cmd=(
    Rscript "${SELECT_BEST_SCRIPT}"
    "--run_dir=${run_dir}"
    "--objective_columns=${objective_columns}"
    "--required_files=${SELECT_REQUIRED_FILES}"
  )
  print_command "Select ${label} best seed" "${cmd[@]}"
  if truthy "${DRY_RUN}"; then
    return
  fi
  if ! "${cmd[@]}" > "${log_path}" 2>&1; then
    echo "Best-seed selection failed for ${label}. Log: ${log_path}" >&2
    tail -40 "${log_path}" >&2 || true
    exit 1
  fi
  cat "${log_path}"
}

prepare_joint_soft_coupling_table() {
  if is_null_value "${INVIVO_BEST_SEED_DIR}" && is_null_value "${INVITRO_BEST_SEED_DIR}"; then
    if truthy "${JOINT_WARMUP_ENABLE}"; then
      echo "joint_warmup_enable=TRUE requires --invivo_best_seed_dir and --invitro_best_seed_dir." >&2
      exit 2
    fi
    JOINT_WARMUP_ENABLE="FALSE"
    return
  fi
  if is_null_value "${INVIVO_BEST_SEED_DIR}" || is_null_value "${INVITRO_BEST_SEED_DIR}"; then
    echo "Both --invivo_best_seed_dir and --invitro_best_seed_dir are required for joint warm start." >&2
    exit 2
  fi
  if ! truthy "${DRY_RUN}"; then
    INVIVO_BEST_SEED_DIR="$(resolve_existing_dir "in vivo best seed directory" "${INVIVO_BEST_SEED_DIR}")"
    INVITRO_BEST_SEED_DIR="$(resolve_existing_dir "in vitro best seed directory" "${INVITRO_BEST_SEED_DIR}")"
  fi
  JOINT_WARMUP_ENABLE="TRUE"

  if is_null_value "${JOINT_WARMUP_SEED_LABEL}"; then
    JOINT_WARMUP_SEED_LABEL="$(derive_joint_warmup_seed_label)"
  else
    JOINT_WARMUP_SEED_LABEL="$(sanitize_label "${JOINT_WARMUP_SEED_LABEL}")"
  fi
  label_joint_run_prefix

  if [[ -z "${JOINT_SOFT_COUPLING_PARAMETERS_TABLE}" ]]; then
    JOINT_SOFT_COUPLING_PARAMETERS_TABLE="${PROJECT_ROOT}/oxygen/data/O2_supply_demand/joint_soft_coupling_parameters_table__${JOINT_WARMUP_SEED_LABEL}.csv"
  fi
  mkdir -p "$(dirname "${JOINT_SOFT_COUPLING_PARAMETERS_TABLE}")"
  JOINT_SOFT_COUPLING_PARAMETERS_TABLE="$(cd "$(dirname "${JOINT_SOFT_COUPLING_PARAMETERS_TABLE}")" && pwd)/$(basename "${JOINT_SOFT_COUPLING_PARAMETERS_TABLE}")"

  local cmd=(
    Rscript "${JOINT_WARM_START_SCRIPT}"
    "--invivo-seed-dir=${INVIVO_BEST_SEED_DIR}"
    "--invitro-seed-dir=${INVITRO_BEST_SEED_DIR}"
    "--seed-label=${JOINT_WARMUP_SEED_LABEL}"
    "--out=${JOINT_SOFT_COUPLING_PARAMETERS_TABLE}"
  )
  if [[ -n "${JOINT_SOFT_COUPLING_DELTA_PARAMS}" ]]; then
    cmd+=("--delta-params=${JOINT_SOFT_COUPLING_DELTA_PARAMS}")
  fi
  run_or_print "Generate joint soft-coupling table" "${cmd[@]}"
}

run_invivo_fit() {
  INVIVO_RUN_DIR="${OUT_ROOT}/${INVIVO_RUN_PREFIX}"
  mkdir -p "${INVIVO_RUN_DIR}"
  local cmd=(
    bash "${FIT_RUNNER_SCRIPT}"
    --fit_invivo
    --mode=run
    "--config=${CONFIG_PATH}"
    "--out_root=${OUT_ROOT}"
    "--run_prefix=${INVIVO_RUN_PREFIX}"
    "--append_run_prefix_timestamp=${APPEND_RUN_PREFIX_TIMESTAMP}"
    "--seeds_csv=${INVIVO_SEEDS_CSV}"
    "--n_cores=${INVIVO_N_CORES}"
    "--auto_viz=${AUTO_VIZ}"
  )
  run_or_print "Run in vivo fit" "${cmd[@]}"
}

run_invitro_fit() {
  INVITRO_RUN_DIR="${OUT_ROOT}/${INVITRO_RUN_PREFIX}"
  mkdir -p "${INVITRO_RUN_DIR}"
  maybe_flow_density_args
  local seed
  IFS=',' read -r -a invitro_seed_items <<< "${INVITRO_SEEDS_CSV}"
  for seed in "${invitro_seed_items[@]}"; do
    seed="$(echo "${seed}" | tr -d '[:space:]')"
    local seed_dir="${INVITRO_RUN_DIR}/seed${seed}"
    mkdir -p "${seed_dir}"
    local cmd=(
      bash "${FIT_RUNNER_SCRIPT}"
      --fit_invitro
      "--config=${CONFIG_PATH}"
      "--seed=${seed}"
      "--itermax=${ITERMAX}"
      "--de_reltol=${DE_RELTOL}"
      "--de_steptol=${DE_STEPTOL}"
      "--NP=${NP}"
      "--n_cores=${INVITRO_N_CORES}"
      "--out_dir=${seed_dir}"
      "--parameter_table=${PARAMETER_TABLE}"
      "--fit_objects_dir=${FIT_OBJECTS_DIR}"
      "--auto_viz=${AUTO_VIZ}"
      "${FLOW_DENSITY_ARGS[@]}"
    )
    run_or_print "Run in vitro seed ${seed}" "${cmd[@]}"
  done
}

build_joint_warmup_args() {
  JOINT_WARMUP_ARGS=()
  if [[ -n "${JOINT_SOFT_COUPLING_SIGMA_DEFAULT}" ]]; then
    JOINT_WARMUP_ARGS+=("--joint_soft_coupling_sigma_default=${JOINT_SOFT_COUPLING_SIGMA_DEFAULT}")
  fi
  if [[ -n "${JOINT_SOFT_COUPLING_WELSCH_C}" ]]; then
    JOINT_WARMUP_ARGS+=("--joint_soft_coupling_welsch_c=${JOINT_SOFT_COUPLING_WELSCH_C}")
  fi
  if truthy "${JOINT_WARMUP_ENABLE}"; then
    JOINT_WARMUP_ARGS+=(
      "--joint_warmup_enable=TRUE"
      "--joint_warmup_seed_label=${JOINT_WARMUP_SEED_LABEL}"
      "--joint_warmup_invivo_seed_dir=${INVIVO_BEST_SEED_DIR}"
      "--joint_warmup_invitro_seed_dir=${INVITRO_BEST_SEED_DIR}"
    )
    if [[ -n "${JOINT_WARMUP_SIGMAN}" ]]; then
      JOINT_WARMUP_ARGS+=("--joint_warmup_sigmaN=${JOINT_WARMUP_SIGMAN}")
    fi
  else
    JOINT_WARMUP_ARGS+=("--joint_warmup_enable=FALSE")
  fi
  if [[ -n "${JOINT_SOFT_COUPLING_PARAMETERS_TABLE}" ]]; then
    JOINT_WARMUP_ARGS+=("--joint_soft_coupling_parameters_table=${JOINT_SOFT_COUPLING_PARAMETERS_TABLE}")
  fi
}

run_joint_fit() {
  JOINT_RUN_DIR="${OUT_ROOT}/${JOINT_RUN_PREFIX}"
  mkdir -p "${JOINT_RUN_DIR}"
  maybe_flow_density_args
  build_joint_warmup_args
  local cmd=(
    bash "${JOINT_RUNNER_SCRIPT}"
    --mode=run
    "--config=${CONFIG_PATH}"
    "--out_root=${OUT_ROOT}"
    "--run_prefix=${JOINT_RUN_PREFIX}"
    "--append_run_prefix_timestamp=${APPEND_RUN_PREFIX_TIMESTAMP}"
    "--seeds_csv=${JOINT_SEEDS_CSV}"
    "--n_cores=${JOINT_N_CORES}"
    "--auto_viz=${AUTO_VIZ}"
    "--itermax=${ITERMAX}"
    "--de_reltol=${DE_RELTOL}"
    "--de_steptol=${DE_STEPTOL}"
    "--NP=${NP}"
    "--invitro_parameter_table=${PARAMETER_TABLE}"
    "--fit_objects_dir=${FIT_OBJECTS_DIR}"
    "${JOINT_WARMUP_ARGS[@]}"
    "${FLOW_DENSITY_ARGS[@]}"
  )
  run_or_print "Run joint fit" "${cmd[@]}"
}

run_best_seed_joint_pipeline() {
  if is_null_value "${INVIVO_BEST_SEED_DIR}"; then
    if is_null_value "${INVIVO_RUN_DIR}"; then
      run_invivo_fit
      run_extra_results "in vivo" "${INVIVO_RUN_DIR}"
      select_best_seed "invivo" "${INVIVO_RUN_DIR}" "${INVIVO_OBJECTIVE_COLUMNS}"
      if ! truthy "${DRY_RUN}"; then
        INVIVO_BEST_SEED_DIR="$(first_line "${INVIVO_RUN_DIR}/best_seed_from_summary.dir")"
      fi
    else
      INVIVO_RUN_DIR="$(resolve_existing_dir "in vivo run directory" "${INVIVO_RUN_DIR}")"
      echo "Skipping in vivo fitting; using existing run directory: ${INVIVO_RUN_DIR}"
      run_extra_results "in vivo" "${INVIVO_RUN_DIR}"
      select_best_seed "invivo" "${INVIVO_RUN_DIR}" "${INVIVO_OBJECTIVE_COLUMNS}"
      if ! truthy "${DRY_RUN}"; then
        INVIVO_BEST_SEED_DIR="$(first_line "${INVIVO_RUN_DIR}/best_seed_from_summary.dir")"
      fi
    fi
  else
    INVIVO_BEST_SEED_DIR="$(resolve_existing_dir "in vivo best seed directory" "${INVIVO_BEST_SEED_DIR}")"
    if is_null_value "${INVIVO_RUN_DIR}"; then
      INVIVO_RUN_DIR="$(cd "$(dirname "${INVIVO_BEST_SEED_DIR}")" && pwd)"
    else
      INVIVO_RUN_DIR="$(resolve_existing_dir "in vivo run directory" "${INVIVO_RUN_DIR}")"
    fi
    echo "Skipping in vivo fitting and best-seed selection; using: ${INVIVO_BEST_SEED_DIR}"
  fi

  if is_null_value "${INVITRO_BEST_SEED_DIR}"; then
    if is_null_value "${INVITRO_RUN_DIR}"; then
      run_invitro_fit
      run_extra_results "in vitro" "${INVITRO_RUN_DIR}"
      select_best_seed "invitro" "${INVITRO_RUN_DIR}" "${INVITRO_OBJECTIVE_COLUMNS}"
      if ! truthy "${DRY_RUN}"; then
        INVITRO_BEST_SEED_DIR="$(first_line "${INVITRO_RUN_DIR}/best_seed_from_summary.dir")"
      fi
    else
      INVITRO_RUN_DIR="$(resolve_existing_dir "in vitro run directory" "${INVITRO_RUN_DIR}")"
      echo "Skipping in vitro fitting; using existing run directory: ${INVITRO_RUN_DIR}"
      run_extra_results "in vitro" "${INVITRO_RUN_DIR}"
      select_best_seed "invitro" "${INVITRO_RUN_DIR}" "${INVITRO_OBJECTIVE_COLUMNS}"
      if ! truthy "${DRY_RUN}"; then
        INVITRO_BEST_SEED_DIR="$(first_line "${INVITRO_RUN_DIR}/best_seed_from_summary.dir")"
      fi
    fi
  else
    INVITRO_BEST_SEED_DIR="$(resolve_existing_dir "in vitro best seed directory" "${INVITRO_BEST_SEED_DIR}")"
    if is_null_value "${INVITRO_RUN_DIR}"; then
      INVITRO_RUN_DIR="$(cd "$(dirname "${INVITRO_BEST_SEED_DIR}")" && pwd)"
    else
      INVITRO_RUN_DIR="$(resolve_existing_dir "in vitro run directory" "${INVITRO_RUN_DIR}")"
    fi
    echo "Skipping in vitro fitting and best-seed selection; using: ${INVITRO_BEST_SEED_DIR}"
  fi

  if truthy "${DRY_RUN}"; then
    echo "DRY_RUN=TRUE; best-seed dirs are not materialized, so joint warm-start table generation is shown with placeholder paths."
    INVIVO_BEST_SEED_DIR="${INVIVO_BEST_SEED_DIR:-${INVIVO_RUN_DIR}/<selected_seed>}"
    INVITRO_BEST_SEED_DIR="${INVITRO_BEST_SEED_DIR:-${INVITRO_RUN_DIR}/<selected_seed>}"
  fi

  prepare_joint_soft_coupling_table
  run_joint_fit
  run_extra_results "joint" "${OUT_ROOT}/${JOINT_RUN_PREFIX}"
}

run_multi_warmup_pipeline() {
  if (( MULTI_WARMUP_INVIVO_TOP_N > 0 )); then
    if ! is_null_value "${INVIVO_RUN_DIR}"; then
      INVIVO_RUN_DIR="$(resolve_existing_dir "in vivo run directory" "${INVIVO_RUN_DIR}")"
      echo "MULTI_WARMUP using existing in vivo run directory: ${INVIVO_RUN_DIR}"
    else
      echo "MULTI_WARMUP no in vivo run directory supplied; running in vivo source fit first."
      run_invivo_fit
      run_extra_results "in vivo" "${INVIVO_RUN_DIR}"
    fi
  else
    INVIVO_RUN_DIR=""
    echo "MULTI_WARMUP invivo_top_n=0; skipping in vivo source fit."
  fi

  if (( MULTI_WARMUP_INVITRO_TOP_N > 0 )); then
    if ! is_null_value "${INVITRO_RUN_DIR}"; then
      INVITRO_RUN_DIR="$(resolve_existing_dir "in vitro run directory" "${INVITRO_RUN_DIR}")"
      echo "MULTI_WARMUP using existing in vitro run directory: ${INVITRO_RUN_DIR}"
    else
      echo "MULTI_WARMUP no in vitro run directory supplied; running in vitro source fit first."
      run_invitro_fit
      run_extra_results "in vitro" "${INVITRO_RUN_DIR}"
    fi
  else
    INVITRO_RUN_DIR=""
    echo "MULTI_WARMUP invitro_top_n=0; skipping in vitro source fit."
  fi

  local cmd=(
    bash "${MULTI_WARMUP_RUNNER_SCRIPT}"
    "--project_root=${PROJECT_ROOT}"
    "--config_path=${CONFIG_PATH}"
    "--out_root=${OUT_ROOT}"
    "--multi_warmup_prefix=${JOINT_RUN_PREFIX}"
    "--invivo_run_dir=${INVIVO_RUN_DIR}"
    "--invitro_run_dir=${INVITRO_RUN_DIR}"
    "--joint_seeds_csv=${JOINT_SEEDS_CSV}"
    "--joint_n_cores=${JOINT_N_CORES}"
    "--r_module=${R_MODULE}"
    "--dry_run=${DRY_RUN}"
    "--run_extra_results=${RUN_EXTRA_RESULTS}"
    "--force_extra_results=${FORCE_EXTRA_RESULTS}"
    "--parameter_table=${PARAMETER_TABLE}"
    "--fit_objects_dir=${FIT_OBJECTS_DIR}"
    "--flow_density_path=${FLOW_DENSITY_PATH}"
    "--itermax=${ITERMAX}"
    "--de_reltol=${DE_RELTOL}"
    "--de_steptol=${DE_STEPTOL}"
    "--NP=${NP}"
    "--auto_viz=${AUTO_VIZ}"
    "--joint_soft_coupling_delta_params=${JOINT_SOFT_COUPLING_DELTA_PARAMS}"
    "--multi_warmup_top_n=${MULTI_WARMUP_TOP_N}"
    "--multi_warmup_invivo_top_n=${MULTI_WARMUP_INVIVO_TOP_N}"
    "--multi_warmup_invitro_top_n=${MULTI_WARMUP_INVITRO_TOP_N}"
    "--multi_warmup_umap_seed=${MULTI_WARMUP_UMAP_SEED}"
    "--multi_warmup_invivo_k=${MULTI_WARMUP_INVIVO_K}"
    "--multi_warmup_invitro_k=${MULTI_WARMUP_INVITRO_K}"
    "--multi_warmup_invitro_anchor_ranks=${MULTI_WARMUP_INVITRO_ANCHOR_RANKS}"
    "--multi_warmup_include_phase2=${MULTI_WARMUP_INCLUDE_PHASE2}"
    "--multi_warmup_phase2_invitro_anchor_ranks=${MULTI_WARMUP_PHASE2_INVITRO_ANCHOR_RANKS}"
  )
  if [[ -n "${JOINT_WARMUP_SIGMAN}" ]]; then
    cmd+=("--joint_warmup_sigmaN=${JOINT_WARMUP_SIGMAN}")
  fi
  if [[ -n "${JOINT_SOFT_COUPLING_SIGMA_DEFAULT}" ]]; then
    cmd+=("--joint_soft_coupling_sigma_default=${JOINT_SOFT_COUPLING_SIGMA_DEFAULT}")
  fi
  if [[ -n "${JOINT_SOFT_COUPLING_WELSCH_C}" ]]; then
    cmd+=("--joint_soft_coupling_welsch_c=${JOINT_SOFT_COUPLING_WELSCH_C}")
  fi
  run_or_print "Run multi-warm-up joint sweep" "${cmd[@]}"
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKFLOW_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
DEFAULT_PROJECT_ROOT="$(cd "${WORKFLOW_ROOT}/../../.." && pwd)"
DEFAULT_R_MODULE="R/4.4"
DEFAULT_INVIVO_RUN_PREFIX="fit_invivo_O2_buffering_local"
DEFAULT_INVITRO_RUN_PREFIX="fit_invitro_O2_buffering_local"
DEFAULT_JOINT_RUN_PREFIX="fit_joint_O2_buffering_local"
DEFAULT_INVIVO_BEST_SEED_REL="oxygen/results/fit_invivo_O2_buffering_500seed/seed50"
DEFAULT_INVITRO_BEST_SEED_REL="oxygen/results/fit_invitro_O2_buffering_500seed/seed350"
DEFAULT_TOTAL_SEEDS="1"
DEFAULT_N_CORES="1"
DEFAULT_AUTO_VIZ="TRUE"
DEFAULT_RUN_EXTRA_RESULTS="TRUE"
DEFAULT_FORCE_EXTRA_RESULTS="FALSE"
DEFAULT_DRY_RUN="FALSE"
DEFAULT_APPEND_RUN_PREFIX_TIMESTAMP="FALSE"
DEFAULT_ITERMAX="500"
DEFAULT_DE_RELTOL="1e-4"
DEFAULT_DE_STEPTOL="25"
DEFAULT_NP="80"
DEFAULT_SELECT_REQUIRED_FILES="best_params.tsv"
DEFAULT_INVIVO_OBJECTIVE_COLUMNS="objective"
DEFAULT_INVITRO_OBJECTIVE_COLUMNS="objective_total,objective"
DEFAULT_JOINT_WARMUP_ENABLE="TRUE"
DEFAULT_JOINT_WARMUP_SEED_LABEL=""
DEFAULT_JOINT_WARMUP_SIGMAN=""
DEFAULT_JOINT_SOFT_COUPLING_SIGMA_DEFAULT=""
DEFAULT_JOINT_SOFT_COUPLING_WELSCH_C=""
DEFAULT_JOINT_SOFT_COUPLING_DELTA_PARAMS="default"
DEFAULT_MULTI_WARMUP_TOP_N="10"
DEFAULT_MULTI_WARMUP_INVIVO_TOP_N=""
DEFAULT_MULTI_WARMUP_INVITRO_TOP_N=""
DEFAULT_MULTI_WARMUP_UMAP_SEED="1"
DEFAULT_MULTI_WARMUP_INVIVO_K="auto"
DEFAULT_MULTI_WARMUP_INVITRO_K="auto"
DEFAULT_MULTI_WARMUP_INVITRO_ANCHOR_RANKS="1"
DEFAULT_MULTI_WARMUP_INCLUDE_PHASE2="FALSE"
DEFAULT_MULTI_WARMUP_PHASE2_INVITRO_ANCHOR_RANKS="auto"

FITTING_MODE="${FITTING_MODE:-}"
JOINT_FITTING_MODE="${JOINT_FITTING_MODE:-}"
PROJECT_ROOT="${PROJECT_ROOT:-}"
R_MODULE="${R_MODULE:-}"
CONFIG_PATH="${CONFIG_PATH:-}"
OUT_ROOT="${OUT_ROOT:-}"
RUN_PREFIX="${RUN_PREFIX:-}"
INVIVO_RUN_PREFIX="${INVIVO_RUN_PREFIX:-}"
INVITRO_RUN_PREFIX="${INVITRO_RUN_PREFIX:-}"
JOINT_RUN_PREFIX="${JOINT_RUN_PREFIX:-}"
TOTAL_SEEDS="${TOTAL_SEEDS:-}"
INVIVO_TOTAL_SEEDS="${INVIVO_TOTAL_SEEDS:-}"
INVITRO_TOTAL_SEEDS="${INVITRO_TOTAL_SEEDS:-}"
JOINT_TOTAL_SEEDS="${JOINT_TOTAL_SEEDS:-}"
SEEDS_CSV="${SEEDS_CSV:-}"
INVIVO_SEEDS_CSV="${INVIVO_SEEDS_CSV:-}"
INVITRO_SEEDS_CSV="${INVITRO_SEEDS_CSV:-}"
JOINT_SEEDS_CSV="${JOINT_SEEDS_CSV:-}"
N_CORES="${N_CORES:-}"
INVIVO_N_CORES="${INVIVO_N_CORES:-}"
INVITRO_N_CORES="${INVITRO_N_CORES:-}"
JOINT_N_CORES="${JOINT_N_CORES:-}"
PARAMETER_TABLE="${PARAMETER_TABLE:-}"
FIT_OBJECTS_DIR="${FIT_OBJECTS_DIR:-}"
FLOW_DENSITY_PATH="${FLOW_DENSITY_PATH:-}"
ITERMAX="${ITERMAX:-}"
DE_RELTOL="${DE_RELTOL:-}"
DE_STEPTOL="${DE_STEPTOL:-}"
NP="${NP:-}"
AUTO_VIZ="${AUTO_VIZ:-}"
RUN_EXTRA_RESULTS="${RUN_EXTRA_RESULTS:-}"
FORCE_EXTRA_RESULTS="${FORCE_EXTRA_RESULTS:-}"
DRY_RUN="${DRY_RUN:-}"
APPEND_RUN_PREFIX_TIMESTAMP="${APPEND_RUN_PREFIX_TIMESTAMP:-}"
INVIVO_RUN_DIR="${INVIVO_RUN_DIR:-}"
INVITRO_RUN_DIR="${INVITRO_RUN_DIR:-}"
INVIVO_BEST_SEED_DIR="${INVIVO_BEST_SEED_DIR:-}"
INVITRO_BEST_SEED_DIR="${INVITRO_BEST_SEED_DIR:-}"
USER_INVIVO_BEST_SEED_DIR="${USER_INVIVO_BEST_SEED_DIR:-FALSE}"
USER_INVITRO_BEST_SEED_DIR="${USER_INVITRO_BEST_SEED_DIR:-FALSE}"
JOINT_WARMUP_ENABLE="${JOINT_WARMUP_ENABLE:-}"
JOINT_WARMUP_SEED_LABEL="${JOINT_WARMUP_SEED_LABEL:-}"
JOINT_WARMUP_SIGMAN="${JOINT_WARMUP_SIGMAN:-}"
JOINT_SOFT_COUPLING_SIGMA_DEFAULT="${JOINT_SOFT_COUPLING_SIGMA_DEFAULT:-}"
JOINT_SOFT_COUPLING_WELSCH_C="${JOINT_SOFT_COUPLING_WELSCH_C:-}"
JOINT_SOFT_COUPLING_PARAMETERS_TABLE="${JOINT_SOFT_COUPLING_PARAMETERS_TABLE:-}"
JOINT_SOFT_COUPLING_DELTA_PARAMS="${JOINT_SOFT_COUPLING_DELTA_PARAMS:-}"
MULTI_WARMUP_TOP_N="${MULTI_WARMUP_TOP_N:-}"
MULTI_WARMUP_INVIVO_TOP_N="${MULTI_WARMUP_INVIVO_TOP_N:-}"
MULTI_WARMUP_INVITRO_TOP_N="${MULTI_WARMUP_INVITRO_TOP_N:-}"
MULTI_WARMUP_UMAP_SEED="${MULTI_WARMUP_UMAP_SEED:-}"
MULTI_WARMUP_INVIVO_K="${MULTI_WARMUP_INVIVO_K:-}"
MULTI_WARMUP_INVITRO_K="${MULTI_WARMUP_INVITRO_K:-}"
MULTI_WARMUP_INVITRO_ANCHOR_RANKS="${MULTI_WARMUP_INVITRO_ANCHOR_RANKS:-}"
MULTI_WARMUP_INCLUDE_PHASE2="${MULTI_WARMUP_INCLUDE_PHASE2:-}"
MULTI_WARMUP_PHASE2_INVITRO_ANCHOR_RANKS="${MULTI_WARMUP_PHASE2_INVITRO_ANCHOR_RANKS:-}"
SELECT_REQUIRED_FILES="${SELECT_REQUIRED_FILES:-}"
INVIVO_OBJECTIVE_COLUMNS="${INVIVO_OBJECTIVE_COLUMNS:-}"
INVITRO_OBJECTIVE_COLUMNS="${INVITRO_OBJECTIVE_COLUMNS:-}"

parse_args "$@"

PROJECT_ROOT="${PROJECT_ROOT:-${DEFAULT_PROJECT_ROOT}}"
R_MODULE="${R_MODULE:-${DEFAULT_R_MODULE}}"
INVIVO_RUN_PREFIX="${INVIVO_RUN_PREFIX:-${RUN_PREFIX:-${DEFAULT_INVIVO_RUN_PREFIX}}}"
INVITRO_RUN_PREFIX="${INVITRO_RUN_PREFIX:-${RUN_PREFIX:-${DEFAULT_INVITRO_RUN_PREFIX}}}"
JOINT_RUN_PREFIX="${JOINT_RUN_PREFIX:-${RUN_PREFIX:-${DEFAULT_JOINT_RUN_PREFIX}}}"
TOTAL_SEEDS="${TOTAL_SEEDS:-${DEFAULT_TOTAL_SEEDS}}"
INVIVO_TOTAL_SEEDS="${INVIVO_TOTAL_SEEDS:-${TOTAL_SEEDS}}"
INVITRO_TOTAL_SEEDS="${INVITRO_TOTAL_SEEDS:-${TOTAL_SEEDS}}"
JOINT_TOTAL_SEEDS="${JOINT_TOTAL_SEEDS:-${TOTAL_SEEDS}}"
INVIVO_SEEDS_CSV="${INVIVO_SEEDS_CSV:-${SEEDS_CSV:-$(csv_from_total "${INVIVO_TOTAL_SEEDS}")}}"
INVITRO_SEEDS_CSV="${INVITRO_SEEDS_CSV:-${SEEDS_CSV:-$(csv_from_total "${INVITRO_TOTAL_SEEDS}")}}"
JOINT_SEEDS_CSV="${JOINT_SEEDS_CSV:-${SEEDS_CSV:-$(csv_from_total "${JOINT_TOTAL_SEEDS}")}}"
N_CORES="${N_CORES:-${DEFAULT_N_CORES}}"
INVIVO_N_CORES="${INVIVO_N_CORES:-${N_CORES}}"
INVITRO_N_CORES="${INVITRO_N_CORES:-${N_CORES}}"
JOINT_N_CORES="${JOINT_N_CORES:-${N_CORES}}"
ITERMAX="${ITERMAX:-${DEFAULT_ITERMAX}}"
DE_RELTOL="${DE_RELTOL:-${DEFAULT_DE_RELTOL}}"
DE_STEPTOL="${DE_STEPTOL:-${DEFAULT_DE_STEPTOL}}"
NP="${NP:-${DEFAULT_NP}}"
AUTO_VIZ="${AUTO_VIZ:-${DEFAULT_AUTO_VIZ}}"
RUN_EXTRA_RESULTS="${RUN_EXTRA_RESULTS:-${DEFAULT_RUN_EXTRA_RESULTS}}"
FORCE_EXTRA_RESULTS="${FORCE_EXTRA_RESULTS:-${DEFAULT_FORCE_EXTRA_RESULTS}}"
DRY_RUN="${DRY_RUN:-${DEFAULT_DRY_RUN}}"
APPEND_RUN_PREFIX_TIMESTAMP="${APPEND_RUN_PREFIX_TIMESTAMP:-${DEFAULT_APPEND_RUN_PREFIX_TIMESTAMP}}"
JOINT_WARMUP_ENABLE="${JOINT_WARMUP_ENABLE:-${DEFAULT_JOINT_WARMUP_ENABLE}}"
if truthy "${JOINT_WARMUP_ENABLE}"; then
  INVIVO_BEST_SEED_DIR="${INVIVO_BEST_SEED_DIR:-${PROJECT_ROOT}/${DEFAULT_INVIVO_BEST_SEED_REL}}"
  INVITRO_BEST_SEED_DIR="${INVITRO_BEST_SEED_DIR:-${PROJECT_ROOT}/${DEFAULT_INVITRO_BEST_SEED_REL}}"
else
  INVIVO_BEST_SEED_DIR="${INVIVO_BEST_SEED_DIR:-}"
  INVITRO_BEST_SEED_DIR="${INVITRO_BEST_SEED_DIR:-}"
fi
JOINT_WARMUP_SEED_LABEL="${JOINT_WARMUP_SEED_LABEL:-${DEFAULT_JOINT_WARMUP_SEED_LABEL}}"
JOINT_WARMUP_SIGMAN="${JOINT_WARMUP_SIGMAN:-${DEFAULT_JOINT_WARMUP_SIGMAN}}"
JOINT_SOFT_COUPLING_SIGMA_DEFAULT="${JOINT_SOFT_COUPLING_SIGMA_DEFAULT:-${DEFAULT_JOINT_SOFT_COUPLING_SIGMA_DEFAULT}}"
JOINT_SOFT_COUPLING_WELSCH_C="${JOINT_SOFT_COUPLING_WELSCH_C:-${DEFAULT_JOINT_SOFT_COUPLING_WELSCH_C}}"
JOINT_SOFT_COUPLING_DELTA_PARAMS="${JOINT_SOFT_COUPLING_DELTA_PARAMS:-${DEFAULT_JOINT_SOFT_COUPLING_DELTA_PARAMS}}"
MULTI_WARMUP_TOP_N="${MULTI_WARMUP_TOP_N:-${DEFAULT_MULTI_WARMUP_TOP_N}}"
MULTI_WARMUP_INVIVO_TOP_N="${MULTI_WARMUP_INVIVO_TOP_N:-${DEFAULT_MULTI_WARMUP_INVIVO_TOP_N:-${MULTI_WARMUP_TOP_N}}}"
MULTI_WARMUP_INVITRO_TOP_N="${MULTI_WARMUP_INVITRO_TOP_N:-${DEFAULT_MULTI_WARMUP_INVITRO_TOP_N:-${MULTI_WARMUP_TOP_N}}}"
MULTI_WARMUP_UMAP_SEED="${MULTI_WARMUP_UMAP_SEED:-${DEFAULT_MULTI_WARMUP_UMAP_SEED}}"
MULTI_WARMUP_INVIVO_K="${MULTI_WARMUP_INVIVO_K:-${DEFAULT_MULTI_WARMUP_INVIVO_K}}"
MULTI_WARMUP_INVITRO_K="${MULTI_WARMUP_INVITRO_K:-${DEFAULT_MULTI_WARMUP_INVITRO_K}}"
MULTI_WARMUP_INVITRO_ANCHOR_RANKS="${MULTI_WARMUP_INVITRO_ANCHOR_RANKS:-${DEFAULT_MULTI_WARMUP_INVITRO_ANCHOR_RANKS}}"
MULTI_WARMUP_INCLUDE_PHASE2="${MULTI_WARMUP_INCLUDE_PHASE2:-${DEFAULT_MULTI_WARMUP_INCLUDE_PHASE2}}"
MULTI_WARMUP_PHASE2_INVITRO_ANCHOR_RANKS="${MULTI_WARMUP_PHASE2_INVITRO_ANCHOR_RANKS:-${DEFAULT_MULTI_WARMUP_PHASE2_INVITRO_ANCHOR_RANKS}}"
SELECT_REQUIRED_FILES="${SELECT_REQUIRED_FILES:-${DEFAULT_SELECT_REQUIRED_FILES}}"
INVIVO_OBJECTIVE_COLUMNS="${INVIVO_OBJECTIVE_COLUMNS:-${DEFAULT_INVIVO_OBJECTIVE_COLUMNS}}"
INVITRO_OBJECTIVE_COLUMNS="${INVITRO_OBJECTIVE_COLUMNS:-${DEFAULT_INVITRO_OBJECTIVE_COLUMNS}}"

FITTING_MODE="$(normalize_fitting_mode "${FITTING_MODE}")"
if [[ -z "${FITTING_MODE}" ]]; then
  echo "--fitting_mode must be one of invivo, invitro, or joint." >&2
  usage >&2
  exit 2
fi

if [[ "${FITTING_MODE}" != "joint" ]]; then
  JOINT_FITTING_MODE="OFF"
elif [[ -z "${JOINT_FITTING_MODE}" ]]; then
  JOINT_FITTING_MODE="DIRECT"
fi
JOINT_FITTING_MODE="$(echo "${JOINT_FITTING_MODE}" | tr '[:lower:]' '[:upper:]')"
case "${JOINT_FITTING_MODE}" in
  SINGLE) JOINT_FITTING_MODE="DIRECT" ;;
esac
case "${JOINT_FITTING_MODE}" in
  OFF|JOINT|DIRECT|MULTI_WARMUP) ;;
  *)
    echo "--joint_fitting_mode must be OFF, JOINT, DIRECT, or MULTI_WARMUP. SINGLE is accepted as a legacy alias for DIRECT." >&2
    exit 2
    ;;
esac

if [[ "${JOINT_FITTING_MODE}" == "MULTI_WARMUP" ]]; then
  if truthy "${USER_INVIVO_BEST_SEED_DIR}" || truthy "${USER_INVITRO_BEST_SEED_DIR}"; then
    echo "MULTI_WARMUP mode does not accept user-specified best seed directories; provide --invivo_run_dir and --invitro_run_dir instead." >&2
    exit 2
  fi
  require_nonnegative_int MULTI_WARMUP_INVIVO_TOP_N "${MULTI_WARMUP_INVIVO_TOP_N}"
  require_nonnegative_int MULTI_WARMUP_INVITRO_TOP_N "${MULTI_WARMUP_INVITRO_TOP_N}"
  if (( MULTI_WARMUP_INVIVO_TOP_N == 0 && MULTI_WARMUP_INVITRO_TOP_N == 0 )); then
    echo "At least one of MULTI_WARMUP_INVIVO_TOP_N or MULTI_WARMUP_INVITRO_TOP_N must be greater than 0." >&2
    exit 2
  fi
  INVIVO_BEST_SEED_DIR=""
  INVITRO_BEST_SEED_DIR=""
  JOINT_WARMUP_ENABLE="FALSE"
fi

require_positive_int INVIVO_TOTAL_SEEDS "${INVIVO_TOTAL_SEEDS}"
require_positive_int INVITRO_TOTAL_SEEDS "${INVITRO_TOTAL_SEEDS}"
require_positive_int JOINT_TOTAL_SEEDS "${JOINT_TOTAL_SEEDS}"
require_positive_int INVIVO_N_CORES "${INVIVO_N_CORES}"
require_positive_int INVITRO_N_CORES "${INVITRO_N_CORES}"
require_positive_int JOINT_N_CORES "${JOINT_N_CORES}"
require_positive_int ITERMAX "${ITERMAX}"
require_positive_int DE_STEPTOL "${DE_STEPTOL}"
require_positive_int NP "${NP}"
validate_seed_csv INVIVO_SEEDS_CSV "${INVIVO_SEEDS_CSV}"
validate_seed_csv INVITRO_SEEDS_CSV "${INVITRO_SEEDS_CSV}"
validate_seed_csv JOINT_SEEDS_CSV "${JOINT_SEEDS_CSV}"

if truthy "${APPEND_RUN_PREFIX_TIMESTAMP}" && truthy "${RUN_EXTRA_RESULTS}"; then
  echo "--append_run_prefix_timestamp=TRUE is not supported with --run_extra_results=TRUE because the run directory is not deterministic." >&2
  echo "Use --append_run_prefix_timestamp=FALSE or --run_extra_results=FALSE." >&2
  exit 2
fi

PROJECT_ROOT="$(cd "${PROJECT_ROOT}" && pwd)"
if [[ -z "${CONFIG_PATH}" ]]; then
  CONFIG_PATH="${PROJECT_ROOT}/oxygen/config/O2_supply_demand.yaml"
fi
if [[ -z "${OUT_ROOT}" ]]; then
  OUT_ROOT="${PROJECT_ROOT}/oxygen/results"
fi
CONFIG_PATH="$(cd "$(dirname "${CONFIG_PATH}")" && pwd)/$(basename "${CONFIG_PATH}")"
mkdir -p "${OUT_ROOT}"
OUT_ROOT="$(cd "${OUT_ROOT}" && pwd)"

FIT_RUNNER_SCRIPT="${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/runner/run_fit_model_O2_supply_demand_MAP.sh"
JOINT_RUNNER_SCRIPT="${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/runner/run_fit_joint_model_O2_supply_demand_MAP.sh"
MULTI_WARMUP_RUNNER_SCRIPT="${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/runner/run_multi_warmup_joint.sh"
EXTRA_RESULTS_SCRIPT="${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/analysis/fit_results/extra_results.R"
SELECT_BEST_SCRIPT="${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/analysis/fit_results/select_best_seed_from_summary.R"
JOINT_WARM_START_SCRIPT="${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/analysis/warm_start/make_joint_soft_coupling_parameters_table.R"

if [[ -z "${PARAMETER_TABLE}" ]]; then
  PARAMETER_TABLE="${PROJECT_ROOT}/oxygen/data/O2_supply_demand/parameter_table_invitro_buffering.csv"
fi
if [[ -z "${FIT_OBJECTS_DIR}" ]]; then
  FIT_OBJECTS_DIR="${PROJECT_ROOT}/oxygen/ploidyOxygen/data/fit_objects"
fi
if [[ -z "${FLOW_DENSITY_PATH}" ]]; then
  FLOW_DENSITY_PATH="${PROJECT_ROOT}/oxygen/data/g0g1_ploidy_density_grid.csv"
fi
PARAMETER_TABLE="$(cd "$(dirname "${PARAMETER_TABLE}")" && pwd)/$(basename "${PARAMETER_TABLE}")"
FIT_OBJECTS_DIR="$(cd "${FIT_OBJECTS_DIR}" && pwd)"
FLOW_DENSITY_PATH="$(cd "$(dirname "${FLOW_DENSITY_PATH}")" && pwd)/$(basename "${FLOW_DENSITY_PATH}")"

for path in "${CONFIG_PATH}" "${FIT_RUNNER_SCRIPT}" "${JOINT_RUNNER_SCRIPT}" "${MULTI_WARMUP_RUNNER_SCRIPT}" \
            "${EXTRA_RESULTS_SCRIPT}" "${SELECT_BEST_SCRIPT}" "${JOINT_WARM_START_SCRIPT}" \
            "${PARAMETER_TABLE}"; do
  if [[ ! -f "${path}" ]]; then
    echo "Missing required file: ${path}" >&2
    exit 1
  fi
done
if [[ ! -d "${FIT_OBJECTS_DIR}" ]]; then
  echo "Missing fit_objects_dir: ${FIT_OBJECTS_DIR}" >&2
  exit 1
fi

ensure_rscript

echo "O2 local runner"
echo "  fitting_mode: ${FITTING_MODE}"
echo "  joint_fitting_mode: ${JOINT_FITTING_MODE}"
echo "  project_root: ${PROJECT_ROOT}"
echo "  out_root: ${OUT_ROOT}"
echo "  config_path: ${CONFIG_PATH}"
echo "  invivo seeds: ${INVIVO_SEEDS_CSV}; n_cores=${INVIVO_N_CORES}; run_prefix=${INVIVO_RUN_PREFIX}"
echo "  invitro seeds: ${INVITRO_SEEDS_CSV}; n_cores=${INVITRO_N_CORES}; run_prefix=${INVITRO_RUN_PREFIX}"
echo "  joint seeds: ${JOINT_SEEDS_CSV}; n_cores=${JOINT_N_CORES}; run_prefix=${JOINT_RUN_PREFIX}"
echo "  run_extra_results: ${RUN_EXTRA_RESULTS}"
echo "  dry_run: ${DRY_RUN}"

case "${FITTING_MODE}" in
  invivo)
    run_invivo_fit
    run_extra_results "in vivo" "${INVIVO_RUN_DIR}"
    ;;
  invitro)
    run_invitro_fit
    run_extra_results "in vitro" "${INVITRO_RUN_DIR}"
    ;;
  joint)
    case "${JOINT_FITTING_MODE}" in
      OFF)
        echo "joint_fitting_mode=OFF; no joint fitting run."
        ;;
      JOINT)
        echo "joint_fitting_mode=JOINT using local best-seed selection pipeline."
        run_best_seed_joint_pipeline
        ;;
      DIRECT)
        if ! is_null_value "${INVIVO_RUN_DIR}" || ! is_null_value "${INVITRO_RUN_DIR}"; then
          echo "Ignoring invivo_run_dir/invitro_run_dir for DIRECT joint mode."
        fi
        prepare_joint_soft_coupling_table
        run_joint_fit
        run_extra_results "joint" "${OUT_ROOT}/${JOINT_RUN_PREFIX}"
        ;;
      MULTI_WARMUP)
        echo "joint_fitting_mode=MULTI_WARMUP using source-run ratio clustering."
        run_multi_warmup_pipeline
        ;;
    esac
    ;;
esac
