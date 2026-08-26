#!/usr/bin/env bash
# Unified local O2 runner for in vivo, in vitro, and joint fitting.

set -euo pipefail

O2SD_SHELL_UTILS="$(cd "$(dirname "${BASH_SOURCE[0]}")/../util" && pwd)/o2_supply_demand_map_shell_utils.sh"
# shellcheck source=../util/o2_supply_demand_map_shell_utils.sh
source "${O2SD_SHELL_UTILS}"

usage() {
  cat <<'EOF'
Usage:
  bash run_o2_fit.sh --fitting_mode=invivo [options]
  bash run_o2_fit.sh --fitting_mode=invitro [options]
  bash run_o2_fit.sh --fitting_mode=joint --invivo_run_dir=DIR --invitro_run_dir=DIR [options]
  bash run_o2_fit.sh --fitting_mode=all [options]

Required modes:
  --fitting_mode=invivo|invitro|joint|all

Joint mode behavior:
  Build a pooled parameter-space t-SNE from the specified in-vivo and in-vitro
  source runs. Cluster the best points separately for each dataset, select the
  objective-minimum seed from every primary cluster on both sides, and pair all
  representatives by Cartesian product. No alternative joint mode, second-level
  cluster, or curve filter is available.

All mode behavior:
  Run in vivo first, then in vitro, then pass those generated result directories
  into the fixed joint primary-cluster workflow.

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
In-vitro and joint options:
  --parameter_table=/path/to/invitro_parameter_table.csv
  --fit_objects_dir=/path/to/fit_objects
  --flow_density_path=/path/to/g0g1_ploidy_density_grid.csv
  --itermax=1000 --itermax_max=1000

Joint options:
  --joint_run_prefix=name
  --joint_seeds_csv=1,2,3
  --joint_total_seeds=1
  --joint_soft_coupling_sigma_default=0.65
  --joint_soft_coupling_welsch_c=0.4
  --joint_warmup_sigmaN=0.0304
  --joint_soft_coupling_delta_params=default|all|none|param1,param2
  --joint_tsne_seed=123
  --joint_cluster_seed=123
  --joint_landscape_max_seeds=N
  --joint_landscape_n_threads=N

Local defaults run one seed per mode. Increase *_total_seeds or pass *_seeds_csv
explicitly when running a multi-seed local fit.
EOF
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

parse_args() {
  for arg in "$@"; do
    case "${arg}" in
      --help|-h)
        usage
        exit 0
        ;;
      --fitting_mode=*) FITTING_MODE="${arg#*=}" ;;
      --joint_fitting_mode=*)
        echo "--joint_fitting_mode has been removed; --fitting_mode=joint always uses the bilateral primary-cluster workflow." >&2
        exit 2
        ;;
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
      --passage_mode=*)
        echo "--passage_mode has been removed; in vitro always uses the fixed v2 passage implementation." >&2
        exit 2
        ;;
      --invivo_best_seed_dir=*|--joint_warmup_invivo_seed_dir=*|--joint_warmup_invivo_best_seed_dir=*|--invitro_best_seed_dir=*|--joint_warmup_invitro_seed_dir=*|--joint_warmup_invitro_best_seed_dir=*|--joint_warmup_vitro_seed_dir=*|--joint_warmup_enable=*|--joint_warmup_seed_label=*|--joint_seed_label=*|--seed_label=*|--joint_soft_coupling_parameters_table=*|--joint_soft_coupling_parameters_table_path=*|--select_required_files=*|--invivo_objective_columns=*|--invitro_objective_columns=*)
        echo "${arg%%=*} has been removed; joint anchors and coupling tables are derived from the specified source runs." >&2
        exit 2
        ;;
      --joint_warmup_sigmaN=*) JOINT_WARMUP_SIGMAN="${arg#*=}" ;;
      --joint_soft_coupling_sigma_default=*) JOINT_SOFT_COUPLING_SIGMA_DEFAULT="${arg#*=}" ;;
      --joint_soft_coupling_welsch_c=*) JOINT_SOFT_COUPLING_WELSCH_C="${arg#*=}" ;;
      --joint_soft_coupling_delta_params=*) JOINT_SOFT_COUPLING_DELTA_PARAMS="${arg#*=}" ;;
      --joint_tsne_seed=*|--multi_warmup_tsne_seed=*|--landscape_tsne_seed=*) JOINT_TSNE_SEED="${arg#*=}" ;;
      --joint_cluster_seed=*|--multi_warmup_cluster_seed=*|--landscape_cluster_seed=*) JOINT_CLUSTER_SEED="${arg#*=}" ;;
      --joint_landscape_max_seeds=*|--multi_warmup_landscape_max_seeds=*|--landscape_max_seeds=*) JOINT_LANDSCAPE_MAX_SEEDS="${arg#*=}" ;;
      --joint_landscape_n_threads=*|--landscape_n_threads=*) JOINT_LANDSCAPE_N_THREADS="${arg#*=}" ;;
      --multi_warmup_pair_method=*|--pair_method=*|--multi_warmup_reductions=*|--landscape_reductions=*|--multi_warmup_subcluster_seed=*|--landscape_subcluster_seed=*|--multi_warmup_pairing_policy=*|--pairing_policy=*|--multi_warmup_deduplicate_pairs=*|--deduplicate_pairs=*|--multi_warmup_reference_subcluster_dir=*|--reference_subcluster_dir=*)
        echo "${arg%%=*} has been removed; joint fitting now has one fixed bilateral primary-cluster workflow." >&2
        exit 2
        ;;
      --itermax=*) ITERMAX="${arg#*=}" ;;
      --itermax_max=*) ITERMAX_MAX="${arg#*=}" ;;
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
    "--itermax=${ITERMAX}"
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
      "--itermax_max=${ITERMAX_MAX}"
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

run_joint_primary_cluster_pipeline() {
  INVIVO_RUN_DIR="$(resolve_existing_dir "in vivo run directory" "${INVIVO_RUN_DIR}")"
  INVITRO_RUN_DIR="$(resolve_existing_dir "in vitro run directory" "${INVITRO_RUN_DIR}")"
  echo "Joint source in vivo run: ${INVIVO_RUN_DIR}"
  echo "Joint source in vitro run: ${INVITRO_RUN_DIR}"

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
    "--joint_cluster_seed=${JOINT_CLUSTER_SEED}"
    "--joint_tsne_seed=${JOINT_TSNE_SEED}"
    "--joint_landscape_n_threads=${JOINT_LANDSCAPE_N_THREADS}"
  )
  if [[ -n "${JOINT_LANDSCAPE_MAX_SEEDS}" ]]; then
    cmd+=("--joint_landscape_max_seeds=${JOINT_LANDSCAPE_MAX_SEEDS}")
  fi
  if [[ -n "${JOINT_WARMUP_SIGMAN}" ]]; then
    cmd+=("--joint_warmup_sigmaN=${JOINT_WARMUP_SIGMAN}")
  fi
  if [[ -n "${JOINT_SOFT_COUPLING_SIGMA_DEFAULT}" ]]; then
    cmd+=("--joint_soft_coupling_sigma_default=${JOINT_SOFT_COUPLING_SIGMA_DEFAULT}")
  fi
  if [[ -n "${JOINT_SOFT_COUPLING_WELSCH_C}" ]]; then
    cmd+=("--joint_soft_coupling_welsch_c=${JOINT_SOFT_COUPLING_WELSCH_C}")
  fi
  run_or_print "Run joint primary-cluster sweep" "${cmd[@]}"
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKFLOW_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
DEFAULT_PROJECT_ROOT="$(cd "${WORKFLOW_ROOT}/../../.." && pwd)"
DEFAULT_R_MODULE="R/4.4"
DEFAULT_INVIVO_RUN_PREFIX="fit_invivo_O2_buffering_local"
DEFAULT_INVITRO_RUN_PREFIX="fit_invitro_O2_buffering_local"
DEFAULT_JOINT_RUN_PREFIX="fit_joint_O2_buffering_local"
DEFAULT_TOTAL_SEEDS="1"
DEFAULT_N_CORES="1"
DEFAULT_AUTO_VIZ="TRUE"
DEFAULT_RUN_EXTRA_RESULTS="TRUE"
DEFAULT_FORCE_EXTRA_RESULTS="FALSE"
DEFAULT_DRY_RUN="FALSE"
DEFAULT_APPEND_RUN_PREFIX_TIMESTAMP="FALSE"
DEFAULT_ITERMAX="1000"
DEFAULT_ITERMAX_MAX="1000"
DEFAULT_DE_RELTOL="1e-4"
DEFAULT_DE_STEPTOL="25"
DEFAULT_NP="80"
DEFAULT_JOINT_WARMUP_SIGMAN=""
DEFAULT_JOINT_SOFT_COUPLING_SIGMA_DEFAULT=""
DEFAULT_JOINT_SOFT_COUPLING_WELSCH_C=""
DEFAULT_JOINT_SOFT_COUPLING_DELTA_PARAMS="default"
DEFAULT_JOINT_LANDSCAPE_MAX_SEEDS=""
DEFAULT_JOINT_CLUSTER_SEED="123"
DEFAULT_JOINT_TSNE_SEED="123"
DEFAULT_JOINT_LANDSCAPE_N_THREADS="8"

FITTING_MODE="${FITTING_MODE:-}"
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
ITERMAX_MAX="${ITERMAX_MAX:-}"
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
JOINT_WARMUP_SIGMAN="${JOINT_WARMUP_SIGMAN:-}"
JOINT_SOFT_COUPLING_SIGMA_DEFAULT="${JOINT_SOFT_COUPLING_SIGMA_DEFAULT:-}"
JOINT_SOFT_COUPLING_WELSCH_C="${JOINT_SOFT_COUPLING_WELSCH_C:-}"
JOINT_SOFT_COUPLING_DELTA_PARAMS="${JOINT_SOFT_COUPLING_DELTA_PARAMS:-}"
JOINT_LANDSCAPE_MAX_SEEDS="${JOINT_LANDSCAPE_MAX_SEEDS:-}"
JOINT_CLUSTER_SEED="${JOINT_CLUSTER_SEED:-}"
JOINT_TSNE_SEED="${JOINT_TSNE_SEED:-}"
JOINT_LANDSCAPE_N_THREADS="${JOINT_LANDSCAPE_N_THREADS:-}"

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
ITERMAX_MAX="${ITERMAX_MAX:-${DEFAULT_ITERMAX_MAX}}"
DE_RELTOL="${DE_RELTOL:-${DEFAULT_DE_RELTOL}}"
DE_STEPTOL="${DE_STEPTOL:-${DEFAULT_DE_STEPTOL}}"
NP="${NP:-${DEFAULT_NP}}"
AUTO_VIZ="${AUTO_VIZ:-${DEFAULT_AUTO_VIZ}}"
RUN_EXTRA_RESULTS="${RUN_EXTRA_RESULTS:-${DEFAULT_RUN_EXTRA_RESULTS}}"
FORCE_EXTRA_RESULTS="${FORCE_EXTRA_RESULTS:-${DEFAULT_FORCE_EXTRA_RESULTS}}"
DRY_RUN="${DRY_RUN:-${DEFAULT_DRY_RUN}}"
APPEND_RUN_PREFIX_TIMESTAMP="${APPEND_RUN_PREFIX_TIMESTAMP:-${DEFAULT_APPEND_RUN_PREFIX_TIMESTAMP}}"
JOINT_WARMUP_SIGMAN="${JOINT_WARMUP_SIGMAN:-${DEFAULT_JOINT_WARMUP_SIGMAN}}"
JOINT_SOFT_COUPLING_SIGMA_DEFAULT="${JOINT_SOFT_COUPLING_SIGMA_DEFAULT:-${DEFAULT_JOINT_SOFT_COUPLING_SIGMA_DEFAULT}}"
JOINT_SOFT_COUPLING_WELSCH_C="${JOINT_SOFT_COUPLING_WELSCH_C:-${DEFAULT_JOINT_SOFT_COUPLING_WELSCH_C}}"
JOINT_SOFT_COUPLING_DELTA_PARAMS="${JOINT_SOFT_COUPLING_DELTA_PARAMS:-${DEFAULT_JOINT_SOFT_COUPLING_DELTA_PARAMS}}"
JOINT_LANDSCAPE_MAX_SEEDS="${JOINT_LANDSCAPE_MAX_SEEDS:-${DEFAULT_JOINT_LANDSCAPE_MAX_SEEDS}}"
JOINT_CLUSTER_SEED="${JOINT_CLUSTER_SEED:-${DEFAULT_JOINT_CLUSTER_SEED}}"
JOINT_TSNE_SEED="${JOINT_TSNE_SEED:-${DEFAULT_JOINT_TSNE_SEED}}"
JOINT_LANDSCAPE_N_THREADS="${JOINT_LANDSCAPE_N_THREADS:-${DEFAULT_JOINT_LANDSCAPE_N_THREADS}}"

FITTING_MODE="$(normalize_fitting_mode "${FITTING_MODE}")"
if [[ -z "${FITTING_MODE}" ]]; then
  echo "--fitting_mode must be one of invivo, invitro, joint, or all." >&2
  usage >&2
  exit 2
fi

if [[ "${FITTING_MODE}" == "joint" ]]; then
  if is_null_value "${INVIVO_RUN_DIR}" || is_null_value "${INVITRO_RUN_DIR}"; then
    echo "--fitting_mode=joint requires both --invivo_run_dir and --invitro_run_dir." >&2
    exit 2
  fi
fi

require_positive_int INVIVO_TOTAL_SEEDS "${INVIVO_TOTAL_SEEDS}"
require_positive_int INVITRO_TOTAL_SEEDS "${INVITRO_TOTAL_SEEDS}"
require_positive_int JOINT_TOTAL_SEEDS "${JOINT_TOTAL_SEEDS}"
require_positive_int INVIVO_N_CORES "${INVIVO_N_CORES}"
require_positive_int INVITRO_N_CORES "${INVITRO_N_CORES}"
require_positive_int JOINT_N_CORES "${JOINT_N_CORES}"
require_positive_int ITERMAX "${ITERMAX}"
require_positive_int ITERMAX_MAX "${ITERMAX_MAX}"
require_positive_int DE_STEPTOL "${DE_STEPTOL}"
require_positive_int NP "${NP}"
require_positive_int JOINT_CLUSTER_SEED "${JOINT_CLUSTER_SEED}"
require_positive_int JOINT_TSNE_SEED "${JOINT_TSNE_SEED}"
require_positive_int JOINT_LANDSCAPE_N_THREADS "${JOINT_LANDSCAPE_N_THREADS}"
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
MULTI_WARMUP_RUNNER_SCRIPT="${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/runner/run_multi_warmup_joint.sh"
EXTRA_RESULTS_SCRIPT="${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/analysis/fit_results/extra_results.R"

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

for path in "${CONFIG_PATH}" "${FIT_RUNNER_SCRIPT}" "${MULTI_WARMUP_RUNNER_SCRIPT}" \
            "${EXTRA_RESULTS_SCRIPT}" "${PARAMETER_TABLE}"; do
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
    echo "Using the fixed bilateral primary-cluster Cartesian joint workflow."
    run_joint_primary_cluster_pipeline
    ;;
  all)
    echo "Running the complete chain: in vivo -> in vitro -> primary clusters -> joint fitting."
    run_invivo_fit
    run_extra_results "in vivo" "${INVIVO_RUN_DIR}"
    run_invitro_fit
    run_extra_results "in vitro" "${INVITRO_RUN_DIR}"
    run_joint_primary_cluster_pipeline
    ;;
esac
