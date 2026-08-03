#!/bin/bash
# Submit one joint array per selected multi-warm-up pair.

set -euo pipefail

O2SD_SHELL_UTILS="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../util" && pwd)/o2_supply_demand_map_shell_utils.sh"
# shellcheck source=../../util/o2_supply_demand_map_shell_utils.sh
source "${O2SD_SHELL_UTILS}"

usage() {
  cat <<'EOF'
Usage:
  bash submit_multi_warmup_joint.sh --invivo_run_dir=DIR --invitro_run_dir=DIR [options]

This launcher:
  1. Submits seed-plan / landscape clustering as a Slurm job by default.
  2. Submits a dependent controller that reads the generated manifest.
  3. Generates one joint soft-coupling parameter table per warm-up pair.
  4. Submits one joint seed array per warm-up pair, in manifest order.
  5. Submits per-pair extra_results jobs.
  6. Submits a final collector/report job that writes multi-warm-up_results.html.

Required:
  --invivo_run_dir=DIR
  --invitro_run_dir=DIR

Common options:
  --project_root=DIR
  --config_path=FILE
  --out_root=DIR
  --multi_warmup_prefix=NAME
  --joint_total_seeds=500
  --joint_array_tasks=500
  --joint_seeds_per_task=1
  --joint_qos=xxlarge
  --joint_time_limit=12:00:00
  --joint_soft_coupling_welsch_c=0.4
  --multi_warmup_pair_method=legacy|landscape_subcluster
  --multi_warmup_invivo_top_n=10  (0 disables in vivo source clustering; not both sides)
  --multi_warmup_invitro_top_n=10 (0 disables in vitro source clustering; not both sides)
  --multi_warmup_umap_seed=1
  --multi_warmup_reductions=tsne,umap
  --multi_warmup_landscape_umap_seed=123
  --multi_warmup_landscape_max_seeds=N
  --multi_warmup_pairing_policy=cartesian_by_method|invitro_best_to_invivo_subclusters
  --multi_warmup_deduplicate_pairs=FALSE
  --multi_warmup_invivo_curve_filter=TRUE
  --multi_warmup_invivo_curve_class=monotone_increasing
  --multi_warmup_parallel_seed_space=TRUE
  --multi_warmup_seed_space_qos=small
  --multi_warmup_seed_space_time_limit=2:00:00
  --multi_warmup_seed_space_mem=4G
  --multi_warmup_seed_space_cpus=1
  --multi_warmup_seed_plan_as_job=TRUE
  --multi_warmup_seed_plan_qos=small
  --multi_warmup_seed_plan_time_limit=12:00:00
  --multi_warmup_seed_plan_mem=32G
  --multi_warmup_seed_plan_cpus=1
  --multi_warmup_monotonicity_qos=small
  --multi_warmup_monotonicity_time_limit=4:00:00
  --multi_warmup_monotonicity_mem=8G
  --multi_warmup_monotonicity_tasks_per_array_task=25
  --multi_warmup_submit_qos=small
  --multi_warmup_submit_time_limit=4:00:00
  --multi_warmup_submit_mem=8G
  --dry_run=TRUE|FALSE
EOF
}

parse_args() {
  for arg in "$@"; do
    case "${arg}" in
      --help|-h) usage; exit 0 ;;
      --project_root=*) PROJECT_ROOT="${arg#*=}" ;;
      --config_path=*|--config=*) CONFIG_PATH="${arg#*=}" ;;
      --out_root=*) OUT_ROOT="${arg#*=}" ;;
      --multi_warmup_prefix=*|--joint_run_prefix=*) MULTI_WARMUP_PREFIX="${arg#*=}" ;;
      --invivo_run_dir=*) INVIVO_RUN_DIR="${arg#*=}" ;;
      --invitro_run_dir=*) INVITRO_RUN_DIR="${arg#*=}" ;;
      --joint_total_seeds=*) JOINT_TOTAL_SEEDS="${arg#*=}" ;;
      --joint_array_tasks=*) JOINT_ARRAY_TASKS="${arg#*=}" ;;
      --joint_seeds_per_task=*) JOINT_SEEDS_PER_TASK="${arg#*=}" ;;
      --joint_n_cores=*|--n_cores=*) JOINT_N_CORES="${arg#*=}" ;;
      --joint_mem=*|--mem=*) JOINT_MEM="${arg#*=}" ;;
      --joint_qos=*|--qos=*) JOINT_QOS="${arg#*=}" ;;
      --joint_time_limit=*|--time=*|--time_limit=*) JOINT_TIME_LIMIT="${arg#*=}" ;;
      --postprocess_qos=*) POSTPROCESS_QOS="${arg#*=}" ;;
      --postprocess_time_limit=*) POSTPROCESS_TIME_LIMIT="${arg#*=}" ;;
      --postprocess_mem=*) POSTPROCESS_MEM="${arg#*=}" ;;
      --report_qos=*) REPORT_QOS="${arg#*=}" ;;
      --report_time_limit=*) REPORT_TIME_LIMIT="${arg#*=}" ;;
      --report_mem=*) REPORT_MEM="${arg#*=}" ;;
      --parameter_table=*|--invitro_parameter_table=*) PARAMETER_TABLE="${arg#*=}" ;;
      --fit_objects_dir=*) FIT_OBJECTS_DIR="${arg#*=}" ;;
      --flow_density_path=*) FLOW_DENSITY_PATH="${arg#*=}" ;;
      --death_data_path=*|--invitro_death_data_path=*) DEATH_DATA_PATH="${arg#*=}" ;;
      --death_weight=*|--joint_invitro_death_weight=*) DEATH_WEIGHT="${arg#*=}" ;;
      --sigma_death_logit=*|--invitro_sigma_death_logit=*) SIGMA_DEATH_LOGIT="${arg#*=}" ;;
      --death_fraction_eps=*|--invitro_death_fraction_eps=*) DEATH_FRACTION_EPS="${arg#*=}" ;;
      --itermax=*) ITERMAX="${arg#*=}" ;;
      --de_reltol=*) DE_RELTOL="${arg#*=}" ;;
      --de_steptol=*) DE_STEPTOL="${arg#*=}" ;;
      --np=*|--NP=*) NP="${arg#*=}" ;;
      --auto_viz=*) AUTO_VIZ="${arg#*=}" ;;
      --r_module=*) R_MODULE="${arg#*=}" ;;
      --force_extra_results=*) FORCE_EXTRA_RESULTS="${arg#*=}" ;;
      --dry_run=*) DRY_RUN="${arg#*=}" ;;
      --joint_warmup_sigmaN=*) JOINT_WARMUP_SIGMAN="${arg#*=}" ;;
      --joint_soft_coupling_sigma_default=*) JOINT_SOFT_COUPLING_SIGMA_DEFAULT="${arg#*=}" ;;
      --joint_soft_coupling_welsch_c=*) JOINT_SOFT_COUPLING_WELSCH_C="${arg#*=}" ;;
      --joint_soft_coupling_delta_params=*) JOINT_SOFT_COUPLING_DELTA_PARAMS="${arg#*=}" ;;
      --multi_warmup_pair_method=*|--pair_method=*) MULTI_WARMUP_PAIR_METHOD="${arg#*=}" ;;
      --multi_warmup_top_n=*) MULTI_WARMUP_TOP_N="${arg#*=}" ;;
      --multi_warmup_invivo_top_n=*|--invivo_top_n=*) MULTI_WARMUP_INVIVO_TOP_N="${arg#*=}" ;;
      --multi_warmup_invitro_top_n=*|--invitro_top_n=*) MULTI_WARMUP_INVITRO_TOP_N="${arg#*=}" ;;
      --multi_warmup_umap_seed=*|--umap_seed=*) MULTI_WARMUP_UMAP_SEED="${arg#*=}" ;;
      --multi_warmup_invivo_k=*) MULTI_WARMUP_INVIVO_K="${arg#*=}" ;;
      --multi_warmup_invitro_k=*) MULTI_WARMUP_INVITRO_K="${arg#*=}" ;;
      --multi_warmup_invitro_anchor_ranks=*) MULTI_WARMUP_INVITRO_ANCHOR_RANKS="${arg#*=}" ;;
      --multi_warmup_include_phase2=*) MULTI_WARMUP_INCLUDE_PHASE2="${arg#*=}" ;;
      --multi_warmup_phase2_invitro_anchor_ranks=*) MULTI_WARMUP_PHASE2_INVITRO_ANCHOR_RANKS="${arg#*=}" ;;
      --multi_warmup_reductions=*|--landscape_reductions=*) MULTI_WARMUP_REDUCTIONS="${arg#*=}" ;;
      --multi_warmup_landscape_umap_seed=*|--landscape_umap_seed=*) MULTI_WARMUP_LANDSCAPE_UMAP_SEED="${arg#*=}" ;;
      --multi_warmup_landscape_max_seeds=*|--landscape_max_seeds=*) MULTI_WARMUP_LANDSCAPE_MAX_SEEDS="${arg#*=}" ;;
      --multi_warmup_cluster_seed=*|--landscape_cluster_seed=*) MULTI_WARMUP_CLUSTER_SEED="${arg#*=}" ;;
      --multi_warmup_subcluster_seed=*|--landscape_subcluster_seed=*) MULTI_WARMUP_SUBCLUSTER_SEED="${arg#*=}" ;;
      --multi_warmup_tsne_seed=*|--landscape_tsne_seed=*) MULTI_WARMUP_TSNE_SEED="${arg#*=}" ;;
      --multi_warmup_pairing_policy=*|--pairing_policy=*) MULTI_WARMUP_PAIRING_POLICY="${arg#*=}" ;;
      --multi_warmup_deduplicate_pairs=*|--deduplicate_pairs=*) MULTI_WARMUP_DEDUPLICATE_PAIRS="${arg#*=}" ;;
      --multi_warmup_invivo_curve_filter=*|--invivo_curve_filter=*) MULTI_WARMUP_INVIVO_CURVE_FILTER="${arg#*=}" ;;
      --multi_warmup_invivo_curve_class=*|--invivo_curve_class=*) MULTI_WARMUP_INVIVO_CURVE_CLASS="${arg#*=}" ;;
      --multi_warmup_parallel_seed_space=*|--parallel_seed_space=*) MULTI_WARMUP_PARALLEL_SEED_SPACE="${arg#*=}" ;;
      --multi_warmup_seed_space_qos=*|--seed_space_qos=*) MULTI_WARMUP_SEED_SPACE_QOS="${arg#*=}" ;;
      --multi_warmup_seed_space_time_limit=*|--seed_space_time_limit=*) MULTI_WARMUP_SEED_SPACE_TIME_LIMIT="${arg#*=}" ;;
      --multi_warmup_seed_space_mem=*|--seed_space_mem=*) MULTI_WARMUP_SEED_SPACE_MEM="${arg#*=}" ;;
      --multi_warmup_seed_space_cpus=*|--seed_space_cpus=*) MULTI_WARMUP_SEED_SPACE_CPUS="${arg#*=}" ;;
      --multi_warmup_invivo_best_csv=*|--invivo_best_csv=*) MULTI_WARMUP_INVIVO_BEST_CSV="${arg#*=}" ;;
      --multi_warmup_invivo_initial_csv=*|--invivo_initial_csv=*) MULTI_WARMUP_INVIVO_INITIAL_CSV="${arg#*=}" ;;
      --multi_warmup_invitro_best_csv=*|--invitro_best_csv=*) MULTI_WARMUP_INVITRO_BEST_CSV="${arg#*=}" ;;
      --multi_warmup_invitro_initial_csv=*|--invitro_initial_csv=*) MULTI_WARMUP_INVITRO_INITIAL_CSV="${arg#*=}" ;;
      --multi_warmup_monotonicity_qos=*|--monotonicity_qos=*) MULTI_WARMUP_MONOTONICITY_QOS="${arg#*=}" ;;
      --multi_warmup_monotonicity_time_limit=*|--monotonicity_time_limit=*) MULTI_WARMUP_MONOTONICITY_TIME_LIMIT="${arg#*=}" ;;
      --multi_warmup_monotonicity_mem=*|--monotonicity_mem=*) MULTI_WARMUP_MONOTONICITY_MEM="${arg#*=}" ;;
      --multi_warmup_monotonicity_cpus=*|--monotonicity_cpus=*) MULTI_WARMUP_MONOTONICITY_CPUS="${arg#*=}" ;;
      --multi_warmup_monotonicity_tasks_per_array_task=*|--monotonicity_tasks_per_array_task=*) MULTI_WARMUP_MONOTONICITY_TASKS_PER_ARRAY_TASK="${arg#*=}" ;;
      --multi_warmup_validation_qos=*|--validation_qos=*) MULTI_WARMUP_VALIDATION_QOS="${arg#*=}" ;;
      --multi_warmup_validation_time_limit=*|--validation_time_limit=*) MULTI_WARMUP_VALIDATION_TIME_LIMIT="${arg#*=}" ;;
      --multi_warmup_validation_mem=*|--validation_mem=*) MULTI_WARMUP_VALIDATION_MEM="${arg#*=}" ;;
      --multi_warmup_reference_subcluster_dir=*|--reference_subcluster_dir=*) MULTI_WARMUP_REFERENCE_SUBCLUSTER_DIR="${arg#*=}" ;;
      --multi_warmup_seed_plan_as_job=*|--seed_plan_as_job=*) MULTI_WARMUP_SEED_PLAN_AS_JOB="${arg#*=}" ;;
      --multi_warmup_seed_plan_qos=*|--seed_plan_qos=*) MULTI_WARMUP_SEED_PLAN_QOS="${arg#*=}" ;;
      --multi_warmup_seed_plan_time_limit=*|--seed_plan_time_limit=*) MULTI_WARMUP_SEED_PLAN_TIME_LIMIT="${arg#*=}" ;;
      --multi_warmup_seed_plan_mem=*|--seed_plan_mem=*) MULTI_WARMUP_SEED_PLAN_MEM="${arg#*=}" ;;
      --multi_warmup_seed_plan_cpus=*|--seed_plan_cpus=*) MULTI_WARMUP_SEED_PLAN_CPUS="${arg#*=}" ;;
      --multi_warmup_submit_qos=*|--submit_qos=*) MULTI_WARMUP_SUBMIT_QOS="${arg#*=}" ;;
      --multi_warmup_submit_time_limit=*|--submit_time_limit=*) MULTI_WARMUP_SUBMIT_TIME_LIMIT="${arg#*=}" ;;
      --multi_warmup_submit_mem=*|--submit_mem=*) MULTI_WARMUP_SUBMIT_MEM="${arg#*=}" ;;
      --internal_stage=*) INTERNAL_STAGE="${arg#*=}" ;;
      --log_root=*|--log_dir=*) LOG_ROOT="${arg#*=}" ;;
      *) echo "Unknown argument: ${arg}" >&2; usage >&2; exit 2 ;;
    esac
  done
}

submit_or_print() {
  local label="$1"
  shift
  {
    printf "%s:" "${label}"
    printf " %q" "$@"
    printf "\n"
  } | tee -a "${PROGRESS_LOG}" >&2
  if truthy "${DRY_RUN}"; then
    echo "DRY_RUN_JOB_ID"
  else
    "$@" | awk '{print $NF}'
  fi
}

ORIGINAL_ARGS=("$@")
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SELF_SCRIPT="${SCRIPT_DIR}/$(basename "${BASH_SOURCE[0]}")"
WORKFLOW_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
DEFAULT_PROJECT_ROOT="$(cd "${WORKFLOW_ROOT}/../../.." && pwd)"

PROJECT_ROOT="${PROJECT_ROOT:-${DEFAULT_PROJECT_ROOT}}"
CONFIG_PATH="${CONFIG_PATH:-}"
OUT_ROOT="${OUT_ROOT:-}"
MULTI_WARMUP_PREFIX="${MULTI_WARMUP_PREFIX:-fit_joint_multi_warmup}"
INVIVO_RUN_DIR="${INVIVO_RUN_DIR:-}"
INVITRO_RUN_DIR="${INVITRO_RUN_DIR:-}"
JOINT_TOTAL_SEEDS="${JOINT_TOTAL_SEEDS:-500}"
JOINT_ARRAY_TASKS="${JOINT_ARRAY_TASKS:-500}"
JOINT_SEEDS_PER_TASK="${JOINT_SEEDS_PER_TASK:-1}"
JOINT_N_CORES="${JOINT_N_CORES:-22}"
JOINT_MEM="${JOINT_MEM:-32G}"
JOINT_QOS="${JOINT_QOS:-xxlarge}"
JOINT_TIME_LIMIT="${JOINT_TIME_LIMIT:-12:00:00}"
POSTPROCESS_QOS="${POSTPROCESS_QOS:-small}"
POSTPROCESS_TIME_LIMIT="${POSTPROCESS_TIME_LIMIT:-4:00:00}"
POSTPROCESS_MEM="${POSTPROCESS_MEM:-8G}"
REPORT_QOS="${REPORT_QOS:-small}"
REPORT_TIME_LIMIT="${REPORT_TIME_LIMIT:-4:00:00}"
REPORT_MEM="${REPORT_MEM:-8G}"
PARAMETER_TABLE="${PARAMETER_TABLE:-}"
FIT_OBJECTS_DIR="${FIT_OBJECTS_DIR:-}"
FLOW_DENSITY_PATH="${FLOW_DENSITY_PATH:-}"
DEATH_DATA_PATH="${DEATH_DATA_PATH:-}"
DEATH_WEIGHT="${DEATH_WEIGHT:-1}"
SIGMA_DEATH_LOGIT="${SIGMA_DEATH_LOGIT:-0.75}"
DEATH_FRACTION_EPS="${DEATH_FRACTION_EPS:-1e-4}"
ITERMAX="${ITERMAX:-500}"
DE_RELTOL="${DE_RELTOL:-1e-4}"
DE_STEPTOL="${DE_STEPTOL:-25}"
NP="${NP:-80}"
AUTO_VIZ="${AUTO_VIZ:-TRUE}"
R_MODULE="${R_MODULE:-R/4.4}"
FORCE_EXTRA_RESULTS="${FORCE_EXTRA_RESULTS:-FALSE}"
DRY_RUN="${DRY_RUN:-FALSE}"
JOINT_WARMUP_SIGMAN="${JOINT_WARMUP_SIGMAN:-}"
JOINT_SOFT_COUPLING_SIGMA_DEFAULT="${JOINT_SOFT_COUPLING_SIGMA_DEFAULT:-}"
JOINT_SOFT_COUPLING_WELSCH_C="${JOINT_SOFT_COUPLING_WELSCH_C:-}"
JOINT_SOFT_COUPLING_DELTA_PARAMS="${JOINT_SOFT_COUPLING_DELTA_PARAMS:-default}"
MULTI_WARMUP_PAIR_METHOD="${MULTI_WARMUP_PAIR_METHOD:-legacy}"
MULTI_WARMUP_TOP_N="${MULTI_WARMUP_TOP_N:-10}"
MULTI_WARMUP_INVIVO_TOP_N="${MULTI_WARMUP_INVIVO_TOP_N:-${MULTI_WARMUP_TOP_N}}"
MULTI_WARMUP_INVITRO_TOP_N="${MULTI_WARMUP_INVITRO_TOP_N:-${MULTI_WARMUP_TOP_N}}"
MULTI_WARMUP_UMAP_SEED="${MULTI_WARMUP_UMAP_SEED:-1}"
MULTI_WARMUP_INVIVO_K="${MULTI_WARMUP_INVIVO_K:-auto}"
MULTI_WARMUP_INVITRO_K="${MULTI_WARMUP_INVITRO_K:-auto}"
MULTI_WARMUP_INVITRO_ANCHOR_RANKS="${MULTI_WARMUP_INVITRO_ANCHOR_RANKS:-1}"
MULTI_WARMUP_INCLUDE_PHASE2="${MULTI_WARMUP_INCLUDE_PHASE2:-FALSE}"
MULTI_WARMUP_PHASE2_INVITRO_ANCHOR_RANKS="${MULTI_WARMUP_PHASE2_INVITRO_ANCHOR_RANKS:-auto}"
MULTI_WARMUP_REDUCTIONS="${MULTI_WARMUP_REDUCTIONS:-tsne,umap}"
MULTI_WARMUP_LANDSCAPE_UMAP_SEED="${MULTI_WARMUP_LANDSCAPE_UMAP_SEED:-123}"
MULTI_WARMUP_LANDSCAPE_MAX_SEEDS="${MULTI_WARMUP_LANDSCAPE_MAX_SEEDS:-}"
MULTI_WARMUP_CLUSTER_SEED="${MULTI_WARMUP_CLUSTER_SEED:-123}"
MULTI_WARMUP_SUBCLUSTER_SEED="${MULTI_WARMUP_SUBCLUSTER_SEED:-1123}"
MULTI_WARMUP_TSNE_SEED="${MULTI_WARMUP_TSNE_SEED:-123}"
MULTI_WARMUP_PAIRING_POLICY="${MULTI_WARMUP_PAIRING_POLICY:-cartesian_by_method}"
MULTI_WARMUP_DEDUPLICATE_PAIRS="${MULTI_WARMUP_DEDUPLICATE_PAIRS:-FALSE}"
MULTI_WARMUP_INVIVO_CURVE_FILTER="${MULTI_WARMUP_INVIVO_CURVE_FILTER:-TRUE}"
MULTI_WARMUP_INVIVO_CURVE_CLASS="${MULTI_WARMUP_INVIVO_CURVE_CLASS:-monotone_increasing}"
MULTI_WARMUP_PARALLEL_SEED_SPACE="${MULTI_WARMUP_PARALLEL_SEED_SPACE:-TRUE}"
MULTI_WARMUP_SEED_SPACE_QOS="${MULTI_WARMUP_SEED_SPACE_QOS:-small}"
MULTI_WARMUP_SEED_SPACE_TIME_LIMIT="${MULTI_WARMUP_SEED_SPACE_TIME_LIMIT:-2:00:00}"
MULTI_WARMUP_SEED_SPACE_MEM="${MULTI_WARMUP_SEED_SPACE_MEM:-4G}"
MULTI_WARMUP_SEED_SPACE_CPUS="${MULTI_WARMUP_SEED_SPACE_CPUS:-1}"
MULTI_WARMUP_INVIVO_BEST_CSV="${MULTI_WARMUP_INVIVO_BEST_CSV:-}"
MULTI_WARMUP_INVIVO_INITIAL_CSV="${MULTI_WARMUP_INVIVO_INITIAL_CSV:-}"
MULTI_WARMUP_INVITRO_BEST_CSV="${MULTI_WARMUP_INVITRO_BEST_CSV:-}"
MULTI_WARMUP_INVITRO_INITIAL_CSV="${MULTI_WARMUP_INVITRO_INITIAL_CSV:-}"
MULTI_WARMUP_MONOTONICITY_QOS="${MULTI_WARMUP_MONOTONICITY_QOS:-small}"
MULTI_WARMUP_MONOTONICITY_TIME_LIMIT="${MULTI_WARMUP_MONOTONICITY_TIME_LIMIT:-4:00:00}"
MULTI_WARMUP_MONOTONICITY_MEM="${MULTI_WARMUP_MONOTONICITY_MEM:-8G}"
MULTI_WARMUP_MONOTONICITY_CPUS="${MULTI_WARMUP_MONOTONICITY_CPUS:-1}"
MULTI_WARMUP_MONOTONICITY_TASKS_PER_ARRAY_TASK="${MULTI_WARMUP_MONOTONICITY_TASKS_PER_ARRAY_TASK:-25}"
MULTI_WARMUP_VALIDATION_QOS="${MULTI_WARMUP_VALIDATION_QOS:-small}"
MULTI_WARMUP_VALIDATION_TIME_LIMIT="${MULTI_WARMUP_VALIDATION_TIME_LIMIT:-4:00:00}"
MULTI_WARMUP_VALIDATION_MEM="${MULTI_WARMUP_VALIDATION_MEM:-16G}"
MULTI_WARMUP_REFERENCE_SUBCLUSTER_DIR="${MULTI_WARMUP_REFERENCE_SUBCLUSTER_DIR:-}"
MULTI_WARMUP_SEED_PLAN_AS_JOB="${MULTI_WARMUP_SEED_PLAN_AS_JOB:-TRUE}"
MULTI_WARMUP_SEED_PLAN_QOS="${MULTI_WARMUP_SEED_PLAN_QOS:-small}"
MULTI_WARMUP_SEED_PLAN_TIME_LIMIT="${MULTI_WARMUP_SEED_PLAN_TIME_LIMIT:-12:00:00}"
MULTI_WARMUP_SEED_PLAN_MEM="${MULTI_WARMUP_SEED_PLAN_MEM:-32G}"
MULTI_WARMUP_SEED_PLAN_CPUS="${MULTI_WARMUP_SEED_PLAN_CPUS:-1}"
MULTI_WARMUP_SUBMIT_QOS="${MULTI_WARMUP_SUBMIT_QOS:-small}"
MULTI_WARMUP_SUBMIT_TIME_LIMIT="${MULTI_WARMUP_SUBMIT_TIME_LIMIT:-4:00:00}"
MULTI_WARMUP_SUBMIT_MEM="${MULTI_WARMUP_SUBMIT_MEM:-8G}"
INTERNAL_STAGE="${INTERNAL_STAGE:-}"
LOG_ROOT="${LOG_ROOT:-}"

parse_args "$@"

PASSTHROUGH_ARGS=()
for arg in "${ORIGINAL_ARGS[@]}"; do
  case "${arg}" in
    --internal_stage=*|--dry_run=*|--multi_warmup_seed_plan_as_job=*|--seed_plan_as_job=*) ;;
    *) PASSTHROUGH_ARGS+=("${arg}") ;;
  esac
done

case "${INTERNAL_STAGE}" in
  ""|seed_plan|curve_class_workflow|submit_pairs) ;;
  *) echo "--internal_stage must be seed_plan, curve_class_workflow, or submit_pairs, got: ${INTERNAL_STAGE}" >&2; exit 2 ;;
esac

MULTI_WARMUP_PAIR_METHOD="$(echo "${MULTI_WARMUP_PAIR_METHOD}" | tr '[:upper:]' '[:lower:]' | tr '-' '_')"
case "${MULTI_WARMUP_PAIR_METHOD}" in
  legacy|landscape_subcluster) ;;
  *) echo "--multi_warmup_pair_method must be legacy or landscape_subcluster, got: ${MULTI_WARMUP_PAIR_METHOD}" >&2; exit 2 ;;
esac
require_nonnegative_int MULTI_WARMUP_INVIVO_TOP_N "${MULTI_WARMUP_INVIVO_TOP_N}"
require_nonnegative_int MULTI_WARMUP_INVITRO_TOP_N "${MULTI_WARMUP_INVITRO_TOP_N}"
require_nonnegative_int MULTI_WARMUP_MONOTONICITY_TASKS_PER_ARRAY_TASK "${MULTI_WARMUP_MONOTONICITY_TASKS_PER_ARRAY_TASK}"
if (( MULTI_WARMUP_MONOTONICITY_TASKS_PER_ARRAY_TASK < 1 )); then
  echo "--multi_warmup_monotonicity_tasks_per_array_task must be >= 1." >&2
  exit 2
fi
if (( MULTI_WARMUP_INVIVO_TOP_N == 0 && MULTI_WARMUP_INVITRO_TOP_N == 0 )); then
  echo "At least one of --invivo_top_n or --invitro_top_n must be greater than 0." >&2
  exit 2
fi
if [[ "${MULTI_WARMUP_PAIR_METHOD}" == "landscape_subcluster" ]] && (( MULTI_WARMUP_INVIVO_TOP_N == 0 || MULTI_WARMUP_INVITRO_TOP_N == 0 )); then
  echo "landscape_subcluster pair method requires both in vivo and in vitro source runs." >&2
  exit 2
fi
if (( MULTI_WARMUP_INVIVO_TOP_N > 0 )) && is_null_value "${INVIVO_RUN_DIR}"; then
  echo "HPC multi-warm-up submission requires existing --invivo_run_dir when --invivo_top_n > 0." >&2
  echo "Submit the source in vivo fit first, then rerun this launcher." >&2
  exit 2
fi
if (( MULTI_WARMUP_INVITRO_TOP_N > 0 )) && is_null_value "${INVITRO_RUN_DIR}"; then
  echo "HPC multi-warm-up submission requires existing --invitro_run_dir when --invitro_top_n > 0." >&2
  echo "Submit the source in vitro fit first, then rerun this launcher." >&2
  exit 2
fi
if (( JOINT_ARRAY_TASKS * JOINT_SEEDS_PER_TASK != JOINT_TOTAL_SEEDS )); then
  echo "joint_array_tasks * joint_seeds_per_task must equal joint_total_seeds." >&2
  exit 2
fi
if ! truthy "${DRY_RUN}" && [[ "${INTERNAL_STAGE}" != "seed_plan" ]] && ! command -v sbatch >/dev/null 2>&1; then
  echo "sbatch not found; run on an HPC login node or use --dry_run=TRUE." >&2
  exit 1
fi

PROJECT_ROOT="$(cd "${PROJECT_ROOT}" && pwd)"
if [[ -z "${CONFIG_PATH}" ]]; then CONFIG_PATH="${PROJECT_ROOT}/oxygen/config/O2_supply_demand.yaml"; fi
if [[ -z "${OUT_ROOT}" ]]; then OUT_ROOT="${PROJECT_ROOT}/oxygen/results"; fi
if [[ -z "${PARAMETER_TABLE}" ]]; then PARAMETER_TABLE="${PROJECT_ROOT}/oxygen/data/O2_supply_demand/parameter_table_invitro_buffering.csv"; fi
if [[ -z "${FIT_OBJECTS_DIR}" ]]; then FIT_OBJECTS_DIR="${PROJECT_ROOT}/oxygen/ploidyOxygen/data/fit_objects"; fi
if [[ -z "${FLOW_DENSITY_PATH}" ]]; then FLOW_DENSITY_PATH="${PROJECT_ROOT}/oxygen/data/g0g1_ploidy_density_grid.csv"; fi
if [[ -z "${DEATH_DATA_PATH}" ]]; then DEATH_DATA_PATH="${PROJECT_ROOT}/data/InVitroData/sum159_dead_cell_endpoint_likelihood_ready_20260731.tsv"; fi
CONFIG_PATH="$(cd "$(dirname "${CONFIG_PATH}")" && pwd)/$(basename "${CONFIG_PATH}")"
OUT_ROOT="$(mkdir -p "${OUT_ROOT}" && cd "${OUT_ROOT}" && pwd)"
if (( MULTI_WARMUP_INVIVO_TOP_N > 0 )); then INVIVO_RUN_DIR="$(cd "${INVIVO_RUN_DIR}" && pwd)"; else INVIVO_RUN_DIR=""; fi
if (( MULTI_WARMUP_INVITRO_TOP_N > 0 )); then INVITRO_RUN_DIR="$(cd "${INVITRO_RUN_DIR}" && pwd)"; else INVITRO_RUN_DIR=""; fi
PARAMETER_TABLE="$(cd "$(dirname "${PARAMETER_TABLE}")" && pwd)/$(basename "${PARAMETER_TABLE}")"
FIT_OBJECTS_DIR="$(cd "${FIT_OBJECTS_DIR}" && pwd)"
FLOW_DENSITY_PATH="$(cd "$(dirname "${FLOW_DENSITY_PATH}")" && pwd)/$(basename "${FLOW_DENSITY_PATH}")"
DEATH_DATA_PATH="$(cd "$(dirname "${DEATH_DATA_PATH}")" && pwd)/$(basename "${DEATH_DATA_PATH}")"
if [[ -z "${LOG_ROOT}" ]]; then LOG_ROOT="${OUT_ROOT}/log"; fi
LOG_ROOT="$(mkdir -p "${LOG_ROOT}" && cd "${LOG_ROOT}" && pwd)"

MULTI_WARMUP_ROOT="${OUT_ROOT}/${MULTI_WARMUP_PREFIX}"
mkdir -p "${MULTI_WARMUP_ROOT}"
PROGRESS_LOG="${MULTI_WARMUP_ROOT}/multi_warmup_progress.log"
JOBS_TSV="${MULTI_WARMUP_ROOT}/multi_warmup_jobs.tsv"
CONTROLLER_JOBS_TSV="${MULTI_WARMUP_ROOT}/multi_warmup_controller_jobs.tsv"
if [[ -z "${INTERNAL_STAGE}" ]]; then
  : > "${PROGRESS_LOG}"
else
  touch "${PROGRESS_LOG}"
fi

LEGACY_SEED_PLAN_SCRIPT="${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/analysis/multi_warmup/build_multi_warmup_seed_plan.R"
LANDSCAPE_SEED_PLAN_SCRIPT="${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/analysis/multi_warmup/build_multi_warmup_pairs_from_landscape_subclusters.R"
if [[ "${MULTI_WARMUP_PAIR_METHOD}" == "landscape_subcluster" ]]; then
  SEED_PLAN_SCRIPT="${LANDSCAPE_SEED_PLAN_SCRIPT}"
else
  SEED_PLAN_SCRIPT="${LEGACY_SEED_PLAN_SCRIPT}"
fi
MAKE_TABLE_SCRIPT="${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/runner/warm_start/make_joint_soft_coupling_parameters_table.R"
COLLECT_SCRIPT="${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/analysis/multi_warmup/collect_multi_warmup_results.R"
REPORT_SCRIPT="${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/analysis/multi_warmup/multi_warmup_results_report.R"
JOINT_ARRAY_SCRIPT="${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/hpc/array_workers/submit_fit_seed_array_joint_buffering.sub"
JOINT_RUNNER_SCRIPT="${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/runner/run_fit_joint_model_O2_supply_demand_MAP.sh"
POSTPROCESS_SCRIPT="${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/hpc/postprocess/postprocess_extra_results.sh"
DENSE_GRID_SCRIPT="${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/runner/dense_grid_monotonicity/run_dense_grid_monotonicity.R"

required_paths=("${SEED_PLAN_SCRIPT}" "${MAKE_TABLE_SCRIPT}" "${COLLECT_SCRIPT}" "${REPORT_SCRIPT}" "${JOINT_ARRAY_SCRIPT}" "${JOINT_RUNNER_SCRIPT}" "${POSTPROCESS_SCRIPT}" "${CONFIG_PATH}" "${PARAMETER_TABLE}")
if [[ "${MULTI_WARMUP_PAIR_METHOD}" == "landscape_subcluster" ]]; then
  required_paths+=("${DENSE_GRID_SCRIPT}")
fi
for path in "${required_paths[@]}"; do
  if [[ ! -f "${path}" ]]; then
    echo "Missing required file: ${path}" >&2
    exit 1
  fi
done

wrap_with_r_body() {
  local body="$1"
  local root_q module_q
  printf -v root_q '%q' "${PROJECT_ROOT}"
  printf -v module_q '%q' "${R_MODULE}"
  shell_join bash -lc "if command -v ml >/dev/null 2>&1; then ml ${module_q}; elif command -v module >/dev/null 2>&1; then module load ${module_q}; fi; cd ${root_q} && ${body}"
}

wrap_with_r_command() {
  local cmd
  cmd="$(shell_join "$@")"
  wrap_with_r_body "${cmd}"
}

count_valid_seed_dirs() {
  local input_dir="$1"
  local max_seeds="${2:-}"
  local count=0
  local path base
  while IFS= read -r -d '' path; do
    base="$(basename "${path}")"
    if [[ "${base}" =~ ^seed[0-9]+$ ]] && [[ -f "${path}/fit_summary.tsv" ]] && [[ -f "${path}/best_params.tsv" ]]; then
      count=$((count + 1))
    fi
  done < <(find "${input_dir}" -maxdepth 1 -mindepth 1 -type d -name 'seed*' -print0 | sort -z)
  if [[ "${max_seeds}" =~ ^[0-9]+$ ]] && (( max_seeds > 0 && count > max_seeds )); then
    count="${max_seeds}"
  fi
  echo "${count}"
}

submit_parallel_seed_space_landscape_workflow() {
  local seed_space_root="${MULTI_WARMUP_ROOT}/seed_space_shards"
  local seed_table_dir="${MULTI_WARMUP_ROOT}/landscape_subcluster/SeedParameterTables"
  local invivo_tasks="${seed_space_root}/invivo_seed_space_tasks.tsv"
  local invitro_tasks="${seed_space_root}/invitro_seed_space_tasks.tsv"
  local invivo_best_csv="${seed_table_dir}/invivo_best_params_by_seed.csv"
  local invivo_initial_csv="${seed_table_dir}/invivo_deoptim_initial_population.csv"
  local invitro_best_csv="${seed_table_dir}/invitro_best_params_by_seed.csv"
  local invitro_initial_csv="${seed_table_dir}/invitro_deoptim_initial_population.csv"
  local invivo_count invitro_count
  invivo_count="$(count_valid_seed_dirs "${INVIVO_RUN_DIR}" "${MULTI_WARMUP_LANDSCAPE_MAX_SEEDS}")"
  invitro_count="$(count_valid_seed_dirs "${INVITRO_RUN_DIR}" "${MULTI_WARMUP_LANDSCAPE_MAX_SEEDS}")"
  if (( invivo_count < 1 || invitro_count < 1 )); then
    echo "Parallel seed-space setup found invalid seed counts: invivo=${invivo_count}, invitro=${invitro_count}" >&2
    exit 1
  fi

  local max_seed_args=()
  if [[ -n "${MULTI_WARMUP_LANDSCAPE_MAX_SEEDS}" ]]; then
    max_seed_args+=("--max_seeds=${MULTI_WARMUP_LANDSCAPE_MAX_SEEDS}")
  fi

  local setup_invivo_cmd setup_invitro_cmd setup_wrap
  setup_invivo_cmd="$(shell_join \
    Rscript "${SEED_PLAN_SCRIPT}" \
    "--stage=build_seed_space_tasks" \
    "--project_root=${PROJECT_ROOT}" \
    "--dataset=invivo" \
    "--input_dir=${INVIVO_RUN_DIR}" \
    "--out_dir=${seed_space_root}" \
    "--tables_dir=${seed_table_dir}" \
    "${max_seed_args[@]}")"
  setup_invitro_cmd="$(shell_join \
    Rscript "${SEED_PLAN_SCRIPT}" \
    "--stage=build_seed_space_tasks" \
    "--project_root=${PROJECT_ROOT}" \
    "--dataset=invitro" \
    "--input_dir=${INVITRO_RUN_DIR}" \
    "--out_dir=${seed_space_root}" \
    "--tables_dir=${seed_table_dir}" \
    "--parameter_table_fallback=${PARAMETER_TABLE}" \
    "${max_seed_args[@]}")"
  setup_wrap="$(wrap_with_r_body "${setup_invivo_cmd} && ${setup_invitro_cmd}")"

  local invivo_array_wrap invitro_array_wrap
  invivo_array_wrap="$(wrap_with_r_command \
    Rscript "${SEED_PLAN_SCRIPT}" \
    "--stage=run_seed_space_task" \
    "--tasks_tsv=${invivo_tasks}")"
  invitro_array_wrap="$(wrap_with_r_command \
    Rscript "${SEED_PLAN_SCRIPT}" \
    "--stage=run_seed_space_task" \
    "--tasks_tsv=${invitro_tasks}")"

  local collect_invivo_cmd collect_invitro_cmd collect_wrap
  collect_invivo_cmd="$(shell_join \
    Rscript "${SEED_PLAN_SCRIPT}" \
    "--stage=collect_seed_space_tables" \
    "--dataset=invivo" \
    "--tasks_tsv=${invivo_tasks}" \
    "--tables_dir=${seed_table_dir}")"
  collect_invitro_cmd="$(shell_join \
    Rscript "${SEED_PLAN_SCRIPT}" \
    "--stage=collect_seed_space_tables" \
    "--dataset=invitro" \
    "--tasks_tsv=${invitro_tasks}" \
    "--tables_dir=${seed_table_dir}")"
  collect_wrap="$(wrap_with_r_body "${collect_invivo_cmd} && ${collect_invitro_cmd}")"

  local prepare_args=(
    bash "${SELF_SCRIPT}"
    "${PASSTHROUGH_ARGS[@]}"
    "--multi_warmup_invivo_best_csv=${invivo_best_csv}"
    "--multi_warmup_invivo_initial_csv=${invivo_initial_csv}"
    "--multi_warmup_invitro_best_csv=${invitro_best_csv}"
    "--multi_warmup_invitro_initial_csv=${invitro_initial_csv}"
    "--internal_stage=seed_plan"
    "--multi_warmup_seed_plan_as_job=FALSE"
    "--dry_run=FALSE"
  )
  local controller_args=(
    bash "${SELF_SCRIPT}"
    "${PASSTHROUGH_ARGS[@]}"
    "--internal_stage=curve_class_workflow"
    "--multi_warmup_seed_plan_as_job=FALSE"
    "--dry_run=FALSE"
  )
  local prepare_wrap controller_wrap
  prepare_wrap="$(shell_join bash -lc "$(shell_join "${prepare_args[@]}")")"
  controller_wrap="$(shell_join bash -lc "$(shell_join "${controller_args[@]}")")"

  printf "stage\tjob_id\tdependency\tqos\twalltime\tmem\tcpus\n" > "${CONTROLLER_JOBS_TSV}"

  local setup_job_id setup_dependency_id
  setup_job_id="$(submit_or_print "Submit seed-space setup" \
    sbatch \
    "--parsable" \
    "--job-name=o2mw_seedspace_setup" \
    "--cpus-per-task=1" \
    "--mem=${MULTI_WARMUP_SEED_PLAN_MEM}" \
    "--qos=${MULTI_WARMUP_SEED_PLAN_QOS}" \
    "--time=${MULTI_WARMUP_SUBMIT_TIME_LIMIT}" \
    "--output=${LOG_ROOT}/o2mw_seedspace_setup_%j.out" \
    "--error=${LOG_ROOT}/o2mw_seedspace_setup_%j.err" \
    "--wrap=${setup_wrap}")"
  setup_dependency_id="${setup_job_id%%;*}"
  printf "seed_space_setup\t%s\t\t%s\t%s\t%s\t1\n" \
    "${setup_job_id}" "${MULTI_WARMUP_SEED_PLAN_QOS}" "${MULTI_WARMUP_SUBMIT_TIME_LIMIT}" "${MULTI_WARMUP_SEED_PLAN_MEM}" >> "${CONTROLLER_JOBS_TSV}"

  local invivo_array_job_id invitro_array_job_id invivo_array_dependency_id invitro_array_dependency_id
  invivo_array_job_id="$(submit_or_print "Submit in vivo seed-space array" \
    sbatch \
    "--parsable" \
    "--job-name=o2mw_seedspace_vi" \
    "--dependency=afterok:${setup_dependency_id}" \
    "--array=1-${invivo_count}" \
    "--cpus-per-task=${MULTI_WARMUP_SEED_SPACE_CPUS}" \
    "--mem=${MULTI_WARMUP_SEED_SPACE_MEM}" \
    "--qos=${MULTI_WARMUP_SEED_SPACE_QOS}" \
    "--time=${MULTI_WARMUP_SEED_SPACE_TIME_LIMIT}" \
    "--output=${LOG_ROOT}/o2mw_seedspace_vi_%A_%a.out" \
    "--error=${LOG_ROOT}/o2mw_seedspace_vi_%A_%a.err" \
    "--wrap=${invivo_array_wrap}")"
  invivo_array_dependency_id="${invivo_array_job_id%%;*}"
  printf "seed_space_invivo_array\t%s\tafterok:%s\t%s\t%s\t%s\t%s\n" \
    "${invivo_array_job_id}" "${setup_dependency_id}" "${MULTI_WARMUP_SEED_SPACE_QOS}" "${MULTI_WARMUP_SEED_SPACE_TIME_LIMIT}" "${MULTI_WARMUP_SEED_SPACE_MEM}" "${MULTI_WARMUP_SEED_SPACE_CPUS}" >> "${CONTROLLER_JOBS_TSV}"

  invitro_array_job_id="$(submit_or_print "Submit in vitro seed-space array" \
    sbatch \
    "--parsable" \
    "--job-name=o2mw_seedspace_vt" \
    "--dependency=afterok:${setup_dependency_id}" \
    "--array=1-${invitro_count}" \
    "--cpus-per-task=${MULTI_WARMUP_SEED_SPACE_CPUS}" \
    "--mem=${MULTI_WARMUP_SEED_SPACE_MEM}" \
    "--qos=${MULTI_WARMUP_SEED_SPACE_QOS}" \
    "--time=${MULTI_WARMUP_SEED_SPACE_TIME_LIMIT}" \
    "--output=${LOG_ROOT}/o2mw_seedspace_vt_%A_%a.out" \
    "--error=${LOG_ROOT}/o2mw_seedspace_vt_%A_%a.err" \
    "--wrap=${invitro_array_wrap}")"
  invitro_array_dependency_id="${invitro_array_job_id%%;*}"
  printf "seed_space_invitro_array\t%s\tafterok:%s\t%s\t%s\t%s\t%s\n" \
    "${invitro_array_job_id}" "${setup_dependency_id}" "${MULTI_WARMUP_SEED_SPACE_QOS}" "${MULTI_WARMUP_SEED_SPACE_TIME_LIMIT}" "${MULTI_WARMUP_SEED_SPACE_MEM}" "${MULTI_WARMUP_SEED_SPACE_CPUS}" >> "${CONTROLLER_JOBS_TSV}"

  local collect_job_id collect_dependency_id
  collect_job_id="$(submit_or_print "Submit seed-space collect" \
    sbatch \
    "--parsable" \
    "--job-name=o2mw_seedspace_collect" \
    "--dependency=afterok:${invivo_array_dependency_id}:${invitro_array_dependency_id}" \
    "--cpus-per-task=1" \
    "--mem=${MULTI_WARMUP_SEED_PLAN_MEM}" \
    "--qos=${MULTI_WARMUP_SEED_PLAN_QOS}" \
    "--time=${MULTI_WARMUP_SUBMIT_TIME_LIMIT}" \
    "--output=${LOG_ROOT}/o2mw_seedspace_collect_%j.out" \
    "--error=${LOG_ROOT}/o2mw_seedspace_collect_%j.err" \
    "--wrap=${collect_wrap}")"
  collect_dependency_id="${collect_job_id%%;*}"
  printf "seed_space_collect\t%s\tafterok:%s:%s\t%s\t%s\t%s\t1\n" \
    "${collect_job_id}" "${invivo_array_dependency_id}" "${invitro_array_dependency_id}" "${MULTI_WARMUP_SEED_PLAN_QOS}" "${MULTI_WARMUP_SUBMIT_TIME_LIMIT}" "${MULTI_WARMUP_SEED_PLAN_MEM}" >> "${CONTROLLER_JOBS_TSV}"

  local prepare_job_id prepare_dependency_id
  prepare_job_id="$(submit_or_print "Submit landscape prepare" \
    sbatch \
    "--parsable" \
    "--job-name=o2mw_landscape_prepare" \
    "--dependency=afterok:${collect_dependency_id}" \
    "--cpus-per-task=${MULTI_WARMUP_SEED_PLAN_CPUS}" \
    "--mem=${MULTI_WARMUP_SEED_PLAN_MEM}" \
    "--qos=${MULTI_WARMUP_SEED_PLAN_QOS}" \
    "--time=${MULTI_WARMUP_SEED_PLAN_TIME_LIMIT}" \
    "--output=${LOG_ROOT}/o2mw_landscape_prepare_%j.out" \
    "--error=${LOG_ROOT}/o2mw_landscape_prepare_%j.err" \
    "--wrap=${prepare_wrap}")"
  prepare_dependency_id="${prepare_job_id%%;*}"
  printf "landscape_prepare\t%s\tafterok:%s\t%s\t%s\t%s\t%s\n" \
    "${prepare_job_id}" "${collect_dependency_id}" "${MULTI_WARMUP_SEED_PLAN_QOS}" "${MULTI_WARMUP_SEED_PLAN_TIME_LIMIT}" "${MULTI_WARMUP_SEED_PLAN_MEM}" "${MULTI_WARMUP_SEED_PLAN_CPUS}" >> "${CONTROLLER_JOBS_TSV}"

  local controller_job_id
  controller_job_id="$(submit_or_print "Submit dense-grid/finalize controller" \
    sbatch \
    "--parsable" \
    "--job-name=o2mw_curve_controller" \
    "--dependency=afterok:${prepare_dependency_id}" \
    "--cpus-per-task=1" \
    "--mem=${MULTI_WARMUP_SUBMIT_MEM}" \
    "--qos=${MULTI_WARMUP_SUBMIT_QOS}" \
    "--time=${MULTI_WARMUP_SUBMIT_TIME_LIMIT}" \
    "--output=${LOG_ROOT}/o2mw_curve_controller_%j.out" \
    "--error=${LOG_ROOT}/o2mw_curve_controller_%j.err" \
    "--wrap=${controller_wrap}")"
  printf "dense_grid_finalize_controller\t%s\tafterok:%s\t%s\t%s\t%s\t1\n" \
    "${controller_job_id}" "${prepare_dependency_id}" "${MULTI_WARMUP_SUBMIT_QOS}" "${MULTI_WARMUP_SUBMIT_TIME_LIMIT}" "${MULTI_WARMUP_SUBMIT_MEM}" >> "${CONTROLLER_JOBS_TSV}"

  log_msg "stage=parallel_seed_space_submitted setup=${setup_job_id} invivo_array=${invivo_array_job_id} invitro_array=${invitro_array_job_id} collect=${collect_job_id} prepare=${prepare_job_id} controller=${controller_job_id}"
  echo "Submitted seed-space setup job: ${setup_job_id}"
  echo "Submitted in vivo seed-space array: ${invivo_array_job_id} (1-${invivo_count})"
  echo "Submitted in vitro seed-space array: ${invitro_array_job_id} (1-${invitro_count})"
  echo "Submitted seed-space collect job: ${collect_job_id}"
  echo "Submitted dependent landscape prepare job: ${prepare_job_id}"
  echo "Submitted dependent dense-grid/finalize controller job: ${controller_job_id}"
  echo "Controller jobs table: ${CONTROLLER_JOBS_TSV}"
}

submit_seed_plan_workflow() {
  local seed_plan_args=(
    bash "${SELF_SCRIPT}"
    "${PASSTHROUGH_ARGS[@]}"
    "--internal_stage=seed_plan"
    "--multi_warmup_seed_plan_as_job=FALSE"
    "--dry_run=FALSE"
  )
  local submit_pairs_args=(
    bash "${SELF_SCRIPT}"
    "${PASSTHROUGH_ARGS[@]}"
    "--internal_stage=submit_pairs"
    "--multi_warmup_seed_plan_as_job=FALSE"
    "--dry_run=FALSE"
  )
  local seed_plan_wrap
  local submit_pairs_wrap
  seed_plan_wrap="$(shell_join bash -lc "$(shell_join "${seed_plan_args[@]}")")"
  submit_pairs_wrap="$(shell_join bash -lc "$(shell_join "${submit_pairs_args[@]}")")"

  printf "stage\tjob_id\tdependency\tqos\twalltime\tmem\tcpus\n" > "${CONTROLLER_JOBS_TSV}"
  local seed_plan_job_id
  seed_plan_job_id="$(submit_or_print "Submit multi-warm-up seed plan" \
    sbatch \
    "--parsable" \
    "--job-name=o2mw_seed_plan" \
    "--cpus-per-task=${MULTI_WARMUP_SEED_PLAN_CPUS}" \
    "--mem=${MULTI_WARMUP_SEED_PLAN_MEM}" \
    "--qos=${MULTI_WARMUP_SEED_PLAN_QOS}" \
    "--time=${MULTI_WARMUP_SEED_PLAN_TIME_LIMIT}" \
    "--output=${LOG_ROOT}/o2mw_seed_plan_%j.out" \
    "--error=${LOG_ROOT}/o2mw_seed_plan_%j.err" \
    "--wrap=${seed_plan_wrap}")"
  local seed_plan_dependency_id="${seed_plan_job_id%%;*}"
  printf "seed_plan\t%s\t\t%s\t%s\t%s\t%s\n" \
    "${seed_plan_job_id}" \
    "${MULTI_WARMUP_SEED_PLAN_QOS}" \
    "${MULTI_WARMUP_SEED_PLAN_TIME_LIMIT}" \
    "${MULTI_WARMUP_SEED_PLAN_MEM}" \
    "${MULTI_WARMUP_SEED_PLAN_CPUS}" >> "${CONTROLLER_JOBS_TSV}"

  local submit_pairs_job_id
  submit_pairs_job_id="$(submit_or_print "Submit multi-warm-up pair submission controller" \
    sbatch \
    "--parsable" \
    "--job-name=o2mw_submit_pairs" \
    "--dependency=afterok:${seed_plan_dependency_id}" \
    "--cpus-per-task=1" \
    "--mem=${MULTI_WARMUP_SUBMIT_MEM}" \
    "--qos=${MULTI_WARMUP_SUBMIT_QOS}" \
    "--time=${MULTI_WARMUP_SUBMIT_TIME_LIMIT}" \
    "--output=${LOG_ROOT}/o2mw_submit_pairs_%j.out" \
    "--error=${LOG_ROOT}/o2mw_submit_pairs_%j.err" \
    "--wrap=${submit_pairs_wrap}")"
  printf "submit_pairs\t%s\tafterok:%s\t%s\t%s\t%s\t1\n" \
    "${submit_pairs_job_id}" \
    "${seed_plan_dependency_id}" \
    "${MULTI_WARMUP_SUBMIT_QOS}" \
    "${MULTI_WARMUP_SUBMIT_TIME_LIMIT}" \
    "${MULTI_WARMUP_SUBMIT_MEM}" >> "${CONTROLLER_JOBS_TSV}"

  log_msg "stage=seed_plan_submitted job=${seed_plan_job_id} qos=${MULTI_WARMUP_SEED_PLAN_QOS} time=${MULTI_WARMUP_SEED_PLAN_TIME_LIMIT} mem=${MULTI_WARMUP_SEED_PLAN_MEM}"
  log_msg "stage=submit_pairs_controller_submitted job=${submit_pairs_job_id} dependency=afterok:${seed_plan_dependency_id}"
  echo "Submitted multi-warm-up seed-plan job: ${seed_plan_job_id}"
  echo "Submitted dependent pair-submission controller job: ${submit_pairs_job_id}"
  echo "Controller jobs table: ${CONTROLLER_JOBS_TSV}"
}

submit_dense_grid_curve_class_workflow() {
  if truthy "${MULTI_WARMUP_PARALLEL_SEED_SPACE}"; then
    submit_parallel_seed_space_landscape_workflow
    return
  fi

  local prepare_args=(
    bash "${SELF_SCRIPT}"
    "${PASSTHROUGH_ARGS[@]}"
    "--internal_stage=seed_plan"
    "--multi_warmup_seed_plan_as_job=FALSE"
    "--dry_run=FALSE"
  )
  local controller_args=(
    bash "${SELF_SCRIPT}"
    "${PASSTHROUGH_ARGS[@]}"
    "--internal_stage=curve_class_workflow"
    "--multi_warmup_seed_plan_as_job=FALSE"
    "--dry_run=FALSE"
  )
  local prepare_wrap controller_wrap
  prepare_wrap="$(shell_join bash -lc "$(shell_join "${prepare_args[@]}")")"
  controller_wrap="$(shell_join bash -lc "$(shell_join "${controller_args[@]}")")"

  printf "stage\tjob_id\tdependency\tqos\twalltime\tmem\tcpus\n" > "${CONTROLLER_JOBS_TSV}"
  local prepare_job_id
  prepare_job_id="$(submit_or_print "Submit landscape prepare" \
    sbatch \
    "--parsable" \
    "--job-name=o2mw_landscape_prepare" \
    "--cpus-per-task=${MULTI_WARMUP_SEED_PLAN_CPUS}" \
    "--mem=${MULTI_WARMUP_SEED_PLAN_MEM}" \
    "--qos=${MULTI_WARMUP_SEED_PLAN_QOS}" \
    "--time=${MULTI_WARMUP_SEED_PLAN_TIME_LIMIT}" \
    "--output=${LOG_ROOT}/o2mw_landscape_prepare_%j.out" \
    "--error=${LOG_ROOT}/o2mw_landscape_prepare_%j.err" \
    "--wrap=${prepare_wrap}")"
  local prepare_dependency_id="${prepare_job_id%%;*}"
  printf "landscape_prepare\t%s\t\t%s\t%s\t%s\t%s\n" \
    "${prepare_job_id}" \
    "${MULTI_WARMUP_SEED_PLAN_QOS}" \
    "${MULTI_WARMUP_SEED_PLAN_TIME_LIMIT}" \
    "${MULTI_WARMUP_SEED_PLAN_MEM}" \
    "${MULTI_WARMUP_SEED_PLAN_CPUS}" >> "${CONTROLLER_JOBS_TSV}"

  local controller_job_id
  controller_job_id="$(submit_or_print "Submit dense-grid/finalize controller" \
    sbatch \
    "--parsable" \
    "--job-name=o2mw_curve_controller" \
    "--dependency=afterok:${prepare_dependency_id}" \
    "--cpus-per-task=1" \
    "--mem=${MULTI_WARMUP_SUBMIT_MEM}" \
    "--qos=${MULTI_WARMUP_SUBMIT_QOS}" \
    "--time=${MULTI_WARMUP_SUBMIT_TIME_LIMIT}" \
    "--output=${LOG_ROOT}/o2mw_curve_controller_%j.out" \
    "--error=${LOG_ROOT}/o2mw_curve_controller_%j.err" \
    "--wrap=${controller_wrap}")"
  printf "dense_grid_finalize_controller\t%s\tafterok:%s\t%s\t%s\t%s\t1\n" \
    "${controller_job_id}" \
    "${prepare_dependency_id}" \
    "${MULTI_WARMUP_SUBMIT_QOS}" \
    "${MULTI_WARMUP_SUBMIT_TIME_LIMIT}" \
    "${MULTI_WARMUP_SUBMIT_MEM}" >> "${CONTROLLER_JOBS_TSV}"

  log_msg "stage=landscape_prepare_submitted job=${prepare_job_id}"
  log_msg "stage=dense_grid_finalize_controller_submitted job=${controller_job_id} dependency=afterok:${prepare_dependency_id}"
  echo "Submitted landscape prepare job: ${prepare_job_id}"
  echo "Submitted dependent dense-grid/finalize controller job: ${controller_job_id}"
  echo "Controller jobs table: ${CONTROLLER_JOBS_TSV}"
}

run_curve_class_workflow() {
  local dense_out="${MULTI_WARMUP_ROOT}/cross_validation/best_fit_parameter_feature/03_dense-grid_monotonicity_classification/monotonicity_classification/dense-grid_monotonicity_classification"
  local max_seed_args=()
  if [[ -n "${MULTI_WARMUP_LANDSCAPE_MAX_SEEDS}" ]]; then
    max_seed_args+=("--max_seeds=${MULTI_WARMUP_LANDSCAPE_MAX_SEEDS}")
  fi

  if [[ ! -s "${CONTROLLER_JOBS_TSV}" ]]; then
    printf "stage\tjob_id\tdependency\tqos\twalltime\tmem\tcpus\n" > "${CONTROLLER_JOBS_TSV}"
  fi
  log_msg "stage=curve_class_workflow_start dense_out=${dense_out}"
  run_or_print "Build dense-grid monotonicity tasks" \
    Rscript "${DENSE_GRID_SCRIPT}" \
    "--mode=build_tasks" \
    "--part=monotonicity" \
    "--run_dir=${INVIVO_RUN_DIR}" \
    "--out_dir=${dense_out}" \
    "--tasks_per_array_task=${MULTI_WARMUP_MONOTONICITY_TASKS_PER_ARRAY_TASK}" \
    "${max_seed_args[@]}"

  local metadata="${dense_out}/hpc/task_lists/monotonicity_task_metadata.tsv"
  local task_file="${dense_out}/hpc/task_lists/monotonicity_seed_o2_tasks.tsv"
  if [[ ! -f "${metadata}" ]]; then
    echo "Missing dense-grid task metadata: ${metadata}" >&2
    exit 1
  fi
  local n_array_tasks
  n_array_tasks="$(awk -F $'\t' '$1 == "n_array_tasks" {print $2; exit}' "${metadata}")"
  if ! [[ "${n_array_tasks}" =~ ^[0-9]+$ ]] || (( n_array_tasks < 1 )); then
    echo "Invalid n_array_tasks in ${metadata}: ${n_array_tasks}" >&2
    exit 1
  fi

  local array_wrap
  array_wrap="$(wrap_with_r_command \
    Rscript "${DENSE_GRID_SCRIPT}" \
    "--mode=run_tasks" \
    "--part=monotonicity" \
    "--run_dir=${INVIVO_RUN_DIR}" \
    "--out_dir=${dense_out}" \
    "--task_file=${task_file}" \
    "--tasks_per_array_task=${MULTI_WARMUP_MONOTONICITY_TASKS_PER_ARRAY_TASK}")"
  local array_job_id
  array_job_id="$(submit_or_print "Submit dense-grid monotonicity array" \
    sbatch \
    "--parsable" \
    "--job-name=o2mw_curve_grid" \
    "--array=1-${n_array_tasks}" \
    "--cpus-per-task=${MULTI_WARMUP_MONOTONICITY_CPUS}" \
    "--mem=${MULTI_WARMUP_MONOTONICITY_MEM}" \
    "--qos=${MULTI_WARMUP_MONOTONICITY_QOS}" \
    "--time=${MULTI_WARMUP_MONOTONICITY_TIME_LIMIT}" \
    "--output=${LOG_ROOT}/o2mw_curve_grid_%A_%a.out" \
    "--error=${LOG_ROOT}/o2mw_curve_grid_%A_%a.err" \
    "--wrap=${array_wrap}")"
  local array_dependency_id="${array_job_id%%;*}"
  printf "monotonicity_array\t%s\t\t%s\t%s\t%s\t%s\n" \
    "${array_job_id}" "${MULTI_WARMUP_MONOTONICITY_QOS}" "${MULTI_WARMUP_MONOTONICITY_TIME_LIMIT}" "${MULTI_WARMUP_MONOTONICITY_MEM}" "${MULTI_WARMUP_MONOTONICITY_CPUS}" >> "${CONTROLLER_JOBS_TSV}"

  local merge_wrap
  merge_wrap="$(wrap_with_r_command Rscript "${DENSE_GRID_SCRIPT}" "--mode=merge" "--part=monotonicity" "--run_dir=${INVIVO_RUN_DIR}" "--out_dir=${dense_out}" "${max_seed_args[@]}")"
  local merge_job_id
  merge_job_id="$(submit_or_print "Submit dense-grid monotonicity merge" \
    sbatch --parsable --job-name=o2mw_curve_merge "--dependency=afterok:${array_dependency_id}" \
    "--cpus-per-task=1" "--mem=${MULTI_WARMUP_VALIDATION_MEM}" "--qos=${MULTI_WARMUP_VALIDATION_QOS}" "--time=${MULTI_WARMUP_VALIDATION_TIME_LIMIT}" \
    "--output=${LOG_ROOT}/o2mw_curve_merge_%j.out" "--error=${LOG_ROOT}/o2mw_curve_merge_%j.err" "--wrap=${merge_wrap}")"
  local merge_dependency_id="${merge_job_id%%;*}"
  printf "monotonicity_merge\t%s\tafterok:%s\t%s\t%s\t%s\t1\n" "${merge_job_id}" "${array_dependency_id}" "${MULTI_WARMUP_VALIDATION_QOS}" "${MULTI_WARMUP_VALIDATION_TIME_LIMIT}" "${MULTI_WARMUP_VALIDATION_MEM}" >> "${CONTROLLER_JOBS_TSV}"

  local finalize_args=(
    Rscript "${SEED_PLAN_SCRIPT}"
    "--stage=finalize_pairs"
    "--project_root=${PROJECT_ROOT}"
    "--invivo_run_dir=${INVIVO_RUN_DIR}"
    "--invitro_run_dir=${INVITRO_RUN_DIR}"
    "--out_dir=${MULTI_WARMUP_ROOT}"
    "--reductions=${MULTI_WARMUP_REDUCTIONS}"
    "--umap_seed=${MULTI_WARMUP_LANDSCAPE_UMAP_SEED}"
    "--tsne_seed=${MULTI_WARMUP_TSNE_SEED}"
    "--cluster_seed=${MULTI_WARMUP_CLUSTER_SEED}"
    "--subcluster_seed=${MULTI_WARMUP_SUBCLUSTER_SEED}"
    "--pairing_policy=${MULTI_WARMUP_PAIRING_POLICY}"
    "--deduplicate_pairs=${MULTI_WARMUP_DEDUPLICATE_PAIRS}"
    "--invivo_curve_filter=${MULTI_WARMUP_INVIVO_CURVE_FILTER}"
    "--invivo_curve_class=${MULTI_WARMUP_INVIVO_CURVE_CLASS}"
  )
  if [[ -n "${MULTI_WARMUP_REFERENCE_SUBCLUSTER_DIR}" ]]; then finalize_args+=("--reference_subcluster_dir=${MULTI_WARMUP_REFERENCE_SUBCLUSTER_DIR}"); fi
  local finalize_wrap
  finalize_wrap="$(wrap_with_r_command "${finalize_args[@]}")"
  local finalize_job_id
  finalize_job_id="$(submit_or_print "Submit final seed-pair selection" \
    sbatch --parsable --job-name=o2mw_finalize_pairs "--dependency=afterok:${merge_dependency_id}" \
    "--cpus-per-task=1" "--mem=${MULTI_WARMUP_VALIDATION_MEM}" "--qos=${MULTI_WARMUP_VALIDATION_QOS}" "--time=${MULTI_WARMUP_VALIDATION_TIME_LIMIT}" \
    "--output=${LOG_ROOT}/o2mw_finalize_pairs_%j.out" "--error=${LOG_ROOT}/o2mw_finalize_pairs_%j.err" "--wrap=${finalize_wrap}")"
  local finalize_dependency_id="${finalize_job_id%%;*}"
  printf "finalize_pairs\t%s\tafterok:%s\t%s\t%s\t%s\t1\n" "${finalize_job_id}" "${merge_dependency_id}" "${MULTI_WARMUP_VALIDATION_QOS}" "${MULTI_WARMUP_VALIDATION_TIME_LIMIT}" "${MULTI_WARMUP_VALIDATION_MEM}" >> "${CONTROLLER_JOBS_TSV}"

  local submit_pairs_args=(
    bash "${SELF_SCRIPT}"
    "${PASSTHROUGH_ARGS[@]}"
    "--internal_stage=submit_pairs"
    "--multi_warmup_seed_plan_as_job=FALSE"
    "--dry_run=FALSE"
  )
  local submit_pairs_wrap submit_pairs_job_id
  submit_pairs_wrap="$(shell_join bash -lc "$(shell_join "${submit_pairs_args[@]}")")"
  submit_pairs_job_id="$(submit_or_print "Submit multi-warm-up pair submission controller" \
    sbatch --parsable --job-name=o2mw_submit_pairs "--dependency=afterok:${finalize_dependency_id}" \
    "--cpus-per-task=1" "--mem=${MULTI_WARMUP_SUBMIT_MEM}" "--qos=${MULTI_WARMUP_SUBMIT_QOS}" "--time=${MULTI_WARMUP_SUBMIT_TIME_LIMIT}" \
    "--output=${LOG_ROOT}/o2mw_submit_pairs_%j.out" "--error=${LOG_ROOT}/o2mw_submit_pairs_%j.err" "--wrap=${submit_pairs_wrap}")"
  printf "submit_pairs\t%s\tafterok:%s\t%s\t%s\t%s\t1\n" "${submit_pairs_job_id}" "${finalize_dependency_id}" "${MULTI_WARMUP_SUBMIT_QOS}" "${MULTI_WARMUP_SUBMIT_TIME_LIMIT}" "${MULTI_WARMUP_SUBMIT_MEM}" >> "${CONTROLLER_JOBS_TSV}"
  log_msg "stage=curve_class_workflow_submitted array=${array_job_id} finalize=${finalize_job_id} submit_pairs=${submit_pairs_job_id}"
}

load_r_module
if [[ "${INTERNAL_STAGE}" == "submit_pairs" ]]; then
  log_msg "stage=submit_pairs_start pair_method=${MULTI_WARMUP_PAIR_METHOD} manifest=${MULTI_WARMUP_ROOT}/multi_warmup_manifest.tsv"
else
  log_msg "stage=seed_plan pair_method=${MULTI_WARMUP_PAIR_METHOD} invivo_run_dir=${INVIVO_RUN_DIR} invitro_run_dir=${INVITRO_RUN_DIR}"
fi
if [[ "${MULTI_WARMUP_PAIR_METHOD}" == "landscape_subcluster" ]]; then
  landscape_stage="prepare_landscape"
  effective_invivo_curve_filter="${MULTI_WARMUP_INVIVO_CURVE_FILTER}"
  seed_plan_cmd=(
    Rscript "${SEED_PLAN_SCRIPT}"
    "--stage=${landscape_stage}"
    "--project_root=${PROJECT_ROOT}"
    "--invivo_run_dir=${INVIVO_RUN_DIR}"
    "--invitro_run_dir=${INVITRO_RUN_DIR}"
    "--out_dir=${MULTI_WARMUP_ROOT}"
    "--reductions=${MULTI_WARMUP_REDUCTIONS}"
    "--umap_seed=${MULTI_WARMUP_LANDSCAPE_UMAP_SEED}"
    "--tsne_seed=${MULTI_WARMUP_TSNE_SEED}"
    "--cluster_seed=${MULTI_WARMUP_CLUSTER_SEED}"
    "--subcluster_seed=${MULTI_WARMUP_SUBCLUSTER_SEED}"
    "--pairing_policy=${MULTI_WARMUP_PAIRING_POLICY}"
    "--deduplicate_pairs=${MULTI_WARMUP_DEDUPLICATE_PAIRS}"
    "--invivo_curve_filter=${effective_invivo_curve_filter}"
    "--invivo_curve_class=${MULTI_WARMUP_INVIVO_CURVE_CLASS}"
  )
  if [[ -n "${MULTI_WARMUP_LANDSCAPE_MAX_SEEDS}" ]]; then seed_plan_cmd+=("--max_seeds=${MULTI_WARMUP_LANDSCAPE_MAX_SEEDS}"); fi
  if [[ -n "${MULTI_WARMUP_REFERENCE_SUBCLUSTER_DIR}" ]]; then seed_plan_cmd+=("--reference_subcluster_dir=${MULTI_WARMUP_REFERENCE_SUBCLUSTER_DIR}"); fi
  if [[ -n "${MULTI_WARMUP_INVIVO_BEST_CSV}" ]]; then seed_plan_cmd+=("--invivo_best_csv=${MULTI_WARMUP_INVIVO_BEST_CSV}"); fi
  if [[ -n "${MULTI_WARMUP_INVIVO_INITIAL_CSV}" ]]; then seed_plan_cmd+=("--invivo_initial_csv=${MULTI_WARMUP_INVIVO_INITIAL_CSV}"); fi
  if [[ -n "${MULTI_WARMUP_INVITRO_BEST_CSV}" ]]; then seed_plan_cmd+=("--invitro_best_csv=${MULTI_WARMUP_INVITRO_BEST_CSV}"); fi
  if [[ -n "${MULTI_WARMUP_INVITRO_INITIAL_CSV}" ]]; then seed_plan_cmd+=("--invitro_initial_csv=${MULTI_WARMUP_INVITRO_INITIAL_CSV}"); fi
else
  seed_plan_cmd=(
    Rscript "${SEED_PLAN_SCRIPT}"
    "--invivo_run_dir=${INVIVO_RUN_DIR}"
    "--invitro_run_dir=${INVITRO_RUN_DIR}"
    "--out_dir=${MULTI_WARMUP_ROOT}"
    "--top_n=${MULTI_WARMUP_TOP_N}"
    "--invivo_top_n=${MULTI_WARMUP_INVIVO_TOP_N}"
    "--invitro_top_n=${MULTI_WARMUP_INVITRO_TOP_N}"
    "--umap_seed=${MULTI_WARMUP_UMAP_SEED}"
    "--invivo_k=${MULTI_WARMUP_INVIVO_K}"
    "--invitro_k=${MULTI_WARMUP_INVITRO_K}"
    "--invitro_anchor_ranks=${MULTI_WARMUP_INVITRO_ANCHOR_RANKS}"
    "--include_phase2=${MULTI_WARMUP_INCLUDE_PHASE2}"
    "--phase2_invitro_anchor_ranks=${MULTI_WARMUP_PHASE2_INVITRO_ANCHOR_RANKS}"
  )
fi

MANIFEST="${MULTI_WARMUP_ROOT}/multi_warmup_manifest.tsv"
PLAN_MODE_FILE="${MULTI_WARMUP_ROOT}/multi_warmup_seed_plan_mode.tsv"
if [[ "${INTERNAL_STAGE}" == "seed_plan" ]]; then
  run_or_print "Generate multi-warmup seed plan" "${seed_plan_cmd[@]}"
  PREPARED_SUMMARY="${MULTI_WARMUP_ROOT}/landscape_subcluster/pooled_invivo_invitro/full_data_in_vivo_clustring/Tables/pooled_invivo_invitro_best_subcluster_summary_by_method.csv"
  if [[ "${MULTI_WARMUP_PAIR_METHOD}" == "landscape_subcluster" ]]; then
    if ! truthy "${DRY_RUN}" && [[ ! -f "${PREPARED_SUMMARY}" ]]; then
      echo "Missing prepared landscape subcluster summary: ${PREPARED_SUMMARY}" >&2
      exit 1
    fi
    log_msg "stage=landscape_prepare_done summary=${PREPARED_SUMMARY}"
    exit 0
  fi
  if ! truthy "${DRY_RUN}" && [[ ! -f "${MANIFEST}" ]]; then
    echo "Missing generated manifest: ${MANIFEST}" >&2
    exit 1
  fi
  log_msg "stage=seed_plan_done manifest=${MANIFEST}"
  exit 0
fi

if [[ "${INTERNAL_STAGE}" == "curve_class_workflow" ]]; then
  run_curve_class_workflow
  exit 0
fi

if [[ -z "${INTERNAL_STAGE}" ]] && truthy "${MULTI_WARMUP_SEED_PLAN_AS_JOB}"; then
  if [[ "${MULTI_WARMUP_PAIR_METHOD}" == "landscape_subcluster" ]]; then
    submit_dense_grid_curve_class_workflow
  else
    submit_seed_plan_workflow
  fi
  exit 0
fi

if [[ -z "${INTERNAL_STAGE}" ]]; then
  run_or_print "Generate multi-warmup seed plan" "${seed_plan_cmd[@]}"
  if [[ "${MULTI_WARMUP_PAIR_METHOD}" == "landscape_subcluster" ]]; then
    echo "landscape_subcluster warmup requires --multi_warmup_seed_plan_as_job=TRUE so dense-grid curve classification can run as dependent HPC jobs." >&2
    exit 2
  fi
else
  log_msg "stage=submit_pairs_after_seed_plan dependency_manifest=${MANIFEST}"
fi

if ! truthy "${DRY_RUN}" && [[ ! -f "${MANIFEST}" ]]; then
  echo "Missing generated manifest: ${MANIFEST}" >&2
  exit 1
fi
if truthy "${DRY_RUN}" && [[ ! -f "${MANIFEST}" ]]; then
  log_msg "stage=dry_run complete_after_seed_plan_command"
  exit 0
fi

printf "warmup_label\tpair_index\ttotal_pairs\tjoint_run_dir\tsubmit_status\tjob_id\tpostprocess_job_id\tarray_spec\tqos\twalltime\n" > "${JOBS_TSV}"
total_pairs=$(( $(wc -l < "${MANIFEST}") - 1 ))
if (( total_pairs < 1 )); then
  plan_mode=""
  if [[ -f "${PLAN_MODE_FILE}" ]]; then
    plan_mode="$(awk -F $'\t' '$1 == "mode" {print $2; exit}' "${PLAN_MODE_FILE}")"
  fi
  case "${plan_mode}" in
    invivo_only|invitro_only)
      log_msg "stage=cluster_only mode=${plan_mode} manifest=${MANIFEST}"
      exit 0
      ;;
    *)
      echo "Generated manifest has no warm-up pairs: ${MANIFEST}" >&2
      exit 1
      ;;
  esac
fi
pair_index=0
post_job_ids=()
POSTPROCESS_JOB_IDS_FILE="${MULTI_WARMUP_ROOT}/.postprocess_job_ids.tmp"
: > "${POSTPROCESS_JOB_IDS_FILE}"
log_msg "stage=submit_pairs total_pairs=${total_pairs} qos=${JOINT_QOS} time=${JOINT_TIME_LIMIT}"
mkdir -p "${MULTI_WARMUP_ROOT}/joint_soft_coupling_tables"

tail -n +2 "${MANIFEST}" | while IFS=$'\t' read -r warmup_label phase invivo_family invivo_rank invivo_seed invivo_seed_dir invitro_family invitro_rank invitro_seed invitro_seed_dir selection_reason joint_run_prefix joint_table; do
  pair_index=$((pair_index + 1))
  log_msg "stage=pair_submit_start pair=${pair_index}/${total_pairs} warmup_label=${warmup_label}"
  if [[ -z "${joint_table}" ]]; then
    echo "Manifest row for ${warmup_label} is missing joint_soft_coupling_parameters_table." >&2
    exit 1
  fi
  case "${joint_table}" in
    /*) ;;
    *) joint_table="${MULTI_WARMUP_ROOT}/${joint_table}" ;;
  esac
  mkdir -p "$(dirname "${joint_table}")"
  run_or_print "Generate pair soft-coupling table" \
    Rscript "${MAKE_TABLE_SCRIPT}" \
    "--invivo-seed-dir=${invivo_seed_dir}" \
    "--invitro-seed-dir=${invitro_seed_dir}" \
    "--seed-label=${warmup_label}" \
    "--out=${joint_table}" \
    "--delta-params=${JOINT_SOFT_COUPLING_DELTA_PARAMS}"
  if ! truthy "${DRY_RUN}" && [[ ! -f "${joint_table}" ]]; then
    echo "Pair soft-coupling table was not created for ${warmup_label}: ${joint_table}" >&2
    exit 1
  fi

  export_arg="ALL,PROJECT_ROOT=${PROJECT_ROOT},RUNNER_SCRIPT=${JOINT_RUNNER_SCRIPT},CONFIG_PATH=${CONFIG_PATH},OUT_ROOT=${MULTI_WARMUP_ROOT},RUN_PREFIX=${joint_run_prefix},TOTAL_SEEDS=${JOINT_TOTAL_SEEDS},ARRAY_TASKS=${JOINT_ARRAY_TASKS},SEEDS_PER_TASK=${JOINT_SEEDS_PER_TASK},N_CORES=${JOINT_N_CORES},R_MODULE=${R_MODULE},PARAMETER_TABLE=${PARAMETER_TABLE},FIT_OBJECTS_DIR=${FIT_OBJECTS_DIR},FLOW_DENSITY_PATH=${FLOW_DENSITY_PATH},DEATH_DATA_PATH=${DEATH_DATA_PATH},DEATH_WEIGHT=${DEATH_WEIGHT},SIGMA_DEATH_LOGIT=${SIGMA_DEATH_LOGIT},DEATH_FRACTION_EPS=${DEATH_FRACTION_EPS},ITERMAX=${ITERMAX},DE_RELTOL=${DE_RELTOL},DE_STEPTOL=${DE_STEPTOL},NP=${NP},AUTO_VIZ=${AUTO_VIZ},JOINT_WARMUP_ENABLE=TRUE,JOINT_WARMUP_SEED_LABEL=${warmup_label},JOINT_WARMUP_INVIVO_SEED_DIR=${invivo_seed_dir},JOINT_WARMUP_INVITRO_SEED_DIR=${invitro_seed_dir},JOINT_WARMUP_SIGMAN=${JOINT_WARMUP_SIGMAN},JOINT_SOFT_COUPLING_SIGMA_DEFAULT=${JOINT_SOFT_COUPLING_SIGMA_DEFAULT},JOINT_SOFT_COUPLING_WELSCH_C=${JOINT_SOFT_COUPLING_WELSCH_C},JOINT_SOFT_COUPLING_PARAMETERS_TABLE=${joint_table}"
  job_id="$(submit_or_print "Submit pair joint array ${warmup_label}" \
    sbatch \
    "--job-name=o2mw_${warmup_label}" \
    "--array=1-${JOINT_ARRAY_TASKS}" \
    "--cpus-per-task=${JOINT_N_CORES}" \
    "--mem=${JOINT_MEM}" \
    "--qos=${JOINT_QOS}" \
    "--time=${JOINT_TIME_LIMIT}" \
    "--output=${LOG_ROOT}/o2mw_${warmup_label}_%A_%a.out" \
    "--error=${LOG_ROOT}/o2mw_${warmup_label}_%A_%a.err" \
    "--export=${export_arg}" \
    "${JOINT_ARRAY_SCRIPT}")"

  post_id="$(submit_or_print "Submit pair extra_results ${warmup_label}" \
    sbatch \
    "--job-name=o2mw_er_${warmup_label}" \
    "--dependency=afterok:${job_id}" \
    "--cpus-per-task=1" \
    "--mem=${POSTPROCESS_MEM}" \
    "--qos=${POSTPROCESS_QOS}" \
    "--time=${POSTPROCESS_TIME_LIMIT}" \
    "--output=${LOG_ROOT}/o2mw_er_${warmup_label}_%j.out" \
    "--error=${LOG_ROOT}/o2mw_er_${warmup_label}_%j.err" \
    "--export=ALL,PROJECT_ROOT=${PROJECT_ROOT},RUN_DIR=${MULTI_WARMUP_ROOT}/${joint_run_prefix},R_MODULE=${R_MODULE},FORCE_EXTRA_RESULTS=${FORCE_EXTRA_RESULTS}" \
    "${POSTPROCESS_SCRIPT}")"
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "${warmup_label}" "${pair_index}" "${total_pairs}" "${MULTI_WARMUP_ROOT}/${joint_run_prefix}" "submitted" "${job_id}" "${post_id}" "1-${JOINT_ARRAY_TASKS}" "${JOINT_QOS}" "${JOINT_TIME_LIMIT}" >> "${JOBS_TSV}"
  echo "${post_id}" >> "${POSTPROCESS_JOB_IDS_FILE}"
  log_msg "stage=pair_submit_done pair=${pair_index}/${total_pairs} warmup_label=${warmup_label} array_job=${job_id} postprocess_job=${post_id}"
done

dependency_ids="$(paste -sd: "${POSTPROCESS_JOB_IDS_FILE}")"
wrap_cmd="if command -v ml >/dev/null 2>&1; then ml ${R_MODULE}; elif command -v module >/dev/null 2>&1; then module load ${R_MODULE}; fi; cd ${PROJECT_ROOT} && Rscript ${COLLECT_SCRIPT} --multi_warmup_root=${MULTI_WARMUP_ROOT} --manifest=${MANIFEST} --out_dir=${MULTI_WARMUP_ROOT} && Rscript ${REPORT_SCRIPT} --multi_warmup_root=${MULTI_WARMUP_ROOT} --out=${MULTI_WARMUP_ROOT}/multi-warm-up_results.html"
report_job_id="$(submit_or_print "Submit final multi-warmup report" \
  sbatch \
  "--job-name=o2mw_report" \
  "--dependency=afterok:${dependency_ids}" \
  "--cpus-per-task=1" \
  "--mem=${REPORT_MEM}" \
  "--qos=${REPORT_QOS}" \
  "--time=${REPORT_TIME_LIMIT}" \
  "--output=${LOG_ROOT}/o2mw_report_%j.out" \
  "--error=${LOG_ROOT}/o2mw_report_%j.err" \
  "--wrap=${wrap_cmd}")"
log_msg "stage=submitted_all final_report_job=${report_job_id} report=${MULTI_WARMUP_ROOT}/multi-warm-up_results.html"
