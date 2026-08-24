#!/usr/bin/env bash
# Run a multi-warm-up joint sweep from source in vivo/in vitro runs.

set -euo pipefail

O2SD_SHELL_UTILS="$(cd "$(dirname "${BASH_SOURCE[0]}")/../util" && pwd)/o2_supply_demand_map_shell_utils.sh"
# shellcheck source=../util/o2_supply_demand_map_shell_utils.sh
source "${O2SD_SHELL_UTILS}"

usage() {
  cat <<'EOF'
Usage:
  bash run_multi_warmup_joint.sh --invivo_run_dir=DIR --invitro_run_dir=DIR [options]

Required:
  --invivo_run_dir=DIR
  --invitro_run_dir=DIR

Common options:
  --project_root=DIR
  --config_path=FILE
  --out_root=DIR
  --multi_warmup_prefix=NAME
  --joint_seeds_csv=1,2,3
  --joint_n_cores=22
  --dry_run=TRUE|FALSE
  --run_extra_results=TRUE|FALSE
  --force_extra_results=TRUE|FALSE

Seed-plan options:
  --multi_warmup_pair_method=legacy|landscape_subcluster
  --multi_warmup_top_n=10
  --multi_warmup_invivo_top_n=10  (0 disables in vivo source clustering; not both sides)
  --multi_warmup_invitro_top_n=10 (0 disables in vitro source clustering; not both sides)
  --multi_warmup_umap_seed=1
  --multi_warmup_invivo_k=auto
  --multi_warmup_invitro_anchor_ranks=1
  --multi_warmup_include_phase2=FALSE
  --multi_warmup_phase2_invitro_anchor_ranks=auto
  --multi_warmup_reductions=tsne,umap
  --multi_warmup_landscape_umap_seed=123
  --multi_warmup_landscape_max_seeds=N
  --multi_warmup_pairing_policy=cartesian_by_method|invitro_best_to_invivo_subclusters
  --multi_warmup_deduplicate_pairs=FALSE
  --multi_warmup_reference_subcluster_dir=DIR
  --multi_warmup_invivo_curve_filter=TRUE|FALSE
  --multi_warmup_invivo_curve_class=monotone_increasing
  --joint_soft_coupling_sigma_default=0.65
  --joint_soft_coupling_welsch_c=0.4
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
      --joint_seeds_csv=*) JOINT_SEEDS_CSV="${arg#*=}" ;;
      --joint_n_cores=*|--n_cores=*) JOINT_N_CORES="${arg#*=}" ;;
      --r_module=*) R_MODULE="${arg#*=}" ;;
      --dry_run=*) DRY_RUN="${arg#*=}" ;;
      --run_extra_results=*) RUN_EXTRA_RESULTS="${arg#*=}" ;;
      --force_extra_results=*) FORCE_EXTRA_RESULTS="${arg#*=}" ;;
      --parameter_table=*|--invitro_parameter_table=*) PARAMETER_TABLE="${arg#*=}" ;;
      --fit_objects_dir=*) FIT_OBJECTS_DIR="${arg#*=}" ;;
      --flow_density_path=*) FLOW_DENSITY_PATH="${arg#*=}" ;;
      --itermax=*) ITERMAX="${arg#*=}" ;;
      --de_reltol=*) DE_RELTOL="${arg#*=}" ;;
      --de_steptol=*) DE_STEPTOL="${arg#*=}" ;;
      --np=*|--NP=*) NP="${arg#*=}" ;;
      --auto_viz=*) AUTO_VIZ="${arg#*=}" ;;
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
      --multi_warmup_reference_subcluster_dir=*|--reference_subcluster_dir=*) MULTI_WARMUP_REFERENCE_SUBCLUSTER_DIR="${arg#*=}" ;;
      *) echo "Unknown argument: ${arg}" >&2; usage >&2; exit 2 ;;
    esac
  done
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKFLOW_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
DEFAULT_PROJECT_ROOT="$(cd "${WORKFLOW_ROOT}/../../.." && pwd)"

PROJECT_ROOT="${PROJECT_ROOT:-${DEFAULT_PROJECT_ROOT}}"
CONFIG_PATH="${CONFIG_PATH:-}"
OUT_ROOT="${OUT_ROOT:-}"
MULTI_WARMUP_PREFIX="${MULTI_WARMUP_PREFIX:-fit_joint_multi_warmup}"
INVIVO_RUN_DIR="${INVIVO_RUN_DIR:-}"
INVITRO_RUN_DIR="${INVITRO_RUN_DIR:-}"
JOINT_SEEDS_CSV="${JOINT_SEEDS_CSV:-1}"
JOINT_N_CORES="${JOINT_N_CORES:-1}"
R_MODULE="${R_MODULE:-R/4.4}"
DRY_RUN="${DRY_RUN:-FALSE}"
RUN_EXTRA_RESULTS="${RUN_EXTRA_RESULTS:-TRUE}"
FORCE_EXTRA_RESULTS="${FORCE_EXTRA_RESULTS:-FALSE}"
PARAMETER_TABLE="${PARAMETER_TABLE:-}"
FIT_OBJECTS_DIR="${FIT_OBJECTS_DIR:-}"
FLOW_DENSITY_PATH="${FLOW_DENSITY_PATH:-}"
ITERMAX="${ITERMAX:-1000}"
DE_RELTOL="${DE_RELTOL:-1e-4}"
DE_STEPTOL="${DE_STEPTOL:-25}"
NP="${NP:-80}"
AUTO_VIZ="${AUTO_VIZ:-TRUE}"
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
MULTI_WARMUP_INVIVO_CURVE_FILTER="${MULTI_WARMUP_INVIVO_CURVE_FILTER:-FALSE}"
MULTI_WARMUP_INVIVO_CURVE_CLASS="${MULTI_WARMUP_INVIVO_CURVE_CLASS:-monotone_increasing}"
MULTI_WARMUP_REFERENCE_SUBCLUSTER_DIR="${MULTI_WARMUP_REFERENCE_SUBCLUSTER_DIR:-}"

parse_args "$@"

MULTI_WARMUP_PAIR_METHOD="$(echo "${MULTI_WARMUP_PAIR_METHOD}" | tr '[:upper:]' '[:lower:]' | tr '-' '_')"
case "${MULTI_WARMUP_PAIR_METHOD}" in
  legacy|landscape_subcluster) ;;
  *) echo "--multi_warmup_pair_method must be legacy or landscape_subcluster, got: ${MULTI_WARMUP_PAIR_METHOD}" >&2; exit 2 ;;
esac
require_nonnegative_int MULTI_WARMUP_INVIVO_TOP_N "${MULTI_WARMUP_INVIVO_TOP_N}"
require_nonnegative_int MULTI_WARMUP_INVITRO_TOP_N "${MULTI_WARMUP_INVITRO_TOP_N}"
if (( MULTI_WARMUP_INVIVO_TOP_N == 0 && MULTI_WARMUP_INVITRO_TOP_N == 0 )); then
  echo "At least one of --invivo_top_n or --invitro_top_n must be greater than 0." >&2
  exit 2
fi
if [[ "${MULTI_WARMUP_PAIR_METHOD}" == "landscape_subcluster" ]] && (( MULTI_WARMUP_INVIVO_TOP_N == 0 || MULTI_WARMUP_INVITRO_TOP_N == 0 )); then
  echo "landscape_subcluster pair method requires both in vivo and in vitro source runs." >&2
  exit 2
fi
if (( MULTI_WARMUP_INVIVO_TOP_N > 0 )) && is_null_value "${INVIVO_RUN_DIR}"; then
  echo "--invivo_run_dir is required when --invivo_top_n > 0." >&2
  exit 2
fi
if (( MULTI_WARMUP_INVITRO_TOP_N > 0 )) && is_null_value "${INVITRO_RUN_DIR}"; then
  echo "--invitro_run_dir is required when --invitro_top_n > 0." >&2
  exit 2
fi

PROJECT_ROOT="$(cd "${PROJECT_ROOT}" && pwd)"
if [[ -z "${CONFIG_PATH}" ]]; then CONFIG_PATH="${PROJECT_ROOT}/oxygen/config/O2_supply_demand.yaml"; fi
if [[ -z "${OUT_ROOT}" ]]; then OUT_ROOT="${PROJECT_ROOT}/oxygen/results"; fi
CONFIG_PATH="$(cd "$(dirname "${CONFIG_PATH}")" && pwd)/$(basename "${CONFIG_PATH}")"
OUT_ROOT="$(mkdir -p "${OUT_ROOT}" && cd "${OUT_ROOT}" && pwd)"
if (( MULTI_WARMUP_INVIVO_TOP_N > 0 )); then INVIVO_RUN_DIR="$(cd "${INVIVO_RUN_DIR}" && pwd)"; else INVIVO_RUN_DIR=""; fi
if (( MULTI_WARMUP_INVITRO_TOP_N > 0 )); then INVITRO_RUN_DIR="$(cd "${INVITRO_RUN_DIR}" && pwd)"; else INVITRO_RUN_DIR=""; fi

if [[ -z "${PARAMETER_TABLE}" ]]; then PARAMETER_TABLE="${PROJECT_ROOT}/oxygen/data/O2_supply_demand/parameter_table_invitro_buffering.csv"; fi
if [[ -z "${FIT_OBJECTS_DIR}" ]]; then FIT_OBJECTS_DIR="${PROJECT_ROOT}/oxygen/ploidyOxygen/data/fit_objects"; fi
if [[ -z "${FLOW_DENSITY_PATH}" ]]; then FLOW_DENSITY_PATH="${PROJECT_ROOT}/oxygen/data/g0g1_ploidy_density_grid.csv"; fi
PARAMETER_TABLE="$(cd "$(dirname "${PARAMETER_TABLE}")" && pwd)/$(basename "${PARAMETER_TABLE}")"
FIT_OBJECTS_DIR="$(cd "${FIT_OBJECTS_DIR}" && pwd)"
FLOW_DENSITY_PATH="$(cd "$(dirname "${FLOW_DENSITY_PATH}")" && pwd)/$(basename "${FLOW_DENSITY_PATH}")"

MULTI_WARMUP_ROOT="${OUT_ROOT}/${MULTI_WARMUP_PREFIX}"
mkdir -p "${MULTI_WARMUP_ROOT}"
PROGRESS_LOG="${MULTI_WARMUP_ROOT}/multi_warmup_progress.log"
JOBS_TSV="${MULTI_WARMUP_ROOT}/multi_warmup_jobs.tsv"
: > "${PROGRESS_LOG}"

RUN_O2_FIT_SCRIPT="${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/runner/run_o2_fit.sh"
LEGACY_SEED_PLAN_SCRIPT="${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/analysis/multi_warmup/build_multi_warmup_seed_plan.R"
LANDSCAPE_SEED_PLAN_SCRIPT="${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/analysis/multi_warmup/build_multi_warmup_pairs_from_landscape_subclusters.R"
if [[ "${MULTI_WARMUP_PAIR_METHOD}" == "landscape_subcluster" ]]; then
  SEED_PLAN_SCRIPT="${LANDSCAPE_SEED_PLAN_SCRIPT}"
else
  SEED_PLAN_SCRIPT="${LEGACY_SEED_PLAN_SCRIPT}"
fi
COLLECT_SCRIPT="${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/analysis/multi_warmup/collect_multi_warmup_results.R"
REPORT_SCRIPT="${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/analysis/multi_warmup/multi_warmup_results_report.R"

for path in "${RUN_O2_FIT_SCRIPT}" "${SEED_PLAN_SCRIPT}" "${COLLECT_SCRIPT}" "${REPORT_SCRIPT}" "${CONFIG_PATH}" "${PARAMETER_TABLE}"; do
  if [[ ! -f "${path}" ]]; then
    echo "Missing required file: ${path}" >&2
    exit 1
  fi
done

log_msg "stage=prepare pair_method=${MULTI_WARMUP_PAIR_METHOD} source_invivo=${INVIVO_RUN_DIR} source_invitro=${INVITRO_RUN_DIR}"
if [[ "${MULTI_WARMUP_PAIR_METHOD}" == "landscape_subcluster" ]]; then
  seed_plan_cmd=(
    Rscript "${SEED_PLAN_SCRIPT}"
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
  if [[ -n "${MULTI_WARMUP_LANDSCAPE_MAX_SEEDS}" ]]; then seed_plan_cmd+=("--max_seeds=${MULTI_WARMUP_LANDSCAPE_MAX_SEEDS}"); fi
  if [[ -n "${MULTI_WARMUP_REFERENCE_SUBCLUSTER_DIR}" ]]; then seed_plan_cmd+=("--reference_subcluster_dir=${MULTI_WARMUP_REFERENCE_SUBCLUSTER_DIR}"); fi
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
run_or_print "Generate multi-warmup seed plan" "${seed_plan_cmd[@]}"

MANIFEST="${MULTI_WARMUP_ROOT}/multi_warmup_manifest.tsv"
if ! truthy "${DRY_RUN}" && [[ ! -f "${MANIFEST}" ]]; then
  echo "Missing generated manifest: ${MANIFEST}" >&2
  exit 1
fi

if truthy "${DRY_RUN}"; then
  log_msg "stage=dry_run complete_after_seed_plan_command"
  exit 0
fi

total_pairs=$(( $(wc -l < "${MANIFEST}") - 1 ))
PLAN_MODE_FILE="${MULTI_WARMUP_ROOT}/multi_warmup_seed_plan_mode.tsv"
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
printf "warmup_label\tpair_index\ttotal_pairs\tjoint_run_dir\trun_status\textra_results_status\tnext_step\n" > "${JOBS_TSV}"
log_msg "stage=run_pairs total_pairs=${total_pairs}"

pair_index=0
tail -n +2 "${MANIFEST}" | while IFS=$'\t' read -r warmup_label phase invivo_family invivo_rank invivo_seed invivo_seed_dir invitro_family invitro_rank invitro_seed invitro_seed_dir selection_reason joint_run_prefix joint_table; do
  pair_index=$((pair_index + 1))
  pair_log="${MULTI_WARMUP_ROOT}/${warmup_label}.run.log"
  config_dir="${MULTI_WARMUP_ROOT}/configs"
  mkdir -p "${config_dir}" "$(dirname "${joint_table}")"
  pair_config="${config_dir}/${warmup_label}__$(basename "${CONFIG_PATH}")"
  cp "${CONFIG_PATH}" "${pair_config}"
  joint_run_dir="${MULTI_WARMUP_ROOT}/${joint_run_prefix}"
  log_msg "stage=pair_start pair=${pair_index}/${total_pairs} warmup_label=${warmup_label} run_dir=${joint_run_dir}"
  cmd=(
    bash "${RUN_O2_FIT_SCRIPT}"
    --fitting_mode=joint
    --joint_fitting_mode=DIRECT
    "--project_root=${PROJECT_ROOT}"
    "--config_path=${pair_config}"
    "--out_root=${MULTI_WARMUP_ROOT}"
    "--joint_run_prefix=${joint_run_prefix}"
    "--joint_seeds_csv=${JOINT_SEEDS_CSV}"
    "--joint_n_cores=${JOINT_N_CORES}"
    "--r_module=${R_MODULE}"
    "--run_extra_results=${RUN_EXTRA_RESULTS}"
    "--force_extra_results=${FORCE_EXTRA_RESULTS}"
    --append_run_prefix_timestamp=FALSE
    "--parameter_table=${PARAMETER_TABLE}"
    "--fit_objects_dir=${FIT_OBJECTS_DIR}"
    "--flow_density_path=${FLOW_DENSITY_PATH}"
    "--itermax=${ITERMAX}"
    "--de_reltol=${DE_RELTOL}"
    "--de_steptol=${DE_STEPTOL}"
    "--NP=${NP}"
    "--auto_viz=${AUTO_VIZ}"
    "--invivo_best_seed_dir=${invivo_seed_dir}"
    "--invitro_best_seed_dir=${invitro_seed_dir}"
    "--joint_warmup_seed_label=${warmup_label}"
    "--joint_soft_coupling_parameters_table=${joint_table}"
    "--joint_soft_coupling_delta_params=${JOINT_SOFT_COUPLING_DELTA_PARAMS}"
  )
  if [[ -n "${JOINT_WARMUP_SIGMAN}" ]]; then cmd+=("--joint_warmup_sigmaN=${JOINT_WARMUP_SIGMAN}"); fi
  if [[ -n "${JOINT_SOFT_COUPLING_SIGMA_DEFAULT}" ]]; then cmd+=("--joint_soft_coupling_sigma_default=${JOINT_SOFT_COUPLING_SIGMA_DEFAULT}"); fi
  if [[ -n "${JOINT_SOFT_COUPLING_WELSCH_C}" ]]; then cmd+=("--joint_soft_coupling_welsch_c=${JOINT_SOFT_COUPLING_WELSCH_C}"); fi
  print_command "Run pair ${pair_index}/${total_pairs}" "${cmd[@]}"
  if "${cmd[@]}" > "${pair_log}" 2>&1; then
    extra_status="done"
    run_status="done"
    log_msg "stage=pair_done pair=${pair_index}/${total_pairs} warmup_label=${warmup_label} log=${pair_log}"
  else
    extra_status="unknown"
    run_status="failed"
    log_msg "stage=pair_failed pair=${pair_index}/${total_pairs} warmup_label=${warmup_label} log=${pair_log}"
    tail -40 "${pair_log}" >&2 || true
    exit 1
  fi
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "${warmup_label}" "${pair_index}" "${total_pairs}" "${joint_run_dir}" "${run_status}" "${extra_status}" "next_pair_or_collect" >> "${JOBS_TSV}"
done

log_msg "stage=collect next=collect_multi_warmup_results"
run_or_print "Collect multi-warmup results" Rscript "${COLLECT_SCRIPT}" "--multi_warmup_root=${MULTI_WARMUP_ROOT}" "--manifest=${MANIFEST}" "--out_dir=${MULTI_WARMUP_ROOT}"
log_msg "stage=report next=multi_warmup_results_report"
run_or_print "Render multi-warmup report" Rscript "${REPORT_SCRIPT}" "--multi_warmup_root=${MULTI_WARMUP_ROOT}" "--out=${MULTI_WARMUP_ROOT}/multi-warm-up_results.html"
log_msg "stage=complete report=${MULTI_WARMUP_ROOT}/multi-warm-up_results.html"
