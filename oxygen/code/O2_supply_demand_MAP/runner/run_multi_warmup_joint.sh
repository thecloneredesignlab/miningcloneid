#!/usr/bin/env bash
# Local joint-fitting pipeline: pooled t-SNE and bilateral primary-cluster representatives.

set -euo pipefail

O2SD_SHELL_UTILS="$(cd "$(dirname "${BASH_SOURCE[0]}")/../util" && pwd)/o2_supply_demand_map_shell_utils.sh"
# shellcheck source=../util/o2_supply_demand_map_shell_utils.sh
source "${O2SD_SHELL_UTILS}"

usage() {
  cat <<'EOF'
Usage:
  bash run_multi_warmup_joint.sh --invivo_run_dir=DIR --invitro_run_dir=DIR [options]

This is the only joint-fitting workflow. It builds one pooled parameter-space
t-SNE, clusters the in-vivo and in-vitro best points separately, selects the
objective-minimum seed from every primary cluster on both sides, and creates
their Cartesian product. No second-level clusters or curve filters run.

Options:
  --project_root=DIR
  --config_path=FILE
  --out_root=DIR
  --joint_run_prefix=NAME
  --joint_seeds_csv=1,2,3
  --joint_n_cores=N
  --joint_tsne_seed=123
  --joint_cluster_seed=123
  --joint_landscape_max_seeds=N
  --joint_landscape_n_threads=N
  --dry_run=TRUE|FALSE
  --run_extra_results=TRUE|FALSE
  --force_extra_results=TRUE|FALSE
EOF
}

parse_args() {
  for arg in "$@"; do
    case "${arg}" in
      --help|-h) usage; exit 0 ;;
      --project_root=*) PROJECT_ROOT="${arg#*=}" ;;
      --config_path=*|--config=*) CONFIG_PATH="${arg#*=}" ;;
      --out_root=*) OUT_ROOT="${arg#*=}" ;;
      --multi_warmup_prefix=*|--joint_run_prefix=*) JOINT_RUN_PREFIX="${arg#*=}" ;;
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
      --joint_tsne_seed=*|--multi_warmup_tsne_seed=*|--landscape_tsne_seed=*) JOINT_TSNE_SEED="${arg#*=}" ;;
      --joint_cluster_seed=*|--multi_warmup_cluster_seed=*|--landscape_cluster_seed=*) JOINT_CLUSTER_SEED="${arg#*=}" ;;
      --joint_landscape_max_seeds=*|--multi_warmup_landscape_max_seeds=*|--landscape_max_seeds=*) JOINT_LANDSCAPE_MAX_SEEDS="${arg#*=}" ;;
      --joint_landscape_n_threads=*|--landscape_n_threads=*) JOINT_LANDSCAPE_N_THREADS="${arg#*=}" ;;
      --joint_fitting_mode=*|--multi_warmup_pair_method=*|--pair_method=*|--multi_warmup_reductions=*|--landscape_reductions=*|--multi_warmup_subcluster_seed=*|--landscape_subcluster_seed=*|--multi_warmup_pairing_policy=*|--pairing_policy=*|--multi_warmup_invivo_curve_filter=*|--invivo_curve_filter=*|--multi_warmup_invivo_curve_class=*|--invivo_curve_class=*)
        echo "${arg%%=*} has been removed; joint fitting now has one fixed bilateral primary-cluster workflow." >&2
        exit 2
        ;;
      *) echo "Unknown argument: ${arg}" >&2; usage >&2; exit 2 ;;
    esac
  done
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKFLOW_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
DEFAULT_PROJECT_ROOT="$(cd "${WORKFLOW_ROOT}/../../.." && pwd)"
DEFAULT_ITERMAX="1000"

PROJECT_ROOT="${PROJECT_ROOT:-${DEFAULT_PROJECT_ROOT}}"
CONFIG_PATH="${CONFIG_PATH:-}"
OUT_ROOT="${OUT_ROOT:-}"
JOINT_RUN_PREFIX="${JOINT_RUN_PREFIX:-fit_joint_primary_clusters}"
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
ITERMAX="${ITERMAX:-${DEFAULT_ITERMAX}}"
DE_RELTOL="${DE_RELTOL:-1e-4}"
DE_STEPTOL="${DE_STEPTOL:-25}"
NP="${NP:-80}"
AUTO_VIZ="${AUTO_VIZ:-TRUE}"
JOINT_WARMUP_SIGMAN="${JOINT_WARMUP_SIGMAN:-}"
JOINT_SOFT_COUPLING_SIGMA_DEFAULT="${JOINT_SOFT_COUPLING_SIGMA_DEFAULT:-}"
JOINT_SOFT_COUPLING_WELSCH_C="${JOINT_SOFT_COUPLING_WELSCH_C:-}"
JOINT_SOFT_COUPLING_DELTA_PARAMS="${JOINT_SOFT_COUPLING_DELTA_PARAMS:-default}"
JOINT_TSNE_SEED="${JOINT_TSNE_SEED:-123}"
JOINT_CLUSTER_SEED="${JOINT_CLUSTER_SEED:-123}"
JOINT_LANDSCAPE_MAX_SEEDS="${JOINT_LANDSCAPE_MAX_SEEDS:-}"
JOINT_LANDSCAPE_N_THREADS="${JOINT_LANDSCAPE_N_THREADS:-8}"

parse_args "$@"

if is_null_value "${INVIVO_RUN_DIR}" || is_null_value "${INVITRO_RUN_DIR}"; then
  echo "--invivo_run_dir and --invitro_run_dir are both required for joint fitting." >&2
  exit 2
fi
require_positive_int JOINT_N_CORES "${JOINT_N_CORES}"
require_positive_int ITERMAX "${ITERMAX}"
require_positive_int DE_STEPTOL "${DE_STEPTOL}"
require_positive_int NP "${NP}"
require_positive_int JOINT_TSNE_SEED "${JOINT_TSNE_SEED}"
require_positive_int JOINT_CLUSTER_SEED "${JOINT_CLUSTER_SEED}"
require_positive_int JOINT_LANDSCAPE_N_THREADS "${JOINT_LANDSCAPE_N_THREADS}"

PROJECT_ROOT="$(cd "${PROJECT_ROOT}" && pwd)"
CONFIG_PATH="${CONFIG_PATH:-${PROJECT_ROOT}/oxygen/config/O2_supply_demand.yaml}"
OUT_ROOT="${OUT_ROOT:-${PROJECT_ROOT}/oxygen/results}"
PARAMETER_TABLE="${PARAMETER_TABLE:-${PROJECT_ROOT}/oxygen/data/O2_supply_demand/parameter_table_invitro_buffering.csv}"
FIT_OBJECTS_DIR="${FIT_OBJECTS_DIR:-${PROJECT_ROOT}/oxygen/ploidyOxygen/data/fit_objects}"
FLOW_DENSITY_PATH="${FLOW_DENSITY_PATH:-${PROJECT_ROOT}/oxygen/data/g0g1_ploidy_density_grid.csv}"
INVIVO_RUN_DIR="$(cd "${INVIVO_RUN_DIR}" && pwd)"
INVITRO_RUN_DIR="$(cd "${INVITRO_RUN_DIR}" && pwd)"
mkdir -p "${OUT_ROOT}"
OUT_ROOT="$(cd "${OUT_ROOT}" && pwd)"
JOINT_ROOT="${OUT_ROOT}/${JOINT_RUN_PREFIX}"
mkdir -p "${JOINT_ROOT}"

PAIR_SCRIPT="${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/analysis/multi_warmup/build_joint_primary_cluster_pairs.R"
JOINT_RUNNER="${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/runner/run_fit_joint_model_O2_supply_demand_MAP.sh"
WARM_START_SCRIPT="${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/runner/warm_start/make_joint_soft_coupling_parameters_table.R"
EXTRA_RESULTS_SCRIPT="${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/analysis/fit_results/extra_results.R"
for path in "${CONFIG_PATH}" "${PARAMETER_TABLE}" "${PAIR_SCRIPT}" "${JOINT_RUNNER}" "${WARM_START_SCRIPT}" "${EXTRA_RESULTS_SCRIPT}"; do
  [[ -f "${path}" ]] || { echo "Missing required file: ${path}" >&2; exit 1; }
done
[[ -d "${FIT_OBJECTS_DIR}" ]] || { echo "Missing fit_objects_dir: ${FIT_OBJECTS_DIR}" >&2; exit 1; }

seed_plan_cmd=(Rscript "${PAIR_SCRIPT}" "--project_root=${PROJECT_ROOT}" "--invivo_run_dir=${INVIVO_RUN_DIR}" "--invitro_run_dir=${INVITRO_RUN_DIR}" "--out_dir=${JOINT_ROOT}" "--tsne_seed=${JOINT_TSNE_SEED}" "--cluster_seed=${JOINT_CLUSTER_SEED}" "--n_threads=${JOINT_LANDSCAPE_N_THREADS}" "--render_figures=FALSE")
if [[ -n "${JOINT_LANDSCAPE_MAX_SEEDS}" ]]; then seed_plan_cmd+=("--max_seeds=${JOINT_LANDSCAPE_MAX_SEEDS}"); fi
run_or_print "Build fixed joint primary-cluster manifest" "${seed_plan_cmd[@]}"
if truthy "${DRY_RUN}"; then exit 0; fi

MANIFEST="${JOINT_ROOT}/multi_warmup_manifest.tsv"
[[ -f "${MANIFEST}" ]] || { echo "Missing generated manifest: ${MANIFEST}" >&2; exit 1; }

tail -n +2 "${MANIFEST}" | while IFS=$'\t' read -r warmup_label phase invivo_family invivo_rank invivo_seed invivo_seed_dir invitro_family invitro_rank invitro_seed invitro_seed_dir selection_reason joint_run_prefix joint_table; do
  mkdir -p "$(dirname "${joint_table}")"
  table_cmd=(Rscript "${WARM_START_SCRIPT}" "--invivo-seed-dir=${invivo_seed_dir}" "--invitro-seed-dir=${invitro_seed_dir}" "--seed-label=${warmup_label}" "--out=${joint_table}" "--delta-params=${JOINT_SOFT_COUPLING_DELTA_PARAMS}")
  run_or_print "Generate ${warmup_label} soft-coupling table" "${table_cmd[@]}"

  warm_args=("--joint_warmup_enable=TRUE" "--joint_warmup_seed_label=${warmup_label}" "--joint_warmup_invivo_seed_dir=${invivo_seed_dir}" "--joint_warmup_invitro_seed_dir=${invitro_seed_dir}" "--joint_soft_coupling_parameters_table=${joint_table}")
  [[ -n "${JOINT_WARMUP_SIGMAN}" ]] && warm_args+=("--joint_warmup_sigmaN=${JOINT_WARMUP_SIGMAN}")
  [[ -n "${JOINT_SOFT_COUPLING_SIGMA_DEFAULT}" ]] && warm_args+=("--joint_soft_coupling_sigma_default=${JOINT_SOFT_COUPLING_SIGMA_DEFAULT}")
  [[ -n "${JOINT_SOFT_COUPLING_WELSCH_C}" ]] && warm_args+=("--joint_soft_coupling_welsch_c=${JOINT_SOFT_COUPLING_WELSCH_C}")
  flow_args=()
  [[ -f "${FLOW_DENSITY_PATH}" ]] && flow_args+=("--flow_density_path=${FLOW_DENSITY_PATH}")
  run_or_print "Run ${warmup_label} joint fit" bash "${JOINT_RUNNER}" --mode=run "--config=${CONFIG_PATH}" "--out_root=${JOINT_ROOT}" "--run_prefix=${joint_run_prefix}" --append_run_prefix_timestamp=FALSE "--seeds_csv=${JOINT_SEEDS_CSV}" "--n_cores=${JOINT_N_CORES}" "--auto_viz=${AUTO_VIZ}" "--itermax=${ITERMAX}" "--de_reltol=${DE_RELTOL}" "--de_steptol=${DE_STEPTOL}" "--NP=${NP}" "--invitro_parameter_table=${PARAMETER_TABLE}" "--fit_objects_dir=${FIT_OBJECTS_DIR}" "${warm_args[@]}" "${flow_args[@]}"

  pair_dir="${JOINT_ROOT}/${joint_run_prefix}"
  summary="${pair_dir}/extra_results/seed_summary.tsv"
  if truthy "${RUN_EXTRA_RESULTS}" && { truthy "${FORCE_EXTRA_RESULTS}" || [[ ! -f "${summary}" ]]; }; then
    run_or_print "Run ${warmup_label} extra_results" Rscript "${EXTRA_RESULTS_SCRIPT}" "--run_dir=${pair_dir}"
  fi
done
