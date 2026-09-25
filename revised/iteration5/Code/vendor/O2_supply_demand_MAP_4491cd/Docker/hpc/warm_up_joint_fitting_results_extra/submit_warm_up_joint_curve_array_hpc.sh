#!/usr/bin/env bash
# Submit joint warm-up curve classification as a Slurm array workflow.

set -euo pipefail

O2SD_DOCKER_HPC_ROOT="${O2SD_DOCKER_HPC_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
# shellcheck source=../util/o2_supply_demand_map_apptainer_runtime.sh
source "${O2SD_DOCKER_HPC_ROOT}/util/o2_supply_demand_map_apptainer_runtime.sh"

O2SD_SHELL_UTILS="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../util" && pwd)/o2_supply_demand_map_shell_utils.sh"
# shellcheck source=../../../util/o2_supply_demand_map_shell_utils.sh
source "${O2SD_SHELL_UTILS}"
O2SD_R_MODULE_OPTIONAL=TRUE

if [[ "${O2_WARMUP_JOINT_CURVE_ARRAY_LOGIN_SHELL:-0}" != "1" ]]; then
  export O2_WARMUP_JOINT_CURVE_ARRAY_LOGIN_SHELL=1
  exec bash -lc "$(shell_join bash "$0" "$@")"
fi

usage() {
  cat <<'EOF'
Usage:
  bash submit_warm_up_joint_curve_array_hpc.sh [options]

This submits the joint warm-up curve workflow as:
  1. dense-grid build_tasks
  2. dense-grid Slurm array run_tasks
  3. dense-grid merge
  4. dependent curve_regression + summary postprocess job

Important paths:
  --project_root=DIR
  --input_root=DIR
  --output_root=DIR
  --analysis_script=FILE
  --dense_submitter=FILE
  --array_backend=FILE
  --log_root=DIR

Task splitting:
  --tasks_per_array_task=N        Default 1000 seed x O2 calculations per Slurm array task.
  --array_max_concurrent=N        Default 500.

Resources:
  --qos=NAME                      Default xxlarge.
  --r_module=MODULE               Compatibility option; ignored by Docker/hpc.
  --array_cpus=N                  Default 1.
  --array_mem=SIZE                Default 2G.
  --array_time=TIME               Default 12:00:00.
  --merge_cpus=N                  Default 8.
  --merge_mem=SIZE                Default 32G.
  --merge_time=TIME               Default 4:00:00.
  --postprocess_cpus=N            Default 4.
  --postprocess_mem=SIZE          Default 32G.
  --postprocess_time=TIME         Default 4:00:00.

Analysis:
  --o2_grid=CSV
  --reporting_o2=CSV
  --overwrite=TRUE|FALSE          Default TRUE.
  --generate_figures=TRUE|FALSE   Default TRUE.
  --run_validation=TRUE|FALSE     Default TRUE.
  --dependency=SPEC               Extra dependency added to the dense-grid array job.
  --dry_run=TRUE|FALSE            Build scripts/task lists but do not submit jobs.
  --skip_postprocess=TRUE|FALSE   Only submit dense-grid array + merge.
EOF
}

resolve_path() {
  local base="$1"
  local path="${2:-}"
  if [[ -z "${path}" ]]; then
    return 0
  fi
  case "${path}" in
    "~") printf "%s" "${HOME}" ;;
    "~/"*) printf "%s/%s" "${HOME}" "${path#"~/"}" ;;
    /*) printf "%s" "${path}" ;;
    *) printf "%s/%s" "${base}" "${path}" ;;
  esac
}

write_postprocess_sbatch() {
  local script_path="$1"
  local stdout_path="$2"
  local stderr_path="$3"
  {
    printf '#!/usr/bin/env bash\n'
    printf '#SBATCH --job-name=warmup_curve_post\n'
    printf '%s\n' '#SBATCH --ntasks=1'
    printf '#SBATCH --cpus-per-task=%s\n' "${POSTPROCESS_CPUS}"
    printf '#SBATCH --mem=%s\n' "${POSTPROCESS_MEM}"
    printf '#SBATCH --time=%s\n' "${POSTPROCESS_TIME}"
    if [[ -n "${QOS}" ]]; then printf '#SBATCH --qos=%s\n' "${QOS}"; fi
    if [[ -n "${PARTITION}" ]]; then printf '#SBATCH --partition=%s\n' "${PARTITION}"; fi
    if [[ -n "${ACCOUNT}" ]]; then printf '#SBATCH --account=%s\n' "${ACCOUNT}"; fi
    printf '#SBATCH --output=%s\n' "${stdout_path}"
    printf '#SBATCH --error=%s\n' "${stderr_path}"
    printf '\n'
    printf 'set -euo pipefail\n'
    printf 'source %q\n' "${O2SD_DOCKER_HPC_ROOT}/util/o2_supply_demand_map_apptainer_runtime.sh"
    printf 'RSCRIPT_BIN="$(command -v Rscript || true)"\n'
    printf 'if [[ -z "${RSCRIPT_BIN}" ]]; then echo "Container-backed Rscript not found" >&2; exit 1; fi\n'
    printf 'cd %q\n' "${PROJECT_ROOT}"
    printf 'echo "Host: $(hostname)"\n'
    printf 'echo "Working directory: $(pwd)"\n'
    printf 'echo "Rscript: ${RSCRIPT_BIN}"\n'
    printf '"${RSCRIPT_BIN}" --version\n'
    printf 'cmd=( "${RSCRIPT_BIN}" %q --stage=curve_regression,summary --input_root=%q --output_root=%q --overwrite=%q --generate_figures=%q )\n' \
      "${ANALYSIS_SCRIPT}" "${INPUT_ROOT}" "${OUTPUT_ROOT}" "${OVERWRITE}" "${GENERATE_FIGURES}"
    printf 'printf "Command:"\n'
    printf 'printf " %%q" "${cmd[@]}"\n'
    printf 'printf "\\n"\n'
    printf '"${cmd[@]}"\n'
  } > "${script_path}"
  chmod +x "${script_path}"
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKFLOW_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
DEFAULT_PROJECT_ROOT="$(cd "${WORKFLOW_ROOT}/../../.." && pwd)"
DEFAULT_RESULT_NAME="fit_joint_multi_warmup_tsne_sigmaN0p1216_objmin_500seed_20260714_021540"

PROJECT_ROOT="${PROJECT_ROOT:-${DEFAULT_PROJECT_ROOT}}"
INPUT_ROOT=""
OUTPUT_ROOT=""
ANALYSIS_SCRIPT=""
DENSE_SUBMITTER=""
ARRAY_BACKEND=""
LOG_ROOT=""
TASKS_PER_ARRAY_TASK="${TASKS_PER_ARRAY_TASK:-1000}"
ARRAY_MAX_CONCURRENT="${ARRAY_MAX_CONCURRENT:-500}"
QOS="${QOS:-xxlarge}"
PARTITION="${PARTITION:-}"
ACCOUNT="${ACCOUNT:-}"
R_MODULE="${R_MODULE:-R/4.4}"
ARRAY_CPUS="${ARRAY_CPUS:-1}"
ARRAY_MEM="${ARRAY_MEM:-2G}"
ARRAY_TIME="${ARRAY_TIME:-12:00:00}"
MERGE_CPUS="${MERGE_CPUS:-8}"
MERGE_MEM="${MERGE_MEM:-32G}"
MERGE_TIME="${MERGE_TIME:-4:00:00}"
POSTPROCESS_CPUS="${POSTPROCESS_CPUS:-4}"
POSTPROCESS_MEM="${POSTPROCESS_MEM:-32G}"
POSTPROCESS_TIME="${POSTPROCESS_TIME:-4:00:00}"
O2_GRID="${O2_GRID:-}"
REPORTING_O2="${REPORTING_O2:-}"
OVERWRITE="${OVERWRITE:-TRUE}"
GENERATE_FIGURES="${GENERATE_FIGURES:-TRUE}"
RUN_VALIDATION="${RUN_VALIDATION:-TRUE}"
EXTRA_DEPENDENCY="${EXTRA_DEPENDENCY:-}"
DRY_RUN="${DRY_RUN:-FALSE}"
SKIP_POSTPROCESS="${SKIP_POSTPROCESS:-FALSE}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --help|-h) usage; exit 0 ;;
    --project_root=*) PROJECT_ROOT="${1#*=}" ;;
    --input_root=*) INPUT_ROOT="${1#*=}" ;;
    --output_root=*) OUTPUT_ROOT="${1#*=}" ;;
    --analysis_script=*|--script=*) ANALYSIS_SCRIPT="${1#*=}" ;;
    --dense_submitter=*) DENSE_SUBMITTER="${1#*=}" ;;
    --array_backend=*) ARRAY_BACKEND="${1#*=}" ;;
    --log_root=*) LOG_ROOT="${1#*=}" ;;
    --tasks_per_array_task=*) TASKS_PER_ARRAY_TASK="${1#*=}" ;;
    --array_max_concurrent=*) ARRAY_MAX_CONCURRENT="${1#*=}" ;;
    --qos=*) QOS="${1#*=}" ;;
    --partition=*) PARTITION="${1#*=}" ;;
    --account=*) ACCOUNT="${1#*=}" ;;
    --r_module=*) R_MODULE="${1#*=}" ;;
    --array_cpus=*|--monotonicity_array_cpus=*) ARRAY_CPUS="${1#*=}" ;;
    --array_mem=*|--monotonicity_array_mem=*) ARRAY_MEM="${1#*=}" ;;
    --array_time=*|--monotonicity_array_time=*) ARRAY_TIME="${1#*=}" ;;
    --merge_cpus=*|--monotonicity_merge_cpus=*) MERGE_CPUS="${1#*=}" ;;
    --merge_mem=*|--monotonicity_merge_mem=*) MERGE_MEM="${1#*=}" ;;
    --merge_time=*|--monotonicity_merge_time=*) MERGE_TIME="${1#*=}" ;;
    --postprocess_cpus=*) POSTPROCESS_CPUS="${1#*=}" ;;
    --postprocess_mem=*) POSTPROCESS_MEM="${1#*=}" ;;
    --postprocess_time=*) POSTPROCESS_TIME="${1#*=}" ;;
    --o2_grid=*) O2_GRID="${1#*=}" ;;
    --reporting_o2=*) REPORTING_O2="${1#*=}" ;;
    --overwrite=*|--force=*) OVERWRITE="${1#*=}" ;;
    --generate_figures=*) GENERATE_FIGURES="${1#*=}" ;;
    --run_validation=*) RUN_VALIDATION="${1#*=}" ;;
    --dependency=*) EXTRA_DEPENDENCY="${1#*=}" ;;
    --dry_run=*) DRY_RUN="${1#*=}" ;;
    --skip_postprocess=*) SKIP_POSTPROCESS="${1#*=}" ;;
    *)
      echo "Unknown argument: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
  shift
done

PROJECT_ROOT="$(cd "${PROJECT_ROOT}" && pwd)"
if [[ -z "${INPUT_ROOT}" ]]; then
  INPUT_ROOT="${PROJECT_ROOT}/oxygen/results/${DEFAULT_RESULT_NAME}"
else
  INPUT_ROOT="$(resolve_path "${PROJECT_ROOT}" "${INPUT_ROOT}")"
fi
if [[ -z "${OUTPUT_ROOT}" ]]; then
  OUTPUT_ROOT="${PROJECT_ROOT}/oxygen/results/analysis/warm_up_joint_fitting_results_extra/$(basename "${INPUT_ROOT}")"
else
  OUTPUT_ROOT="$(resolve_path "${PROJECT_ROOT}" "${OUTPUT_ROOT}")"
fi
if [[ -z "${ANALYSIS_SCRIPT}" ]]; then
  ANALYSIS_SCRIPT="${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/analysis/warm_up_joint_fitting_results_extra/warm_up_joint_fitting_results_extra.R"
else
  ANALYSIS_SCRIPT="$(resolve_path "${PROJECT_ROOT}" "${ANALYSIS_SCRIPT}")"
fi
if [[ -z "${DENSE_SUBMITTER}" ]]; then
  DENSE_SUBMITTER="${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/Docker/hpc/dense_grid_monotonicity_classification/submit_dense_grid_monotonicity_classification.sh"
else
  DENSE_SUBMITTER="$(resolve_path "${PROJECT_ROOT}" "${DENSE_SUBMITTER}")"
fi
if [[ -z "${ARRAY_BACKEND}" ]]; then
  ARRAY_BACKEND="${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/runner/dense_grid_monotonicity/run_dense_grid_monotonicity.R"
else
  ARRAY_BACKEND="$(resolve_path "${PROJECT_ROOT}" "${ARRAY_BACKEND}")"
fi
if [[ -z "${LOG_ROOT}" ]]; then
  LOG_ROOT="${OUTPUT_ROOT}_hpc_logs/curve_array"
else
  LOG_ROOT="$(resolve_path "${PROJECT_ROOT}" "${LOG_ROOT}")"
fi

SEED_MANIFEST="${OUTPUT_ROOT}/tables/joint_best_curve_seed_manifest.tsv"
DENSE_OUT_DIR="${OUTPUT_ROOT}/curve_classification/dense-grid_monotonicity_classification"

if [[ ! -f "${SEED_MANIFEST}" ]]; then
  echo "Missing seed manifest. Run warm_up_joint_fitting_results_extra.R --stage=prepare first: ${SEED_MANIFEST}" >&2
  exit 1
fi
if [[ ! -f "${ANALYSIS_SCRIPT}" ]]; then
  echo "Missing analysis script: ${ANALYSIS_SCRIPT}" >&2
  exit 1
fi
if [[ ! -f "${DENSE_SUBMITTER}" ]]; then
  echo "Missing dense-grid submitter: ${DENSE_SUBMITTER}" >&2
  exit 1
fi
if [[ ! -f "${ARRAY_BACKEND}" ]]; then
  echo "Missing array backend: ${ARRAY_BACKEND}" >&2
  exit 1
fi
if ! truthy "${DRY_RUN}" && ! command -v sbatch >/dev/null 2>&1; then
  echo "sbatch not found. Run this submitter on the HPC login node, or rerun with --dry_run=TRUE." >&2
  exit 1
fi

load_r_module
if ! command -v Rscript >/dev/null 2>&1; then
  echo "Container-backed Rscript not found; SIF: ${O2SD_CONTAINER_IMAGE}." >&2
  exit 1
fi

STAMP="$(date '+%Y%m%d_%H%M%S')"
RUN_LOG_DIR="${LOG_ROOT}/warm_up_joint_curve_array_${STAMP}"
mkdir -p "${RUN_LOG_DIR}"
SUBMIT_CAPTURE="${RUN_LOG_DIR}/dense_submit_output.log"
SUMMARY_MANIFEST="${RUN_LOG_DIR}/warm_up_joint_curve_array_manifest.tsv"

dense_cmd=(
  bash "${DENSE_SUBMITTER}"
	  "--project_root=${PROJECT_ROOT}"
	  "--run_parts=monotonicity"
	  "--run_dir=${INPUT_ROOT}"
	  "--seed_manifest=${SEED_MANIFEST}"
	  "--result_root=${OUTPUT_ROOT}/curve_classification"
  "--monotonicity_out_dir=${DENSE_OUT_DIR}"
  "--array_backend=${ARRAY_BACKEND}"
  "--log_root=${RUN_LOG_DIR}/dense_grid_logs"
  "--tasks_per_array_task=${TASKS_PER_ARRAY_TASK}"
  "--array_max_concurrent=${ARRAY_MAX_CONCURRENT}"
  "--qos=${QOS}"
  "--r_module=${R_MODULE}"
  "--monotonicity_array_cpus=${ARRAY_CPUS}"
  "--monotonicity_array_mem=${ARRAY_MEM}"
  "--monotonicity_array_time=${ARRAY_TIME}"
  "--monotonicity_merge_cpus=${MERGE_CPUS}"
  "--monotonicity_merge_mem=${MERGE_MEM}"
  "--monotonicity_merge_time=${MERGE_TIME}"
  "--overwrite=${OVERWRITE}"
  "--generate_figures=${GENERATE_FIGURES}"
  "--run_validation=${RUN_VALIDATION}"
  "--dry_run=${DRY_RUN}"
)
if [[ -n "${PARTITION}" ]]; then dense_cmd+=("--partition=${PARTITION}"); fi
if [[ -n "${ACCOUNT}" ]]; then dense_cmd+=("--account=${ACCOUNT}"); fi
if [[ -n "${O2_GRID}" ]]; then dense_cmd+=("--o2_grid=${O2_GRID}"); fi
if [[ -n "${REPORTING_O2}" ]]; then dense_cmd+=("--reporting_o2=${REPORTING_O2}"); fi
if [[ -n "${EXTRA_DEPENDENCY}" ]]; then dense_cmd+=("--dependency=${EXTRA_DEPENDENCY}"); fi

echo "Submitting dense-grid array workflow:"
echo "  $(shell_join "${dense_cmd[@]}")"
"${dense_cmd[@]}" 2>&1 | tee "${SUBMIT_CAPTURE}"

dense_manifest="$(awk -F'Submission manifest: ' '/Submission manifest:/ {print $2}' "${SUBMIT_CAPTURE}" | tail -1)"
if [[ -z "${dense_manifest}" || ! -f "${dense_manifest}" ]]; then
  echo "Could not find dense-grid submission manifest in: ${SUBMIT_CAPTURE}" >&2
  exit 1
fi
merge_job_id="$(awk -F '\t' '$1 == "monotonicity_merge" {print $2}' "${dense_manifest}" | tail -1)"
if [[ -z "${merge_job_id}" ]]; then
  echo "Could not find monotonicity_merge job id in: ${dense_manifest}" >&2
  exit 1
fi

postprocess_job_id=""
postprocess_sbatch=""
if ! truthy "${SKIP_POSTPROCESS}"; then
  postprocess_sbatch="${RUN_LOG_DIR}/curve_regression_summary.sbatch"
  postprocess_out="${RUN_LOG_DIR}/curve_regression_summary_%j.out"
  postprocess_err="${RUN_LOG_DIR}/curve_regression_summary_%j.err"
  write_postprocess_sbatch "${postprocess_sbatch}" "${postprocess_out}" "${postprocess_err}"
  post_cmd=(sbatch --parsable)
  if [[ ! "${merge_job_id}" =~ ^DRY_RUN_ ]]; then
    post_cmd+=("--dependency=afterok:${merge_job_id%%;*}")
  fi
  post_cmd+=(
    "--job-name=warmup_curve_post"
    "--ntasks=1"
    "--cpus-per-task=${POSTPROCESS_CPUS}"
    "--mem=${POSTPROCESS_MEM}"
    "--time=${POSTPROCESS_TIME}"
  )
  if [[ -n "${QOS}" ]]; then post_cmd+=("--qos=${QOS}"); fi
  if [[ -n "${PARTITION}" ]]; then post_cmd+=("--partition=${PARTITION}"); fi
  if [[ -n "${ACCOUNT}" ]]; then post_cmd+=("--account=${ACCOUNT}"); fi
  post_cmd+=("--output=${postprocess_out}" "--error=${postprocess_err}" "${postprocess_sbatch}")
  echo "Submitting curve regression + summary postprocess:"
  echo "  $(shell_join "${post_cmd[@]}")"
  if truthy "${DRY_RUN}"; then
    postprocess_job_id="DRY_RUN_curve_postprocess"
  else
    postprocess_job_id="$("${post_cmd[@]}")"
  fi
fi

{
  printf 'key\tvalue\n'
	  printf 'project_root\t%s\n' "${PROJECT_ROOT}"
	  printf 'input_root\t%s\n' "${INPUT_ROOT}"
	  printf 'output_root\t%s\n' "${OUTPUT_ROOT}"
	  printf 'seed_manifest\t%s\n' "${SEED_MANIFEST}"
	  printf 'dense_out_dir\t%s\n' "${DENSE_OUT_DIR}"
  printf 'dense_manifest\t%s\n' "${dense_manifest}"
  printf 'merge_job_id\t%s\n' "${merge_job_id}"
  printf 'postprocess_job_id\t%s\n' "${postprocess_job_id}"
  printf 'postprocess_sbatch\t%s\n' "${postprocess_sbatch}"
  printf 'tasks_per_array_task\t%s\n' "${TASKS_PER_ARRAY_TASK}"
  printf 'array_max_concurrent\t%s\n' "${ARRAY_MAX_CONCURRENT}"
} > "${SUMMARY_MANIFEST}"

echo "Dense-grid manifest: ${dense_manifest}"
echo "Merge job id: ${merge_job_id}"
if [[ -n "${postprocess_job_id}" ]]; then echo "Postprocess job id: ${postprocess_job_id}"; fi
echo "Warm-up joint curve array manifest: ${SUMMARY_MANIFEST}"
