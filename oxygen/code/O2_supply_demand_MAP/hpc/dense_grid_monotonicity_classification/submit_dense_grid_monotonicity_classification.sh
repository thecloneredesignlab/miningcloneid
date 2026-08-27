#!/usr/bin/env bash
# Submit dense-grid fixed-O2 monotonicity-classification workflows to Slurm.

set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  bash submit_dense_grid_monotonicity_classification.sh [options]

Examples:
  bash submit_dense_grid_monotonicity_classification.sh --run_parts=monotonicity
  bash submit_dense_grid_monotonicity_classification.sh --run_parts=initial_ploidy
  bash submit_dense_grid_monotonicity_classification.sh --run_parts=all

Workflow parts:
  --run_parts=PARTS                  monotonicity, initial_ploidy, or all.
                                     all submits the two workflows independently.

Paths:
  --project_root=DIR                 Repo checkout root on HPC.
  --run_dir=DIR                      500-seed in vivo fit result root.
  --result_root=DIR                  Parent result root. Defaults to
                                     oxygen/results/analysis/best_fit_parameter_feature/03_dense-grid_monotonicity_classification/monotonicity_classification.
  --array_backend=FILE               dense_grid_monotonicity_array_backend.R.
  --log_root=DIR                     Slurm log/script root.

Analysis options:
  --o2_grid=CSV                      O2 grid passed to task-list builder and merge.
  --reporting_o2=CSV                 Reporting O2 values for monotonicity merge.
  --simulation_times=CSV             Initial-ploidy selected/classification times.
                                     Defaults 0,1,2,5,10,20,50,100,200,300,500,700,1000.
  --selected_times=CSV               Backward-compatible alias for --simulation_times.
  --plot_times=CSV                   Initial-ploidy figure times.
                                     Defaults 100,200,300,500,700,1000.
  --time_end=N                       Initial-ploidy full daily endpoint. Defaults 1000.
  --max_seeds=N                      Optional smoke-test seed limit.
  --overwrite=TRUE|FALSE             Rebuild task chunks/final outputs. Defaults TRUE.
  --generate_figures=TRUE|FALSE      Run figure generation in merge jobs. Defaults TRUE.
  --run_validation=TRUE|FALSE        Run workflow validation in merge jobs. Defaults TRUE.

Task splitting:
  --tasks_per_array_task=N           Default chunk size for both workflows.
  --monotonicity_tasks_per_array_task=N
  --initial_tasks_per_array_task=N
  --array_max_concurrent=N           Default Slurm array concurrency cap.
  --monotonicity_array_max_concurrent=N
  --initial_array_max_concurrent=N

Resources:
  --r_module=MODULE                  Defaults R/4.4.2-gfbf-2024a.
  --qos=NAME                         Default QoS.
  --partition=NAME                   Optional Slurm partition.
  --account=NAME                     Optional Slurm account.
  --monotonicity_array_cpus=N        Defaults 1.
  --monotonicity_array_mem=SIZE      Defaults 2G.
  --monotonicity_array_time=TIME     Defaults 1:00:00.
  --monotonicity_merge_cpus=N        Defaults 4. Alias: --monotonicity_cpus.
  --monotonicity_merge_mem=SIZE      Defaults 16G. Alias: --monotonicity_mem.
  --monotonicity_merge_time=TIME     Defaults 2:00:00. Alias: --monotonicity_time.
  --initial_array_cpus=N             Defaults 1.
  --initial_array_mem=SIZE           Defaults 4G.
  --initial_array_time=TIME          Defaults 4:00:00.
  --initial_daily_array_cpus=N       Per-seed daily merge array CPUs. Defaults 1.
  --initial_daily_array_mem=SIZE     Per-seed daily merge array memory. Defaults 2G.
  --initial_daily_array_time=TIME    Per-seed daily merge array time. Defaults 4:00:00.
  --initial_daily_array_max_concurrent=N
                                     Optional per-seed daily merge concurrency cap.
                                     Empty, 0, or none means no Slurm array cap.
  --initial_merge_cpus=N             Defaults 16. Alias: --initial_cpus.
  --initial_merge_mem=SIZE           Defaults 128G. Alias: --initial_mem.
  --initial_merge_time=TIME          Defaults 1-00:00:00. Alias: --initial_time.
  --dependency=SPEC                  Extra Slurm dependency added to each array job.
  --dry_run=TRUE|FALSE               Build task lists and sbatch scripts, but do not submit jobs.
EOF
}

truthy() {
  case "${1:-FALSE}" in
    TRUE|true|True|1|yes|YES|y|Y|on|ON) return 0 ;;
    *) return 1 ;;
  esac
}

shell_join() {
  local out=""
  local token
  local quoted
  for token in "$@"; do
    printf -v quoted "%q" "${token}"
    out+="${quoted} "
  done
  printf "%s" "${out% }"
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

normalize_parts() {
  local raw="${1:-all}"
  local part
  RUN_MONOTONICITY="FALSE"
  RUN_INITIAL_PLOIDY="FALSE"
  IFS=',' read -r -a parts <<< "${raw}"
  for part in "${parts[@]}"; do
    part="$(echo "${part}" | tr '[:upper:]' '[:lower:]' | tr -d '[:space:]')"
    case "${part}" in
      all)
        RUN_MONOTONICITY="TRUE"
        RUN_INITIAL_PLOIDY="TRUE"
        ;;
      monotonicity|dense_grid|dense-grid|classification|ploidy_monotonicity)
        RUN_MONOTONICITY="TRUE"
        ;;
      initial_ploidy|initial-ploidy|trajectory|initial_ploidy_trajectory)
        RUN_INITIAL_PLOIDY="TRUE"
        ;;
      "")
        ;;
      *)
        echo "Unknown run_parts value: ${part}" >&2
        exit 2
        ;;
    esac
  done
  if [[ "${RUN_MONOTONICITY}" != "TRUE" && "${RUN_INITIAL_PLOIDY}" != "TRUE" ]]; then
    echo "No workflow parts selected." >&2
    exit 2
  fi
}

load_r_for_submitter() {
  if [[ -f /etc/profile.d/modules.sh ]]; then source /etc/profile.d/modules.sh; fi
  if command -v module >/dev/null 2>&1; then module use /app/eb/modules/all >/dev/null 2>&1 || true; fi
  if command -v ml >/dev/null 2>&1; then
    ml "${R_MODULE}" >/dev/null 2>&1 || true
  elif command -v module >/dev/null 2>&1; then
    module load "${R_MODULE}" >/dev/null 2>&1 || true
  fi
  RSCRIPT_BIN="$(command -v Rscript || true)"
  if [[ -z "${RSCRIPT_BIN}" ]]; then
    echo "Rscript not found after loading ${R_MODULE}" >&2
    exit 1
  fi
}

append_common_sbatch_headers() {
  local file="$1"
  local job_name="$2"
  local cpus="$3"
  local mem="$4"
  local time_limit="$5"
  local qos="$6"
  local out_file="$7"
  local err_file="$8"
  local array_spec="${9:-}"

  {
    printf '#!/usr/bin/env bash\n'
    printf '#SBATCH --job-name=%s\n' "${job_name}"
    printf '%s\n' '#SBATCH --ntasks=1'
    printf '#SBATCH --cpus-per-task=%s\n' "${cpus}"
    printf '#SBATCH --mem=%s\n' "${mem}"
    printf '#SBATCH --time=%s\n' "${time_limit}"
    if [[ -n "${array_spec}" ]]; then printf '#SBATCH --array=%s\n' "${array_spec}"; fi
    if [[ -n "${qos}" ]]; then printf '#SBATCH --qos=%s\n' "${qos}"; fi
    if [[ -n "${PARTITION}" ]]; then printf '#SBATCH --partition=%s\n' "${PARTITION}"; fi
    if [[ -n "${ACCOUNT}" ]]; then printf '#SBATCH --account=%s\n' "${ACCOUNT}"; fi
    printf '#SBATCH --output=%s\n' "${out_file}"
    printf '#SBATCH --error=%s\n' "${err_file}"
    printf '\n'
    printf 'set -euo pipefail\n'
    printf 'if [[ -f /etc/profile.d/modules.sh ]]; then source /etc/profile.d/modules.sh; fi\n'
    printf 'if command -v module >/dev/null 2>&1; then module use /app/eb/modules/all >/dev/null 2>&1 || true; fi\n'
    printf 'if command -v ml >/dev/null 2>&1; then ml %q; elif command -v module >/dev/null 2>&1; then module load %q; fi\n' "${R_MODULE}" "${R_MODULE}"
    printf 'RSCRIPT_BIN="$(command -v Rscript || true)"\n'
    printf 'if [[ -z "${RSCRIPT_BIN}" ]]; then echo "Rscript not found after loading %s" >&2; exit 1; fi\n' "${R_MODULE}"
    printf 'cd %q\n' "${PROJECT_ROOT}"
    printf 'echo "Host: $(hostname)"\n'
    printf 'echo "Working directory: $(pwd)"\n'
    printf 'echo "Rscript: ${RSCRIPT_BIN}"\n'
    printf '"${RSCRIPT_BIN}" --version\n'
    printf '\n'
  } > "${file}"
}

slurm_array_spec() {
  local n_tasks="$1"
  local concurrency="${2:-}"
  case "${concurrency}" in
    ""|0|none|NONE|None|unlimited|UNLIMITED|Unlimited)
      printf "1-%s" "${n_tasks}"
      ;;
    *)
      printf "1-%s%%%s" "${n_tasks}" "${concurrency}"
      ;;
  esac
}

write_command_block() {
  local file="$1"
  shift
  {
    printf 'cmd=( "${RSCRIPT_BIN}"'
    local arg
    for arg in "$@"; do
      printf ' %q' "${arg}"
    done
    printf ' )\n'
    printf 'printf "Command:"\n'
    printf 'printf " %%q" "${cmd[@]}"\n'
    printf 'printf "\\n"\n'
    printf '"${cmd[@]}"\n'
  } >> "${file}"
}

submit_job() {
  local label="$1"
  local script="$2"
  local dependency="${3:-}"
  shift 3 || true
  local cmd=(sbatch --parsable)
  if [[ -n "${dependency}" ]]; then
    cmd+=(--dependency="${dependency}")
  fi
  if [[ "$#" -gt 0 ]]; then
    cmd+=("$@")
  fi
  cmd+=("${script}")

  printf "%s:" "${label}" | tee -a "${SUBMIT_LOG}" >&2
  printf " %q" "${cmd[@]}" | tee -a "${SUBMIT_LOG}" >&2
  printf "\n" | tee -a "${SUBMIT_LOG}" >&2

  if truthy "${DRY_RUN}"; then
    echo "DRY_RUN_${label}"
  else
    "${cmd[@]}"
  fi
}

task_metadata_value() {
  local metadata="$1"
  local key="$2"
  awk -F '\t' -v k="${key}" '$1 == k {print $2}' "${metadata}"
}

append_optional_common_analysis_args() {
  local -n arr_ref="$1"
  if [[ -n "${O2_GRID}" ]]; then arr_ref+=("--o2_grid=${O2_GRID}"); fi
  if [[ -n "${MAX_SEEDS}" ]]; then arr_ref+=("--max_seeds=${MAX_SEEDS}"); fi
}

parse_args() {
  local arg
  for arg in "$@"; do
    case "${arg}" in
      --help|-h) usage; exit 0 ;;
      --run_parts=*|--workflow_parts=*) RUN_PARTS="${arg#*=}" ;;
      --project_root=*) PROJECT_ROOT="${arg#*=}" ;;
      --run_dir=*|--fit_root=*|--fit_dir=*) RUN_DIR="${arg#*=}" ;;
      --result_root=*|--out_root=*) RESULT_ROOT="${arg#*=}" ;;
      --monotonicity_out_dir=*) MONOTONICITY_OUT_DIR="${arg#*=}" ;;
      --initial_out_dir=*|--initial_output_root=*) INITIAL_OUT_DIR="${arg#*=}" ;;
      --array_backend=*) ARRAY_BACKEND="${arg#*=}" ;;
      --log_root=*) LOG_ROOT="${arg#*=}" ;;
      --o2_grid=*) O2_GRID="${arg#*=}" ;;
      --reporting_o2=*) REPORTING_O2="${arg#*=}" ;;
      --simulation_times=*) SIMULATION_TIMES="${arg#*=}" ;;
      --selected_times=*) SIMULATION_TIMES="${arg#*=}" ;;
      --plot_times=*) PLOT_TIMES="${arg#*=}" ;;
      --time_end=*) TIME_END="${arg#*=}" ;;
      --max_seeds=*) MAX_SEEDS="${arg#*=}" ;;
      --overwrite=*|--force=*) OVERWRITE="${arg#*=}" ;;
      --generate_figures=*) GENERATE_FIGURES="${arg#*=}" ;;
      --run_validation=*) RUN_VALIDATION="${arg#*=}" ;;
      --tasks_per_array_task=*) TASKS_PER_ARRAY_TASK="${arg#*=}" ;;
      --monotonicity_tasks_per_array_task=*) MONOTONICITY_TASKS_PER_ARRAY_TASK="${arg#*=}" ;;
      --initial_tasks_per_array_task=*) INITIAL_TASKS_PER_ARRAY_TASK="${arg#*=}" ;;
      --array_max_concurrent=*) ARRAY_MAX_CONCURRENT="${arg#*=}" ;;
      --monotonicity_array_max_concurrent=*) MONOTONICITY_ARRAY_MAX_CONCURRENT="${arg#*=}" ;;
      --initial_array_max_concurrent=*) INITIAL_ARRAY_MAX_CONCURRENT="${arg#*=}" ;;
      --r_module=*) R_MODULE="${arg#*=}" ;;
      --qos=*) QOS="${arg#*=}" ;;
      --partition=*) PARTITION="${arg#*=}" ;;
      --account=*) ACCOUNT="${arg#*=}" ;;
      --monotonicity_array_cpus=*) MONOTONICITY_ARRAY_CPUS="${arg#*=}" ;;
      --monotonicity_array_mem=*) MONOTONICITY_ARRAY_MEM="${arg#*=}" ;;
      --monotonicity_array_time=*) MONOTONICITY_ARRAY_TIME="${arg#*=}" ;;
      --monotonicity_merge_cpus=*|--monotonicity_cpus=*|--cpus=*) MONOTONICITY_MERGE_CPUS="${arg#*=}" ;;
      --monotonicity_merge_mem=*|--monotonicity_mem=*|--mem=*) MONOTONICITY_MERGE_MEM="${arg#*=}" ;;
      --monotonicity_merge_time=*|--monotonicity_time=*|--time=*|--time_limit=*) MONOTONICITY_MERGE_TIME="${arg#*=}" ;;
      --monotonicity_qos=*) MONOTONICITY_QOS="${arg#*=}" ;;
      --initial_array_cpus=*) INITIAL_ARRAY_CPUS="${arg#*=}" ;;
      --initial_array_mem=*) INITIAL_ARRAY_MEM="${arg#*=}" ;;
      --initial_array_time=*) INITIAL_ARRAY_TIME="${arg#*=}" ;;
      --initial_daily_array_cpus=*) INITIAL_DAILY_ARRAY_CPUS="${arg#*=}" ;;
      --initial_daily_array_mem=*) INITIAL_DAILY_ARRAY_MEM="${arg#*=}" ;;
      --initial_daily_array_time=*) INITIAL_DAILY_ARRAY_TIME="${arg#*=}" ;;
      --initial_daily_array_max_concurrent=*) INITIAL_DAILY_ARRAY_MAX_CONCURRENT="${arg#*=}" ;;
      --initial_merge_cpus=*|--initial_cpus=*) INITIAL_MERGE_CPUS="${arg#*=}" ;;
      --initial_merge_mem=*|--initial_mem=*) INITIAL_MERGE_MEM="${arg#*=}" ;;
      --initial_merge_time=*|--initial_time=*) INITIAL_MERGE_TIME="${arg#*=}" ;;
      --initial_qos=*) INITIAL_QOS="${arg#*=}" ;;
      --dependency=*) EXTRA_DEPENDENCY="${arg#*=}" ;;
      --dry_run=*) DRY_RUN="${arg#*=}" ;;
      *)
        echo "Unknown argument: ${arg}" >&2
        usage >&2
        exit 2
        ;;
    esac
  done
}

build_task_list() {
  local part="$1"
  local out_dir="$2"
  local tasks_per_array_task="$3"
  local task_file="$4"
  local metadata_file="$5"

  local build_cmd=(
    "${ARRAY_BACKEND}"
    "--mode=build_tasks"
    "--part=${part}"
    "--run_dir=${RUN_DIR}"
    "--out_dir=${out_dir}"
    "--tasks_per_array_task=${tasks_per_array_task}"
    "--overwrite=${OVERWRITE}"
  )
  append_optional_common_analysis_args build_cmd
  printf "%s\t%s\t%s\t%s\t%s\t%s\n" "${part}_build_tasks" "local" "" "" "" "$(shell_join "${RSCRIPT_BIN}" "${build_cmd[@]}")" >> "${MANIFEST}"
  printf "Building %s task list:\n  %s\n" "${part}" "$(shell_join "${RSCRIPT_BIN}" "${build_cmd[@]}")" | tee -a "${SUBMIT_LOG}" >&2
  "${RSCRIPT_BIN}" "${build_cmd[@]}"
  if [[ ! -f "${task_file}" || ! -f "${metadata_file}" ]]; then
    echo "Task-list build did not produce expected files: ${task_file}, ${metadata_file}" >&2
    exit 1
  fi
}

submit_workflow_part() {
  local part="$1"
  local out_dir="$2"
  local tasks_per_array_task="$3"
  local array_cpus="$4"
  local array_mem="$5"
  local array_time="$6"
  local array_qos="$7"
  local array_concurrency="$8"
  local merge_cpus="$9"
  local merge_mem="${10}"
  local merge_time="${11}"
  local merge_qos="${12}"

  local task_dir="${out_dir}/hpc/task_lists"
  local task_file="${task_dir}/${part}_seed_o2_tasks.tsv"
  local metadata_file="${task_dir}/${part}_task_metadata.tsv"
  build_task_list "${part}" "${out_dir}" "${tasks_per_array_task}" "${task_file}" "${metadata_file}"

  local n_array_tasks
  n_array_tasks="$(task_metadata_value "${metadata_file}" "n_array_tasks")"
  if [[ -z "${n_array_tasks}" || "${n_array_tasks}" -lt 1 ]]; then
    echo "Invalid n_array_tasks in ${metadata_file}: ${n_array_tasks}" >&2
    exit 1
  fi

  local array_spec
  array_spec="$(slurm_array_spec "${n_array_tasks}" "${array_concurrency}")"
  local label_prefix
  if [[ "${part}" == "initial_ploidy" ]]; then
    label_prefix="initploidy"
  else
    label_prefix="mono"
  fi

  local array_sbatch="${RUN_LOG_DIR}/${label_prefix}_array.sbatch"
  local array_out="${RUN_LOG_DIR}/${label_prefix}_array_%A_%a.out"
  local array_err="${RUN_LOG_DIR}/${label_prefix}_array_%A_%a.err"
  local array_cmd=(
    "${ARRAY_BACKEND}"
    "--mode=run_tasks"
    "--part=${part}"
    "--run_dir=${RUN_DIR}"
    "--out_dir=${out_dir}"
    "--task_file=${task_file}"
    "--tasks_per_array_task=${tasks_per_array_task}"
    "--overwrite=${OVERWRITE}"
  )
  append_optional_common_analysis_args array_cmd
  if [[ "${part}" == "initial_ploidy" ]]; then
    if [[ -n "${SIMULATION_TIMES}" ]]; then array_cmd+=("--simulation_times=${SIMULATION_TIMES}"); fi
    if [[ -n "${TIME_END}" ]]; then array_cmd+=("--time_end=${TIME_END}"); fi
  fi

  append_common_sbatch_headers "${array_sbatch}" "dg_${label_prefix}_arr" "${array_cpus}" "${array_mem}" "${array_time}" "${array_qos}" "${array_out}" "${array_err}" "${array_spec}"
  write_command_block "${array_sbatch}" "${array_cmd[@]}"
  chmod +x "${array_sbatch}"
  local array_sbatch_args=(
    "--job-name=dg_${label_prefix}_arr"
    "--ntasks=1"
    "--array=${array_spec}"
    "--cpus-per-task=${array_cpus}"
    "--mem=${array_mem}"
    "--qos=${array_qos}"
    "--time=${array_time}"
    "--output=${array_out}"
    "--error=${array_err}"
  )
  if [[ -n "${PARTITION}" ]]; then array_sbatch_args+=("--partition=${PARTITION}"); fi
  if [[ -n "${ACCOUNT}" ]]; then array_sbatch_args+=("--account=${ACCOUNT}"); fi
  local array_job_id
  array_job_id="$(submit_job "${part}_array" "${array_sbatch}" "${EXTRA_DEPENDENCY}" "${array_sbatch_args[@]}")"
  printf "%s\t%s\t%s\t%s\t%s\t%s\n" "${part}_array" "${array_job_id}" "${array_sbatch}" "${array_out}" "${array_err}" "$(shell_join Rscript "${array_cmd[@]}")" >> "${MANIFEST}"

  local merge_upstream_job_id="${array_job_id}"
  if [[ "${part}" == "initial_ploidy" ]]; then
    local n_seed
    n_seed="$(task_metadata_value "${metadata_file}" "n_seed")"
    if [[ -z "${n_seed}" || "${n_seed}" -lt 1 ]]; then
      echo "Invalid n_seed in ${metadata_file}: ${n_seed}" >&2
      exit 1
    fi
    local daily_array_spec
    daily_array_spec="$(slurm_array_spec "${n_seed}" "${INITIAL_DAILY_ARRAY_MAX_CONCURRENT}")"
    local daily_sbatch="${RUN_LOG_DIR}/${label_prefix}_daily_seed_array.sbatch"
    local daily_out="${RUN_LOG_DIR}/${label_prefix}_daily_seed_%A_%a.out"
    local daily_err="${RUN_LOG_DIR}/${label_prefix}_daily_seed_%A_%a.err"
    local daily_cmd=(
      "${ARRAY_BACKEND}"
      "--mode=merge_daily_seed"
      "--part=${part}"
      "--run_dir=${RUN_DIR}"
      "--out_dir=${out_dir}"
      "--task_file=${task_file}"
      "--overwrite=${OVERWRITE}"
    )
    append_optional_common_analysis_args daily_cmd
    if [[ -n "${SIMULATION_TIMES}" ]]; then daily_cmd+=("--simulation_times=${SIMULATION_TIMES}"); fi
    if [[ -n "${TIME_END}" ]]; then daily_cmd+=("--time_end=${TIME_END}"); fi

    append_common_sbatch_headers "${daily_sbatch}" "dg_${label_prefix}_daily" "${INITIAL_DAILY_ARRAY_CPUS}" "${INITIAL_DAILY_ARRAY_MEM}" "${INITIAL_DAILY_ARRAY_TIME}" "${array_qos}" "${daily_out}" "${daily_err}" "${daily_array_spec}"
    write_command_block "${daily_sbatch}" "${daily_cmd[@]}"
    chmod +x "${daily_sbatch}"
    local daily_dependency=""
    if [[ ! "${array_job_id}" =~ ^DRY_RUN_ ]]; then
      daily_dependency="afterok:${array_job_id%%;*}"
    fi
    local daily_sbatch_args=(
      "--job-name=dg_${label_prefix}_daily"
      "--ntasks=1"
      "--array=${daily_array_spec}"
      "--cpus-per-task=${INITIAL_DAILY_ARRAY_CPUS}"
      "--mem=${INITIAL_DAILY_ARRAY_MEM}"
      "--qos=${array_qos}"
      "--time=${INITIAL_DAILY_ARRAY_TIME}"
      "--output=${daily_out}"
      "--error=${daily_err}"
    )
    if [[ -n "${PARTITION}" ]]; then daily_sbatch_args+=("--partition=${PARTITION}"); fi
    if [[ -n "${ACCOUNT}" ]]; then daily_sbatch_args+=("--account=${ACCOUNT}"); fi
    local daily_array_job_id
    daily_array_job_id="$(submit_job "${part}_daily_seed_array" "${daily_sbatch}" "${daily_dependency}" "${daily_sbatch_args[@]}")"
    printf "%s\t%s\t%s\t%s\t%s\t%s\n" "${part}_daily_seed_array" "${daily_array_job_id}" "${daily_sbatch}" "${daily_out}" "${daily_err}" "$(shell_join Rscript "${daily_cmd[@]}")" >> "${MANIFEST}"
    merge_upstream_job_id="${daily_array_job_id}"
  fi

  local merge_sbatch="${RUN_LOG_DIR}/${label_prefix}_merge.sbatch"
  local merge_out="${RUN_LOG_DIR}/${label_prefix}_merge_%j.out"
  local merge_err="${RUN_LOG_DIR}/${label_prefix}_merge_%j.err"
  local merge_cmd=(
    "${ARRAY_BACKEND}"
    "--mode=merge"
    "--part=${part}"
    "--run_dir=${RUN_DIR}"
    "--out_dir=${out_dir}"
    "--overwrite=${OVERWRITE}"
    "--generate_figures=${GENERATE_FIGURES}"
    "--run_validation=${RUN_VALIDATION}"
  )
  append_optional_common_analysis_args merge_cmd
  if [[ "${part}" == "monotonicity" ]]; then
    if [[ -n "${REPORTING_O2}" ]]; then merge_cmd+=("--reporting_o2=${REPORTING_O2}"); fi
  else
    if [[ -n "${SIMULATION_TIMES}" ]]; then merge_cmd+=("--simulation_times=${SIMULATION_TIMES}"); fi
    if [[ -n "${PLOT_TIMES}" ]]; then merge_cmd+=("--plot_times=${PLOT_TIMES}"); fi
    if [[ -n "${TIME_END}" ]]; then merge_cmd+=("--time_end=${TIME_END}"); fi
  fi

  append_common_sbatch_headers "${merge_sbatch}" "dg_${label_prefix}_merge" "${merge_cpus}" "${merge_mem}" "${merge_time}" "${merge_qos}" "${merge_out}" "${merge_err}"
  write_command_block "${merge_sbatch}" "${merge_cmd[@]}"
  chmod +x "${merge_sbatch}"
  local merge_dependency=""
  if [[ ! "${merge_upstream_job_id}" =~ ^DRY_RUN_ ]]; then
    merge_dependency="afterok:${merge_upstream_job_id%%;*}"
  fi
  local merge_sbatch_args=(
    "--job-name=dg_${label_prefix}_merge"
    "--ntasks=1"
    "--cpus-per-task=${merge_cpus}"
    "--mem=${merge_mem}"
    "--qos=${merge_qos}"
    "--time=${merge_time}"
    "--output=${merge_out}"
    "--error=${merge_err}"
  )
  if [[ -n "${PARTITION}" ]]; then merge_sbatch_args+=("--partition=${PARTITION}"); fi
  if [[ -n "${ACCOUNT}" ]]; then merge_sbatch_args+=("--account=${ACCOUNT}"); fi
  local merge_job_id
  merge_job_id="$(submit_job "${part}_merge" "${merge_sbatch}" "${merge_dependency}" "${merge_sbatch_args[@]}")"
  printf "%s\t%s\t%s\t%s\t%s\t%s\n" "${part}_merge" "${merge_job_id}" "${merge_sbatch}" "${merge_out}" "${merge_err}" "$(shell_join Rscript "${merge_cmd[@]}")" >> "${MANIFEST}"
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKFLOW_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
DEFAULT_PROJECT_ROOT="$(cd "${WORKFLOW_ROOT}/../../.." && pwd)"

RUN_PARTS="${RUN_PARTS:-all}"
PROJECT_ROOT="${PROJECT_ROOT:-${DEFAULT_PROJECT_ROOT}}"
RESULT_ROOT="${RESULT_ROOT:-oxygen/results/analysis/best_fit_parameter_feature/03_dense-grid_monotonicity_classification/monotonicity_classification}"
RUN_DIR="${RUN_DIR:-oxygen/results/fit_invivo_O2_buffering_500seed}"
MONOTONICITY_OUT_DIR="${MONOTONICITY_OUT_DIR:-}"
INITIAL_OUT_DIR="${INITIAL_OUT_DIR:-}"
ARRAY_BACKEND="${ARRAY_BACKEND:-oxygen/code/O2_supply_demand_MAP/analysis/best_fit_parameter_feature/03_dense-grid_monotonicity_classification/dense_grid_monotonicity_array_backend.R}"
LOG_ROOT="${LOG_ROOT:-}"
O2_GRID="${O2_GRID:-}"
REPORTING_O2="${REPORTING_O2:-}"
SIMULATION_TIMES="${SIMULATION_TIMES:-}"
PLOT_TIMES="${PLOT_TIMES:-}"
TIME_END="${TIME_END:-}"
MAX_SEEDS="${MAX_SEEDS:-}"
OVERWRITE="${OVERWRITE:-TRUE}"
GENERATE_FIGURES="${GENERATE_FIGURES:-TRUE}"
RUN_VALIDATION="${RUN_VALIDATION:-TRUE}"
TASKS_PER_ARRAY_TASK="${TASKS_PER_ARRAY_TASK:-20}"
MONOTONICITY_TASKS_PER_ARRAY_TASK="${MONOTONICITY_TASKS_PER_ARRAY_TASK:-}"
INITIAL_TASKS_PER_ARRAY_TASK="${INITIAL_TASKS_PER_ARRAY_TASK:-}"
ARRAY_MAX_CONCURRENT="${ARRAY_MAX_CONCURRENT:-200}"
MONOTONICITY_ARRAY_MAX_CONCURRENT="${MONOTONICITY_ARRAY_MAX_CONCURRENT:-}"
INITIAL_ARRAY_MAX_CONCURRENT="${INITIAL_ARRAY_MAX_CONCURRENT:-}"
R_MODULE="${R_MODULE:-R/4.4.2-gfbf-2024a}"
QOS="${QOS:-xxlarge}"
PARTITION="${PARTITION:-}"
ACCOUNT="${ACCOUNT:-}"
MONOTONICITY_ARRAY_CPUS="${MONOTONICITY_ARRAY_CPUS:-1}"
MONOTONICITY_ARRAY_MEM="${MONOTONICITY_ARRAY_MEM:-2G}"
MONOTONICITY_ARRAY_TIME="${MONOTONICITY_ARRAY_TIME:-1:00:00}"
MONOTONICITY_MERGE_CPUS="${MONOTONICITY_MERGE_CPUS:-4}"
MONOTONICITY_MERGE_MEM="${MONOTONICITY_MERGE_MEM:-16G}"
MONOTONICITY_MERGE_TIME="${MONOTONICITY_MERGE_TIME:-2:00:00}"
MONOTONICITY_QOS="${MONOTONICITY_QOS:-${QOS}}"
INITIAL_ARRAY_CPUS="${INITIAL_ARRAY_CPUS:-1}"
INITIAL_ARRAY_MEM="${INITIAL_ARRAY_MEM:-4G}"
INITIAL_ARRAY_TIME="${INITIAL_ARRAY_TIME:-4:00:00}"
INITIAL_DAILY_ARRAY_CPUS="${INITIAL_DAILY_ARRAY_CPUS:-1}"
INITIAL_DAILY_ARRAY_MEM="${INITIAL_DAILY_ARRAY_MEM:-2G}"
INITIAL_DAILY_ARRAY_TIME="${INITIAL_DAILY_ARRAY_TIME:-4:00:00}"
INITIAL_DAILY_ARRAY_MAX_CONCURRENT="${INITIAL_DAILY_ARRAY_MAX_CONCURRENT:-}"
INITIAL_MERGE_CPUS="${INITIAL_MERGE_CPUS:-16}"
INITIAL_MERGE_MEM="${INITIAL_MERGE_MEM:-128G}"
INITIAL_MERGE_TIME="${INITIAL_MERGE_TIME:-1-00:00:00}"
INITIAL_QOS="${INITIAL_QOS:-${QOS}}"
EXTRA_DEPENDENCY="${EXTRA_DEPENDENCY:-}"
DRY_RUN="${DRY_RUN:-FALSE}"

parse_args "$@"

MONOTONICITY_TASKS_PER_ARRAY_TASK="${MONOTONICITY_TASKS_PER_ARRAY_TASK:-${TASKS_PER_ARRAY_TASK}}"
INITIAL_TASKS_PER_ARRAY_TASK="${INITIAL_TASKS_PER_ARRAY_TASK:-${TASKS_PER_ARRAY_TASK}}"
MONOTONICITY_ARRAY_MAX_CONCURRENT="${MONOTONICITY_ARRAY_MAX_CONCURRENT:-${ARRAY_MAX_CONCURRENT}}"
INITIAL_ARRAY_MAX_CONCURRENT="${INITIAL_ARRAY_MAX_CONCURRENT:-${ARRAY_MAX_CONCURRENT}}"

PROJECT_ROOT="$(cd "${PROJECT_ROOT}" && pwd)"
RUN_DIR="$(resolve_path "${PROJECT_ROOT}" "${RUN_DIR}")"
RESULT_ROOT="$(resolve_path "${PROJECT_ROOT}" "${RESULT_ROOT}")"
ARRAY_BACKEND="$(resolve_path "${PROJECT_ROOT}" "${ARRAY_BACKEND}")"
if [[ -z "${MONOTONICITY_OUT_DIR}" ]]; then
  MONOTONICITY_OUT_DIR="${RESULT_ROOT}/dense-grid_monotonicity_classification"
else
  MONOTONICITY_OUT_DIR="$(resolve_path "${PROJECT_ROOT}" "${MONOTONICITY_OUT_DIR}")"
fi
if [[ -z "${INITIAL_OUT_DIR}" ]]; then
  INITIAL_OUT_DIR="${RESULT_ROOT}/dense-grid_initial-ploidy_trajectory"
else
  INITIAL_OUT_DIR="$(resolve_path "${PROJECT_ROOT}" "${INITIAL_OUT_DIR}")"
fi
if [[ -z "${LOG_ROOT}" ]]; then
  LOG_ROOT="${RESULT_ROOT}/logs"
else
  LOG_ROOT="$(resolve_path "${PROJECT_ROOT}" "${LOG_ROOT}")"
fi

normalize_parts "${RUN_PARTS}"

if [[ ! -d "${RUN_DIR}" ]]; then
  echo "Missing input run_dir: ${RUN_DIR}" >&2
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

load_r_for_submitter

STAMP="$(date '+%Y%m%d_%H%M%S')"
RUN_LOG_DIR="${LOG_ROOT}/dense_grid_monotonicity_${STAMP}"
mkdir -p "${RUN_LOG_DIR}"
SUBMIT_LOG="${RUN_LOG_DIR}/submit.log"
MANIFEST="${RUN_LOG_DIR}/submitted_jobs.tsv"
: > "${SUBMIT_LOG}"
printf "part\tjob_id\tsbatch_script\tstdout\tstderr\tcommand\n" > "${MANIFEST}"

if [[ "${RUN_MONOTONICITY}" == "TRUE" ]]; then
  submit_workflow_part \
    "monotonicity" \
    "${MONOTONICITY_OUT_DIR}" \
    "${MONOTONICITY_TASKS_PER_ARRAY_TASK}" \
    "${MONOTONICITY_ARRAY_CPUS}" \
    "${MONOTONICITY_ARRAY_MEM}" \
    "${MONOTONICITY_ARRAY_TIME}" \
    "${MONOTONICITY_QOS}" \
    "${MONOTONICITY_ARRAY_MAX_CONCURRENT}" \
    "${MONOTONICITY_MERGE_CPUS}" \
    "${MONOTONICITY_MERGE_MEM}" \
    "${MONOTONICITY_MERGE_TIME}" \
    "${MONOTONICITY_QOS}"
fi

if [[ "${RUN_INITIAL_PLOIDY}" == "TRUE" ]]; then
  submit_workflow_part \
    "initial_ploidy" \
    "${INITIAL_OUT_DIR}" \
    "${INITIAL_TASKS_PER_ARRAY_TASK}" \
    "${INITIAL_ARRAY_CPUS}" \
    "${INITIAL_ARRAY_MEM}" \
    "${INITIAL_ARRAY_TIME}" \
    "${INITIAL_QOS}" \
    "${INITIAL_ARRAY_MAX_CONCURRENT}" \
    "${INITIAL_MERGE_CPUS}" \
    "${INITIAL_MERGE_MEM}" \
    "${INITIAL_MERGE_TIME}" \
    "${INITIAL_QOS}"
fi

echo "Submission manifest: ${MANIFEST}"
if truthy "${DRY_RUN}"; then
  echo "DRY_RUN=TRUE; task lists and sbatch scripts were written, but no jobs were submitted."
fi
