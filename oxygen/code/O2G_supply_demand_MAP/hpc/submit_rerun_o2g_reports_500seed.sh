#!/bin/bash
# Rebuild O2G per-seed viz/report outputs and per-run summary reports on SLURM.
#
# Default run roots:
#   /share/lab_crd/lab_crd/taoli/Project/miningcloneid/oxygen/results/fit_invitro_O2G_buffering_500seed
#   /share/lab_crd/lab_crd/taoli/Project/miningcloneid/oxygen/results/fit_invivo_O2G_buffering_500seed
#   /share/lab_crd/lab_crd/taoli/Project/miningcloneid/oxygen/results/fit_joint_O2G_buffering_500seed
#
# Usage on the HPC login node:
#   bash submit_rerun_o2g_reports_500seed.sh --dry-run
#   bash submit_rerun_o2g_reports_500seed.sh
#
# Optional environment overrides:
#   PROJECT_ROOT=/share/... TOTAL_SEEDS=500 R_MODULE=R/4.4 bash submit_rerun_o2g_reports_500seed.sh

set -euo pipefail

DEFAULT_PROJECT_ROOT="/share/lab_crd/lab_crd/taoli/Project/miningcloneid"
DEFAULT_TOTAL_SEEDS="500"
DEFAULT_SEED_OFFSET="1"
DEFAULT_R_MODULE="R/4.4"
DEFAULT_CPUS_PER_TASK="1"
DEFAULT_MEM="2G"
DEFAULT_TIME_LIMIT="01:00:00"
DEFAULT_REPORT_DT="1"
DEFAULT_TOP_N="6"
DEFAULT_N_CORES="1"
DEFAULT_REPORT_BASENAME="fit_report"
DEFAULT_RENDER_PDF="FALSE"
DEFAULT_NEAR_THRESH="0.05"
DEFAULT_DELETE_LEGACY_REPROT="TRUE"
DEFAULT_DELETE_SUMMARY_EXTRA_RESULTS="TRUE"
DEFAULT_DRY_RUN="FALSE"

INVITRO_RUN_DIR_USER_SET="${INVITRO_RUN_DIR+x}"
INVIVO_RUN_DIR_USER_SET="${INVIVO_RUN_DIR+x}"
JOINT_RUN_DIR_USER_SET="${JOINT_RUN_DIR+x}"
VIZ_INVIVO_SCRIPT_USER_SET="${VIZ_INVIVO_SCRIPT+x}"
VIZ_INVITRO_SCRIPT_USER_SET="${VIZ_INVITRO_SCRIPT+x}"
REPORT_FIT_SCRIPT_USER_SET="${REPORT_FIT_SCRIPT+x}"
REPORT_INVITRO_SCRIPT_USER_SET="${REPORT_INVITRO_SCRIPT+x}"
EXTRA_RESULTS_SCRIPT_USER_SET="${EXTRA_RESULTS_SCRIPT+x}"
DATA_DIR_USER_SET="${DATA_DIR+x}"
OUT_ROOT_USER_SET="${OUT_ROOT+x}"
TASK_ROOT_USER_SET="${TASK_ROOT+x}"
TASK_DIR_USER_SET="${TASK_DIR+x}"

MODE="submit"
if [[ "${1:-}" == "run-seed" || "${1:-}" == "run-summary" || "${1:-}" == "help" || "${1:-}" == "--help" || "${1:-}" == "-h" ]]; then
  MODE="${1}"
  shift || true
fi

PROJECT_ROOT="${PROJECT_ROOT:-${DEFAULT_PROJECT_ROOT}}"
TOTAL_SEEDS="${TOTAL_SEEDS:-${DEFAULT_TOTAL_SEEDS}}"
SEED_OFFSET="${SEED_OFFSET:-${DEFAULT_SEED_OFFSET}}"
R_MODULE="${R_MODULE:-${DEFAULT_R_MODULE}}"
CPUS_PER_TASK="${CPUS_PER_TASK:-${DEFAULT_CPUS_PER_TASK}}"
MEM="${MEM:-${DEFAULT_MEM}}"
TIME_LIMIT="${TIME_LIMIT:-${DEFAULT_TIME_LIMIT}}"
REPORT_DT="${REPORT_DT:-${DEFAULT_REPORT_DT}}"
TOP_N="${TOP_N:-${DEFAULT_TOP_N}}"
N_CORES="${N_CORES:-${DEFAULT_N_CORES}}"
REPORT_BASENAME="${REPORT_BASENAME:-${DEFAULT_REPORT_BASENAME}}"
RENDER_PDF="${RENDER_PDF:-${DEFAULT_RENDER_PDF}}"
NEAR_THRESH="${NEAR_THRESH:-${DEFAULT_NEAR_THRESH}}"
DELETE_LEGACY_REPROT="${DELETE_LEGACY_REPROT:-${DEFAULT_DELETE_LEGACY_REPROT}}"
DELETE_SUMMARY_EXTRA_RESULTS="${DELETE_SUMMARY_EXTRA_RESULTS:-${DEFAULT_DELETE_SUMMARY_EXTRA_RESULTS}}"
DRY_RUN="${DRY_RUN:-${DEFAULT_DRY_RUN}}"

INVITRO_RUN_DIR="${INVITRO_RUN_DIR:-${PROJECT_ROOT}/oxygen/results/fit_invitro_O2G_buffering_500seed}"
INVIVO_RUN_DIR="${INVIVO_RUN_DIR:-${PROJECT_ROOT}/oxygen/results/fit_invivo_O2G_buffering_500seed}"
JOINT_RUN_DIR="${JOINT_RUN_DIR:-${PROJECT_ROOT}/oxygen/results/fit_joint_O2G_buffering_500seed}"

VIZ_INVIVO_SCRIPT="${VIZ_INVIVO_SCRIPT:-${PROJECT_ROOT}/oxygen/code/O2G_supply_demand_MAP/vis/viz_invivo_model_O2G_supply_demand_MAP_results.R}"
VIZ_INVITRO_SCRIPT="${VIZ_INVITRO_SCRIPT:-${PROJECT_ROOT}/oxygen/code/O2G_supply_demand_MAP/vis/viz_invitro_model_O2G_supply_demand_MAP_results.R}"
REPORT_FIT_SCRIPT="${REPORT_FIT_SCRIPT:-${PROJECT_ROOT}/oxygen/code/O2G_supply_demand_MAP/report/render_fit_report.R}"
REPORT_INVITRO_SCRIPT="${REPORT_INVITRO_SCRIPT:-${PROJECT_ROOT}/oxygen/code/O2G_supply_demand_MAP/report/render_invitro_fit_report.R}"
EXTRA_RESULTS_SCRIPT="${EXTRA_RESULTS_SCRIPT:-${PROJECT_ROOT}/oxygen/code/O2G_supply_demand_MAP/analysis/extra_results.R}"
DATA_DIR="${DATA_DIR:-${PROJECT_ROOT}/data/InVivoData_Gemcitabine}"

OUT_ROOT="${OUT_ROOT:-${PROJECT_ROOT}/oxygen/results}"
TASK_ROOT="${TASK_ROOT:-${OUT_ROOT}/report_rerun_tasks}"
STAMP="${STAMP:-$(date '+%Y%m%d_%H%M%S')}"
TASK_DIR="${TASK_DIR:-${TASK_ROOT}/${STAMP}}"

usage() {
  printf "%s\n" \
    "Usage:" \
    "  bash submit_rerun_o2g_reports_500seed.sh [--dry-run]" \
    "" \
    "Internal SLURM modes:" \
    "  submit_rerun_o2g_reports_500seed.sh run-seed <scene> <tasks_tsv>" \
    "  submit_rerun_o2g_reports_500seed.sh run-summary <scene> <run_dir>" \
    "" \
    "Scenes:" \
    "  invitro, invivo, joint" \
    "" \
    "Key overrides:" \
    "  PROJECT_ROOT, INVITRO_RUN_DIR, INVIVO_RUN_DIR, JOINT_RUN_DIR" \
    "  TOTAL_SEEDS, SEED_OFFSET" \
    "  CPUS_PER_TASK, MEM, TIME_LIMIT" \
    "  R_MODULE, DATA_DIR" \
    "  REPORT_DT, TOP_N, N_CORES, REPORT_BASENAME, RENDER_PDF" \
    "  NEAR_THRESH, DELETE_SUMMARY_EXTRA_RESULTS" \
    "" \
    "This script intentionally does not pass --qos to sbatch."
}

truthy() {
  case "${1:-}" in
    TRUE|true|True|1|yes|YES|y|Y|on|ON) return 0 ;;
    *) return 1 ;;
  esac
}

require_positive_integer() {
  local name="$1"
  local value="$2"
  if ! [[ "${value}" =~ ^[0-9]+$ ]]; then
    echo "${name} must be a positive integer, got: ${value}" >&2
    exit 1
  fi
  if (( value <= 0 )); then
    echo "${name} must be > 0, got: ${value}" >&2
    exit 1
  fi
}

validate_common_numeric() {
  require_positive_integer "TOTAL_SEEDS" "${TOTAL_SEEDS}"
  require_positive_integer "SEED_OFFSET" "${SEED_OFFSET}"
  require_positive_integer "CPUS_PER_TASK" "${CPUS_PER_TASK}"
  require_positive_integer "N_CORES" "${N_CORES}"
  require_positive_integer "TOP_N" "${TOP_N}"
}

script_path() {
  local src
  src="${BASH_SOURCE[0]}"
  (
    cd "$(dirname "${src}")" >/dev/null 2>&1
    printf "%s/%s\n" "$(pwd)" "$(basename "${src}")"
  )
}

apply_project_root_defaults() {
  if [[ -z "${INVITRO_RUN_DIR_USER_SET}" ]]; then
    INVITRO_RUN_DIR="${PROJECT_ROOT}/oxygen/results/fit_invitro_O2G_buffering_500seed"
  fi
  if [[ -z "${INVIVO_RUN_DIR_USER_SET}" ]]; then
    INVIVO_RUN_DIR="${PROJECT_ROOT}/oxygen/results/fit_invivo_O2G_buffering_500seed"
  fi
  if [[ -z "${JOINT_RUN_DIR_USER_SET}" ]]; then
    JOINT_RUN_DIR="${PROJECT_ROOT}/oxygen/results/fit_joint_O2G_buffering_500seed"
  fi
  if [[ -z "${VIZ_INVIVO_SCRIPT_USER_SET}" ]]; then
    VIZ_INVIVO_SCRIPT="${PROJECT_ROOT}/oxygen/code/O2G_supply_demand_MAP/vis/viz_invivo_model_O2G_supply_demand_MAP_results.R"
  fi
  if [[ -z "${VIZ_INVITRO_SCRIPT_USER_SET}" ]]; then
    VIZ_INVITRO_SCRIPT="${PROJECT_ROOT}/oxygen/code/O2G_supply_demand_MAP/vis/viz_invitro_model_O2G_supply_demand_MAP_results.R"
  fi
  if [[ -z "${REPORT_FIT_SCRIPT_USER_SET}" ]]; then
    REPORT_FIT_SCRIPT="${PROJECT_ROOT}/oxygen/code/O2G_supply_demand_MAP/report/render_fit_report.R"
  fi
  if [[ -z "${REPORT_INVITRO_SCRIPT_USER_SET}" ]]; then
    REPORT_INVITRO_SCRIPT="${PROJECT_ROOT}/oxygen/code/O2G_supply_demand_MAP/report/render_invitro_fit_report.R"
  fi
  if [[ -z "${EXTRA_RESULTS_SCRIPT_USER_SET}" ]]; then
    EXTRA_RESULTS_SCRIPT="${PROJECT_ROOT}/oxygen/code/O2G_supply_demand_MAP/analysis/extra_results.R"
  fi
  if [[ -z "${DATA_DIR_USER_SET}" ]]; then
    DATA_DIR="${PROJECT_ROOT}/data/InVivoData_Gemcitabine"
  fi
  if [[ -z "${OUT_ROOT_USER_SET}" ]]; then
    OUT_ROOT="${PROJECT_ROOT}/oxygen/results"
  fi
  if [[ -z "${TASK_ROOT_USER_SET}" ]]; then
    TASK_ROOT="${OUT_ROOT}/report_rerun_tasks"
  fi
  if [[ -z "${TASK_DIR_USER_SET}" ]]; then
    TASK_DIR="${TASK_ROOT}/${STAMP}"
  fi
}

load_r_module() {
  if command -v ml >/dev/null 2>&1; then
    ml "${R_MODULE}"
  elif command -v module >/dev/null 2>&1; then
    module load "${R_MODULE}"
  fi
  if ! command -v Rscript >/dev/null 2>&1; then
    echo "Rscript not found after loading module: ${R_MODULE}" >&2
    exit 1
  fi
}

set_thread_limits() {
  export OMP_NUM_THREADS=1
  export OPENBLAS_NUM_THREADS=1
  export MKL_NUM_THREADS=1
  export VECLIB_MAXIMUM_THREADS=1
}

scene_run_dir() {
  case "$1" in
    invitro) printf "%s\n" "${INVITRO_RUN_DIR}" ;;
    invivo) printf "%s\n" "${INVIVO_RUN_DIR}" ;;
    joint) printf "%s\n" "${JOINT_RUN_DIR}" ;;
    *) echo "Unknown scene: $1" >&2; exit 1 ;;
  esac
}

scene_required_files() {
  case "$1" in
    invitro)
      printf "%s\n" fit_summary.tsv best_params.tsv fit_result.rds
      ;;
    invivo)
      printf "%s\n" fit_summary.tsv best_params.tsv fit_config.rds burden_fit.tsv terminal_ploidy_fit.tsv
      ;;
    joint)
      printf "%s\n" fit_summary.tsv best_params.tsv fit_config.rds fit_result.rds joint_components.tsv
      ;;
    *)
      echo "Unknown scene: $1" >&2
      exit 1
      ;;
  esac
}

missing_required_files() {
  local scene="$1"
  local seed_dir="$2"
  local missing=()
  local file_name
  while IFS= read -r file_name; do
    if [[ ! -f "${seed_dir}/${file_name}" ]]; then
      missing+=("${file_name}")
    fi
  done < <(scene_required_files "${scene}")
  if (( ${#missing[@]} > 0 )); then
    local IFS=","
    printf "%s\n" "${missing[*]}"
  fi
}

seed_is_complete() {
  local scene="$1"
  local seed_dir="$2"
  [[ -d "${seed_dir}" ]] || return 1
  [[ -z "$(missing_required_files "${scene}" "${seed_dir}")" ]]
}

validate_submit_paths() {
  local path
  for path in \
    "${VIZ_INVIVO_SCRIPT}" \
    "${VIZ_INVITRO_SCRIPT}" \
    "${REPORT_FIT_SCRIPT}" \
    "${REPORT_INVITRO_SCRIPT}" \
    "${EXTRA_RESULTS_SCRIPT}"; do
    if [[ ! -f "${path}" ]]; then
      echo "Missing required script: ${path}" >&2
      exit 1
    fi
  done

  if [[ ! -d "${DATA_DIR}" ]]; then
    echo "Missing in vivo data directory: ${DATA_DIR}" >&2
    exit 1
  fi

  for path in "${INVITRO_RUN_DIR}" "${INVIVO_RUN_DIR}" "${JOINT_RUN_DIR}"; do
    if [[ ! -d "${path}" ]]; then
      echo "Missing run directory: ${path}" >&2
      exit 1
    fi
  done
}

build_scene_tasks() {
  local scene="$1"
  local run_dir="$2"
  local tasks_tsv="$3"
  local skipped_tsv="$4"
  local seed
  local seed_dir
  local task_id=0
  local missing

  {
    printf "task_id\tscene\tseed\tseed_dir\n"
  } > "${tasks_tsv}"
  {
    printf "scene\tseed\tseed_dir\tmissing_required_files\n"
  } > "${skipped_tsv}"

  for (( seed = SEED_OFFSET; seed < SEED_OFFSET + TOTAL_SEEDS; seed++ )); do
    seed_dir="${run_dir}/seed${seed}"
    if seed_is_complete "${scene}" "${seed_dir}"; then
      task_id=$((task_id + 1))
      printf "%s\t%s\t%s\t%s\n" "${task_id}" "${scene}" "${seed}" "${seed_dir}" >> "${tasks_tsv}"
    else
      missing="$(missing_required_files "${scene}" "${seed_dir}")"
      if [[ -z "${missing}" && ! -d "${seed_dir}" ]]; then
        missing="seed_dir"
      fi
      printf "%s\t%s\t%s\t%s\n" "${scene}" "${seed}" "${seed_dir}" "${missing}" >> "${skipped_tsv}"
    fi
  done

  printf "%s\n" "${task_id}"
}

task_count_from_tsv() {
  local tasks_tsv="$1"
  local n
  n="$(awk 'NR > 1 {n++} END {print n + 0}' "${tasks_tsv}")"
  printf "%s\n" "${n}"
}

export_arg_common() {
  local export_arg="ALL"
  export_arg+=",PROJECT_ROOT=${PROJECT_ROOT}"
  export_arg+=",INVITRO_RUN_DIR=${INVITRO_RUN_DIR}"
  export_arg+=",INVIVO_RUN_DIR=${INVIVO_RUN_DIR}"
  export_arg+=",JOINT_RUN_DIR=${JOINT_RUN_DIR}"
  export_arg+=",VIZ_INVIVO_SCRIPT=${VIZ_INVIVO_SCRIPT}"
  export_arg+=",VIZ_INVITRO_SCRIPT=${VIZ_INVITRO_SCRIPT}"
  export_arg+=",REPORT_FIT_SCRIPT=${REPORT_FIT_SCRIPT}"
  export_arg+=",REPORT_INVITRO_SCRIPT=${REPORT_INVITRO_SCRIPT}"
  export_arg+=",EXTRA_RESULTS_SCRIPT=${EXTRA_RESULTS_SCRIPT}"
  export_arg+=",DATA_DIR=${DATA_DIR}"
  export_arg+=",R_MODULE=${R_MODULE}"
  export_arg+=",REPORT_DT=${REPORT_DT}"
  export_arg+=",TOP_N=${TOP_N}"
  export_arg+=",N_CORES=${N_CORES}"
  export_arg+=",REPORT_BASENAME=${REPORT_BASENAME}"
  export_arg+=",RENDER_PDF=${RENDER_PDF}"
  export_arg+=",NEAR_THRESH=${NEAR_THRESH}"
  export_arg+=",DELETE_LEGACY_REPROT=${DELETE_LEGACY_REPROT}"
  export_arg+=",DELETE_SUMMARY_EXTRA_RESULTS=${DELETE_SUMMARY_EXTRA_RESULTS}"
  printf "%s\n" "${export_arg}"
}

submit_or_print() {
  if truthy "${DRY_RUN}"; then
    printf "DRY-RUN:"
    printf " %q" "$@"
    printf "\n"
  else
    "$@"
  fi
}

run_to_log() {
  local label="$1"
  local log_path="$2"
  shift 2

  {
    echo "[$(date '+%F %T')] ${label} start"
    printf "Command:"
    printf " %q" "$@"
    printf "\n"
  } > "${log_path}"

  echo "[$(date '+%F %T')] ${label} start"
  if "$@" >> "${log_path}" 2>&1; then
    echo "[$(date '+%F %T')] ${label} done"
  else
    local status=$?
    echo "[$(date '+%F %T')] ${label} failed with exit status ${status}" >&2
    tail -n 80 "${log_path}" >&2 || true
    exit "${status}"
  fi
}

run_seed_task() {
  local scene="$1"
  local tasks_tsv="$2"
  local task_id="${SLURM_ARRAY_TASK_ID:-1}"
  local line
  local row_task_id
  local row_scene
  local seed
  local seed_dir
  local status_log
  local missing

  validate_common_numeric
  load_r_module
  set_thread_limits

  if [[ ! -f "${tasks_tsv}" ]]; then
    echo "Missing task table: ${tasks_tsv}" >&2
    exit 1
  fi

  line="$(awk -F '\t' -v id="${task_id}" 'NR > 1 && $1 == id {print; found=1; exit} END {if (!found) exit 2}' "${tasks_tsv}")" || {
    echo "No task row found for SLURM_ARRAY_TASK_ID=${task_id} in ${tasks_tsv}" >&2
    exit 1
  }
  IFS=$'\t' read -r row_task_id row_scene seed seed_dir <<< "${line}"

  if [[ "${row_scene}" != "${scene}" ]]; then
    echo "Task scene mismatch: requested ${scene}, row has ${row_scene}" >&2
    exit 1
  fi

  missing="$(missing_required_files "${scene}" "${seed_dir}")"
  if [[ -n "${missing}" ]]; then
    echo "Skipping seed${seed}; missing required files after submission: ${missing}"
    exit 0
  fi

  status_log="${seed_dir}/report_rerun_status.log"
  : > "${status_log}"
  {
    echo "[$(date '+%F %T')] scene=${scene}"
    echo "[$(date '+%F %T')] seed=${seed}"
    echo "[$(date '+%F %T')] seed_dir=${seed_dir}"
    echo "[$(date '+%F %T')] slurm_job_id=${SLURM_JOB_ID:-NA}"
    echo "[$(date '+%F %T')] slurm_array_task_id=${SLURM_ARRAY_TASK_ID:-NA}"
    echo "[$(date '+%F %T')] deleting old viz/report"
  } | tee -a "${status_log}"

  rm -rf "${seed_dir}/viz" "${seed_dir}/report"
  if truthy "${DELETE_LEGACY_REPROT}"; then
    rm -rf "${seed_dir}/reprot"
  fi

  case "${scene}" in
    invitro)
      run_to_log "in vitro viz" "${seed_dir}/invitro_viz_status.log" \
        Rscript "${VIZ_INVITRO_SCRIPT}" \
        "--fit_dir=${seed_dir}" \
        "--out_dir=${seed_dir}/viz"
      run_to_log "in vitro report" "${seed_dir}/report_status.log" \
        Rscript "${REPORT_INVITRO_SCRIPT}" \
        "--fit_dir=${seed_dir}" \
        "--out_subdir=report" \
        "--report_basename=${REPORT_BASENAME}"
      ;;
    invivo)
      run_to_log "in vivo viz" "${seed_dir}/viz_status.log" \
        Rscript "${VIZ_INVIVO_SCRIPT}" \
        "--fit_dir=${seed_dir}" \
        "--out_dir=${seed_dir}/viz" \
        "--data_dir=${DATA_DIR}" \
        "--report_dt=${REPORT_DT}" \
        "--top_n=${TOP_N}" \
        "--n_cores=${N_CORES}"
      run_to_log "in vivo report" "${seed_dir}/report_status.log" \
        Rscript "${REPORT_FIT_SCRIPT}" \
        "--fit_dir=${seed_dir}" \
        "--out_subdir=report" \
        "--report_basename=${REPORT_BASENAME}" \
        "--render_pdf=${RENDER_PDF}"
      ;;
    joint)
      run_to_log "joint in vivo viz" "${seed_dir}/viz_status.log" \
        Rscript "${VIZ_INVIVO_SCRIPT}" \
        "--fit_dir=${seed_dir}" \
        "--out_dir=${seed_dir}/viz/invivo" \
        "--data_dir=${DATA_DIR}" \
        "--report_dt=${REPORT_DT}" \
        "--top_n=${TOP_N}" \
        "--n_cores=${N_CORES}"
      run_to_log "joint in vitro viz" "${seed_dir}/invitro_viz_status.log" \
        Rscript "${VIZ_INVITRO_SCRIPT}" \
        "--fit_dir=${seed_dir}" \
        "--out_dir=${seed_dir}/viz/invitro"
      run_to_log "joint report" "${seed_dir}/report_status.log" \
        Rscript "${REPORT_FIT_SCRIPT}" \
        "--fit_dir=${seed_dir}" \
        "--out_subdir=report" \
        "--report_basename=${REPORT_BASENAME}" \
        "--render_pdf=${RENDER_PDF}"
      ;;
    *)
      echo "Unknown scene: ${scene}" >&2
      exit 1
      ;;
  esac

  echo "[$(date '+%F %T')] seed task completed" | tee -a "${status_log}"
}

run_summary_task() {
  local scene="$1"
  local run_dir="$2"
  local out_dir="${run_dir}/extra_results"
  local log_path="${run_dir}/summary_report_rerun_status.log"

  validate_common_numeric
  load_r_module
  set_thread_limits

  if [[ ! -d "${run_dir}" ]]; then
    echo "Missing run directory: ${run_dir}" >&2
    exit 1
  fi
  if [[ ! -f "${EXTRA_RESULTS_SCRIPT}" ]]; then
    echo "Missing extra results script: ${EXTRA_RESULTS_SCRIPT}" >&2
    exit 1
  fi

  : > "${log_path}"
  {
    echo "[$(date '+%F %T')] scene=${scene}"
    echo "[$(date '+%F %T')] run_dir=${run_dir}"
    echo "[$(date '+%F %T')] out_dir=${out_dir}"
    echo "[$(date '+%F %T')] slurm_job_id=${SLURM_JOB_ID:-NA}"
  } | tee -a "${log_path}"

  if truthy "${DELETE_SUMMARY_EXTRA_RESULTS}"; then
    echo "[$(date '+%F %T')] deleting old summary extra_results directory" | tee -a "${log_path}"
    rm -rf "${out_dir}"
  fi

  run_to_log "${scene} summary report" "${run_dir}/extra_results_status.log" \
    Rscript "${EXTRA_RESULTS_SCRIPT}" \
    "--run_dir=${run_dir}" \
    "--out_dir=${out_dir}" \
    "--near_thresh=${NEAR_THRESH}"

  echo "[$(date '+%F %T')] summary task completed" | tee -a "${log_path}"
}

parse_submit_args() {
  while (( $# > 0 )); do
    case "$1" in
      --dry-run)
        DRY_RUN="TRUE"
        shift
        ;;
      --project-root=*)
        PROJECT_ROOT="${1#*=}"
        shift
        ;;
      --task-dir=*)
        TASK_DIR="${1#*=}"
        TASK_DIR_USER_SET="x"
        shift
        ;;
      --total-seeds=*)
        TOTAL_SEEDS="${1#*=}"
        shift
        ;;
      --seed-offset=*)
        SEED_OFFSET="${1#*=}"
        shift
        ;;
      --help|-h|help)
        usage
        exit 0
        ;;
      *)
        echo "Unknown argument: $1" >&2
        usage >&2
        exit 1
        ;;
    esac
  done
}

submit_scene() {
  local scene="$1"
  local run_dir="$2"
  local tasks_tsv="${TASK_DIR}/${scene}_per_seed_tasks.tsv"
  local skipped_tsv="${TASK_DIR}/${scene}_skipped_seeds.tsv"
  local task_count
  local script
  local export_arg
  local array_job_id=""
  local summary_job_id=""
  local output_dir="${TASK_DIR}/slurm_logs"
  local dependency_arg=()

  task_count="$(build_scene_tasks "${scene}" "${run_dir}" "${tasks_tsv}" "${skipped_tsv}")"
  echo "Built ${scene} tasks: ${task_count} per-seed task(s)"
  echo "  tasks: ${tasks_tsv}"
  echo "  skipped: ${skipped_tsv}"

  if (( task_count <= 0 )); then
    echo "No complete ${scene} seeds found; skipping ${scene} per-seed and summary submissions."
    return 0
  fi

  script="$(script_path)"
  export_arg="$(export_arg_common)"

  if truthy "${DRY_RUN}"; then
    submit_or_print \
      sbatch \
      --parsable \
      "--job-name=o2g_${scene}_rpt" \
      "--time=${TIME_LIMIT}" \
      "--mem=${MEM}" \
      "--cpus-per-task=${CPUS_PER_TASK}" \
      "--array=1-${task_count}" \
      "--output=${output_dir}/o2g_${scene}_rpt_%A_%a.out" \
      "--error=${output_dir}/o2g_${scene}_rpt_%A_%a.err" \
      "--export=${export_arg}" \
      "${script}" run-seed "${scene}" "${tasks_tsv}"
    submit_or_print \
      sbatch \
      --parsable \
      "--job-name=o2g_${scene}_sum" \
      "--time=${TIME_LIMIT}" \
      "--mem=${MEM}" \
      "--cpus-per-task=${CPUS_PER_TASK}" \
      "--output=${output_dir}/o2g_${scene}_sum_%j.out" \
      "--error=${output_dir}/o2g_${scene}_sum_%j.err" \
      "--export=${export_arg}" \
      "${script}" run-summary "${scene}" "${run_dir}"
    return 0
  fi

  array_job_id="$(
    sbatch \
      --parsable \
      "--job-name=o2g_${scene}_rpt" \
      "--time=${TIME_LIMIT}" \
      "--mem=${MEM}" \
      "--cpus-per-task=${CPUS_PER_TASK}" \
      "--array=1-${task_count}" \
      "--output=${output_dir}/o2g_${scene}_rpt_%A_%a.out" \
      "--error=${output_dir}/o2g_${scene}_rpt_%A_%a.err" \
      "--export=${export_arg}" \
      "${script}" run-seed "${scene}" "${tasks_tsv}"
  )"
  array_job_id="${array_job_id%%;*}"
  echo "Submitted ${scene} per-seed array job: ${array_job_id}"

  dependency_arg=("--dependency=afterok:${array_job_id}")
  summary_job_id="$(
    sbatch \
      --parsable \
      "${dependency_arg[@]}" \
      "--job-name=o2g_${scene}_sum" \
      "--time=${TIME_LIMIT}" \
      "--mem=${MEM}" \
      "--cpus-per-task=${CPUS_PER_TASK}" \
      "--output=${output_dir}/o2g_${scene}_sum_%j.out" \
      "--error=${output_dir}/o2g_${scene}_sum_%j.err" \
      "--export=${export_arg}" \
      "${script}" run-summary "${scene}" "${run_dir}"
  )"
  summary_job_id="${summary_job_id%%;*}"
  echo "Submitted ${scene} summary job: ${summary_job_id} afterok:${array_job_id}"

  printf "%s\t%s\t%s\t%s\n" "${scene}" "${run_dir}" "${array_job_id}" "${summary_job_id}" >> "${TASK_DIR}/submitted_jobs.tsv"
}

submit_main() {
  parse_submit_args "$@"
  apply_project_root_defaults
  validate_common_numeric
  validate_submit_paths

  mkdir -p "${TASK_DIR}/slurm_logs"
  {
    printf "scene\trun_dir\tper_seed_array_job_id\tsummary_job_id\n"
  } > "${TASK_DIR}/submitted_jobs.tsv"

  echo "Submitting O2G report rerun tasks"
  echo "  project_root: ${PROJECT_ROOT}"
  echo "  task_dir: ${TASK_DIR}"
  echo "  total_seeds: ${TOTAL_SEEDS}"
  echo "  seed_offset: ${SEED_OFFSET}"
  echo "  resources: cpus=${CPUS_PER_TASK}, mem=${MEM}, time=${TIME_LIMIT}"
  echo "  qos: default (no --qos passed)"
  echo "  dry_run: ${DRY_RUN}"

  if ! truthy "${DRY_RUN}" && ! command -v sbatch >/dev/null 2>&1; then
    echo "sbatch not found; run this launcher on the HPC login node." >&2
    exit 1
  fi

  submit_scene "invitro" "${INVITRO_RUN_DIR}"
  submit_scene "invivo" "${INVIVO_RUN_DIR}"
  submit_scene "joint" "${JOINT_RUN_DIR}"

  echo "Submission manifest: ${TASK_DIR}/submitted_jobs.tsv"
}

case "${MODE}" in
  submit)
    submit_main "$@"
    ;;
  run-seed)
    if (( $# != 2 )); then
      echo "Usage: $0 run-seed <scene> <tasks_tsv>" >&2
      exit 1
    fi
    run_seed_task "$1" "$2"
    ;;
  run-summary)
    if (( $# != 2 )); then
      echo "Usage: $0 run-summary <scene> <run_dir>" >&2
      exit 1
    fi
    run_summary_task "$1" "$2"
    ;;
  help|--help|-h)
    usage
    ;;
  *)
    echo "Unknown mode: ${MODE}" >&2
    usage >&2
    exit 1
    ;;
esac
