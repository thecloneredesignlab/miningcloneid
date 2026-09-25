#!/usr/bin/env bash
# Submit FixO2 eigen-attractor feature computation as a Slurm array workflow.
#
# Default task granularity is one parameter vector per Slurm array task:
#   task = one best-fit or one full-initial parameter vector evaluated across
#          the requested fixed-O2 grid.

set -euo pipefail

O2SD_DOCKER_HPC_ROOT="${O2SD_DOCKER_HPC_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
# shellcheck source=../util/o2_supply_demand_map_apptainer_runtime.sh
source "${O2SD_DOCKER_HPC_ROOT}/util/o2_supply_demand_map_apptainer_runtime.sh"

O2SD_SHELL_UTILS="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../util" && pwd)/o2_supply_demand_map_shell_utils.sh"
# shellcheck source=../../../util/o2_supply_demand_map_shell_utils.sh
source "${O2SD_SHELL_UTILS}"

usage() {
  cat <<'EOF'
Usage:
  bash submit_fixo2_eigen_attractor_array.sh [options]

Primary options:
  --project_root=/share/lab_crd/lab_crd/taoli/Project/miningcloneid_soft_coupling
  --result_root=oxygen/results/analysis/best_fit_parameter_feature/05_FixO2_eigen_attractor_based_clustering
  --source_root=oxygen/results/analysis/best_fit_parameter_feature/02_parameter_landscape_clustering
  --invivo_input=oxygen/results/fit_invivo_O2_buffering_500seed
  --invitro_input=oxygen/results/fit_invitro_O2_buffering_500seed
  --datasets=invivo,invitro
  --point_types=best,initial
  --o2_n=201
  --o2_min=0
  --o2_max=5
  --o2_values=0,1,5              Optional explicit O2 list; overrides o2_n/min/max.
  --points_per_task=1            Default: one parameter vector per Slurm task.
  --array_max_tasks=1000         Split large task counts into multiple sbatch arrays.
  --array_concurrency=500
  --max_seeds=N
  --r_module=R/4.4.2-gfbf-2024a
  --qos=xxlarge
  --partition=NAME
  --account=NAME
  --task_cpus=1
  --task_mem=2G
  --task_time=00:30:00
  --merge_cpus=4
  --merge_mem=32G
  --merge_time=04:00:00
  --run_reductions=TRUE|FALSE
  --reduction_cpus=20
  --reduction_mem=96G
  --reduction_time=12:00:00
  --reductions=pca,umap,tsne
  --run_clustered=TRUE|FALSE
  --dry_run=TRUE|FALSE

Outputs:
  Feature task table:
    RESULT_ROOT/FixO2EigenAttractors/HPC/fixo2_eigen_parameter_tasks.tsv
  Per-task rows:
    RESULT_ROOT/FixO2EigenAttractors/HPC/task_rows/shard_XXXX/task_XXXXXX.csv
  Merged wide tables:
    RESULT_ROOT/FixO2EigenAttractors/Tables/*_fixo2_eigen_attractor_wide.csv
EOF
}

resolve_path() {
  local base="$1"
  local path="$2"
  case "${path}" in
    "~") path="${HOME}" ;;
    "~/"*) path="${HOME}/${path#"~/"}" ;;
  esac
  case "${path}" in
    /*) printf "%s" "${path}" ;;
    *) printf "%s/%s" "${base}" "${path}" ;;
  esac
}

build_sbatch_common_args() {
  COMMON_ARGS=()
  [[ -n "${QOS}" ]] && COMMON_ARGS+=(--qos="${QOS}")
  [[ -n "${PARTITION}" ]] && COMMON_ARGS+=(--partition="${PARTITION}")
  [[ -n "${ACCOUNT}" ]] && COMMON_ARGS+=(--account="${ACCOUNT}")
  return 0
}

submit_or_print() {
  if truthy "${DRY_RUN}"; then
    printf 'DRY_RUN:'
    printf ' %q' "$@"
    printf '\n'
    return 0
  fi
  "$@"
}

DEFAULT_PROJECT_ROOT="/share/lab_crd/lab_crd/taoli/Project/miningcloneid_soft_coupling"
DEFAULT_RESULT_ROOT="oxygen/results/analysis/best_fit_parameter_feature/05_FixO2_eigen_attractor_based_clustering"
DEFAULT_SOURCE_ROOT="oxygen/results/analysis/best_fit_parameter_feature/02_parameter_landscape_clustering"
DEFAULT_INVIVO_INPUT="oxygen/results/fit_invivo_O2_buffering_500seed"
DEFAULT_INVITRO_INPUT="oxygen/results/fit_invitro_O2_buffering_500seed"

PROJECT_ROOT="${PROJECT_ROOT:-${DEFAULT_PROJECT_ROOT}}"
RESULT_ROOT="${RESULT_ROOT:-${DEFAULT_RESULT_ROOT}}"
SOURCE_ROOT="${SOURCE_ROOT:-${DEFAULT_SOURCE_ROOT}}"
INVIVO_INPUT="${INVIVO_INPUT:-${DEFAULT_INVIVO_INPUT}}"
INVITRO_INPUT="${INVITRO_INPUT:-${DEFAULT_INVITRO_INPUT}}"
DATASETS="${DATASETS:-invivo,invitro}"
POINT_TYPES="${POINT_TYPES:-best,initial}"
O2_N="${O2_N:-201}"
O2_MIN="${O2_MIN:-0}"
O2_MAX="${O2_MAX:-5}"
O2_VALUES="${O2_VALUES:-}"
POINTS_PER_TASK="${POINTS_PER_TASK:-1}"
ARRAY_MAX_TASKS="${ARRAY_MAX_TASKS:-1000}"
ARRAY_CONCURRENCY="${ARRAY_CONCURRENCY:-500}"
MAX_SEEDS="${MAX_SEEDS:-}"
R_MODULE="${R_MODULE:-R/4.4.2-gfbf-2024a}"
QOS="${QOS:-xxlarge}"
PARTITION="${PARTITION:-}"
ACCOUNT="${ACCOUNT:-}"
TASK_CPUS="${TASK_CPUS:-1}"
TASK_MEM="${TASK_MEM:-2G}"
TASK_TIME="${TASK_TIME:-00:30:00}"
MERGE_CPUS="${MERGE_CPUS:-4}"
MERGE_MEM="${MERGE_MEM:-32G}"
MERGE_TIME="${MERGE_TIME:-04:00:00}"
RUN_REDUCTIONS="${RUN_REDUCTIONS:-TRUE}"
REDUCTION_CPUS="${REDUCTION_CPUS:-20}"
REDUCTION_MEM="${REDUCTION_MEM:-96G}"
REDUCTION_TIME="${REDUCTION_TIME:-12:00:00}"
REDUCTIONS="${REDUCTIONS:-pca,umap,tsne}"
RUN_CLUSTERED="${RUN_CLUSTERED:-TRUE}"
FORCE_BUILD="${FORCE_BUILD:-FALSE}"
FORCE_TASKS="${FORCE_TASKS:-FALSE}"
FORCE_MERGE="${FORCE_MERGE:-TRUE}"
DRY_RUN="${DRY_RUN:-FALSE}"

for arg in "$@"; do
  case "${arg}" in
    --help|-h) usage; exit 0 ;;
    --project_root=*) PROJECT_ROOT="${arg#*=}" ;;
    --result_root=*) RESULT_ROOT="${arg#*=}" ;;
    --source_root=*) SOURCE_ROOT="${arg#*=}" ;;
    --invivo_input=*) INVIVO_INPUT="${arg#*=}" ;;
    --invitro_input=*) INVITRO_INPUT="${arg#*=}" ;;
    --datasets=*) DATASETS="${arg#*=}" ;;
    --point_types=*) POINT_TYPES="${arg#*=}" ;;
    --o2_n=*) O2_N="${arg#*=}" ;;
    --o2_min=*) O2_MIN="${arg#*=}" ;;
    --o2_max=*) O2_MAX="${arg#*=}" ;;
    --o2_values=*) O2_VALUES="${arg#*=}" ;;
    --points_per_task=*) POINTS_PER_TASK="${arg#*=}" ;;
    --array_max_tasks=*) ARRAY_MAX_TASKS="${arg#*=}" ;;
    --array_concurrency=*) ARRAY_CONCURRENCY="${arg#*=}" ;;
    --max_seeds=*) MAX_SEEDS="${arg#*=}" ;;
    --r_module=*) R_MODULE="${arg#*=}" ;;
    --qos=*) QOS="${arg#*=}" ;;
    --partition=*) PARTITION="${arg#*=}" ;;
    --account=*) ACCOUNT="${arg#*=}" ;;
    --task_cpus=*) TASK_CPUS="${arg#*=}" ;;
    --task_mem=*) TASK_MEM="${arg#*=}" ;;
    --task_time=*) TASK_TIME="${arg#*=}" ;;
    --merge_cpus=*) MERGE_CPUS="${arg#*=}" ;;
    --merge_mem=*) MERGE_MEM="${arg#*=}" ;;
    --merge_time=*) MERGE_TIME="${arg#*=}" ;;
    --run_reductions=*) RUN_REDUCTIONS="${arg#*=}" ;;
    --reduction_cpus=*) REDUCTION_CPUS="${arg#*=}" ;;
    --reduction_mem=*) REDUCTION_MEM="${arg#*=}" ;;
    --reduction_time=*) REDUCTION_TIME="${arg#*=}" ;;
    --reductions=*) REDUCTIONS="${arg#*=}" ;;
    --run_clustered=*) RUN_CLUSTERED="${arg#*=}" ;;
    --force_build=*) FORCE_BUILD="${arg#*=}" ;;
    --force_tasks=*) FORCE_TASKS="${arg#*=}" ;;
    --force_merge=*) FORCE_MERGE="${arg#*=}" ;;
    --dry_run=*) DRY_RUN="${arg#*=}" ;;
    *) echo "Unknown option: ${arg}" >&2; usage >&2; exit 2 ;;
  esac
done

PROJECT_ROOT="$(resolve_path "$(pwd)" "${PROJECT_ROOT}")"
RESULT_ROOT="$(resolve_path "${PROJECT_ROOT}" "${RESULT_ROOT}")"
SOURCE_ROOT="$(resolve_path "${PROJECT_ROOT}" "${SOURCE_ROOT}")"
INVIVO_INPUT="$(resolve_path "${PROJECT_ROOT}" "${INVIVO_INPUT}")"
INVITRO_INPUT="$(resolve_path "${PROJECT_ROOT}" "${INVITRO_INPUT}")"

HPC_TASK_R="${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/Docker/hpc/fixo2_eigen_attractor/fixo2_eigen_attractor_hpc_task.R"
RUNNER_R="${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/runner/fixed_o2_eigen/run_fixo2_eigen_attractor_pipeline.R"
WORKER_SUB="${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/Docker/hpc/array_workers/run_fixo2_eigen_parameter_task.sub"
TASK_TABLE="${RESULT_ROOT}/FixO2EigenAttractors/HPC/fixo2_eigen_parameter_tasks.tsv"
LOG_DIR="${RESULT_ROOT}/logs/hpc_fixo2_eigen"
mkdir -p "${LOG_DIR}"

if [[ ! -f "${HPC_TASK_R}" ]]; then echo "Missing HPC task R script: ${HPC_TASK_R}" >&2; exit 1; fi
if [[ ! -f "${RUNNER_R}" ]]; then echo "Missing runner R script: ${RUNNER_R}" >&2; exit 1; fi
if [[ ! -f "${WORKER_SUB}" ]]; then echo "Missing Slurm worker script: ${WORKER_SUB}" >&2; exit 1; fi

echo "Project root: ${PROJECT_ROOT}"
echo "Result root:  ${RESULT_ROOT}"
echo "Task table:   ${TASK_TABLE}"
echo "Granularity:  ${POINTS_PER_TASK} parameter vector(s) per Slurm task"

if ! truthy "${DRY_RUN}"; then
  cd "${PROJECT_ROOT}"
  load_r_module
  build_args=(
    "${HPC_TASK_R}"
    --mode=build_tasks
    --result_root="${RESULT_ROOT}"
    --source_root="${SOURCE_ROOT}"
    --invivo_input="${INVIVO_INPUT}"
    --invitro_input="${INVITRO_INPUT}"
    --datasets="${DATASETS}"
    --point_types="${POINT_TYPES}"
    --o2_n="${O2_N}"
    --o2_min="${O2_MIN}"
    --o2_max="${O2_MAX}"
    --o2_values="${O2_VALUES}"
    --force="${FORCE_BUILD}"
  )
  [[ -n "${MAX_SEEDS}" ]] && build_args+=(--max_seeds="${MAX_SEEDS}")
  Rscript "${build_args[@]}"
else
  echo "DRY_RUN: would build task table with ${HPC_TASK_R}"
fi

if [[ ! -f "${TASK_TABLE}" ]] && ! truthy "${DRY_RUN}"; then
  echo "Task table was not created: ${TASK_TABLE}" >&2
  exit 1
fi

if truthy "${DRY_RUN}"; then
  TASK_COUNT=1
else
  TASK_LINE_COUNT="$(wc -l < "${TASK_TABLE}")"
  if [[ "${TASK_LINE_COUNT}" -gt 0 ]]; then
    TASK_COUNT=$(( TASK_LINE_COUNT - 1 ))
  else
    TASK_COUNT=0
  fi
fi
if [[ "${TASK_COUNT}" -lt 1 ]]; then echo "No tasks found in ${TASK_TABLE}" >&2; exit 1; fi

SLURM_TASK_COUNT=$(( (TASK_COUNT + POINTS_PER_TASK - 1) / POINTS_PER_TASK ))
ARRAY_MAX_TASKS=$(( ARRAY_MAX_TASKS < 1 ? 1 : ARRAY_MAX_TASKS ))
ARRAY_CONCURRENCY=$(( ARRAY_CONCURRENCY < 1 ? 1 : ARRAY_CONCURRENCY ))

build_sbatch_common_args
array_job_ids=()
array_start=1
while [[ "${array_start}" -le "${SLURM_TASK_COUNT}" ]]; do
  array_end=$(( array_start + ARRAY_MAX_TASKS - 1 ))
  if [[ "${array_end}" -gt "${SLURM_TASK_COUNT}" ]]; then array_end="${SLURM_TASK_COUNT}"; fi
  array_offset=$(( array_start - 1 ))
  chunk_task_count=$(( array_end - array_start + 1 ))
  array_spec="1-${chunk_task_count}%${ARRAY_CONCURRENCY}"
  echo "Submitting FixO2 eigen array ${array_spec} with ARRAY_TASK_OFFSET=${array_offset}"
  if truthy "${DRY_RUN}"; then
    submit_or_print sbatch "${COMMON_ARGS[@]}" --array="${array_spec}" --cpus-per-task="${TASK_CPUS}" --mem="${TASK_MEM}" --time="${TASK_TIME}" "${WORKER_SUB}"
    jid="DRY${array_start}"
  else
    jid="$(sbatch --parsable \
      "${COMMON_ARGS[@]}" \
      --array="${array_spec}" \
      --cpus-per-task="${TASK_CPUS}" \
      --mem="${TASK_MEM}" \
      --time="${TASK_TIME}" \
      --output="${LOG_DIR}/fixo2_eigen_array_%A_%a.out" \
      --error="${LOG_DIR}/fixo2_eigen_array_%A_%a.err" \
      --export=ALL,R_MODULE="${R_MODULE}",PROJECT_ROOT="${PROJECT_ROOT}",RESULT_ROOT="${RESULT_ROOT}",SOURCE_ROOT="${SOURCE_ROOT}",INVIVO_INPUT="${INVIVO_INPUT}",INVITRO_INPUT="${INVITRO_INPUT}",TASK_TABLE="${TASK_TABLE}",HPC_TASK_R="${HPC_TASK_R}",POINTS_PER_TASK="${POINTS_PER_TASK}",MAX_TASK_ID="${TASK_COUNT}",ARRAY_TASK_OFFSET="${array_offset}",O2_N="${O2_N}",O2_MIN="${O2_MIN}",O2_MAX="${O2_MAX}",O2_VALUES="${O2_VALUES}",FORCE_TASKS="${FORCE_TASKS}" \
      "${WORKER_SUB}")"
  fi
  array_job_ids+=("${jid}")
  array_start=$(( array_end + 1 ))
done

dependency=""
if [[ "${#array_job_ids[@]}" -gt 0 ]]; then
  dependency="afterok:$(IFS=':'; echo "${array_job_ids[*]}")"
fi

merge_script="${LOG_DIR}/merge_fixo2_eigen_tasks.sh"
cat > "${merge_script}" <<EOF
#!/usr/bin/env bash
set -euo pipefail
source "${O2SD_DOCKER_HPC_ROOT}/util/o2_supply_demand_map_apptainer_runtime.sh"
cd "${PROJECT_ROOT}"
o2sd_container_ignore_r_module "${R_MODULE}"
Rscript "${HPC_TASK_R}" --mode=merge_tasks --result_root="${RESULT_ROOT}" --source_root="${SOURCE_ROOT}" --invivo_input="${INVIVO_INPUT}" --invitro_input="${INVITRO_INPUT}" --task_table="${TASK_TABLE}" --force="${FORCE_MERGE}"
EOF
chmod +x "${merge_script}"

merge_cmd=(sbatch --parsable "${COMMON_ARGS[@]}" --dependency="${dependency}" --cpus-per-task="${MERGE_CPUS}" --mem="${MERGE_MEM}" --time="${MERGE_TIME}" --job-name=fixo2eig_merge --output="${LOG_DIR}/fixo2_eigen_merge_%j.out" --error="${LOG_DIR}/fixo2_eigen_merge_%j.err" "${merge_script}")
if truthy "${DRY_RUN}"; then
  submit_or_print "${merge_cmd[@]}"
  merge_job_id="DRYMERGE"
else
  merge_job_id="$("${merge_cmd[@]}")"
fi

echo "Merge job id: ${merge_job_id}"

if truthy "${RUN_REDUCTIONS}"; then
  reduction_script="${LOG_DIR}/run_fixo2_eigen_reductions.sh"
  cat > "${reduction_script}" <<EOF
#!/usr/bin/env bash
set -euo pipefail
source "${O2SD_DOCKER_HPC_ROOT}/util/o2_supply_demand_map_apptainer_runtime.sh"
cd "${PROJECT_ROOT}"
o2sd_container_ignore_r_module "${R_MODULE}"
Rscript "${RUNNER_R}" --run_parts=reductions,report --result_root="${RESULT_ROOT}" --reductions="${REDUCTIONS}" --n_threads="${REDUCTION_CPUS}" --run_clustered="${RUN_CLUSTERED}"
EOF
  chmod +x "${reduction_script}"
  reduction_cmd=(sbatch --parsable "${COMMON_ARGS[@]}" --dependency="afterok:${merge_job_id}" --cpus-per-task="${REDUCTION_CPUS}" --mem="${REDUCTION_MEM}" --time="${REDUCTION_TIME}" --job-name=fixo2eig_reduce --output="${LOG_DIR}/fixo2_eigen_reductions_%j.out" --error="${LOG_DIR}/fixo2_eigen_reductions_%j.err" "${reduction_script}")
  if truthy "${DRY_RUN}"; then
    submit_or_print "${reduction_cmd[@]}"
  else
    reduction_job_id="$("${reduction_cmd[@]}")"
    echo "Reduction/report job id: ${reduction_job_id}"
  fi
fi

echo "Array job ids: ${array_job_ids[*]}"
