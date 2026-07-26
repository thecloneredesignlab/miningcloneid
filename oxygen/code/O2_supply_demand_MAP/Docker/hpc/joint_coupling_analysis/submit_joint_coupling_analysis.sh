#!/usr/bin/env bash

# Submit the complete joint soft-coupling and ploidy-category analysis as one
# Slurm job. The executable worker remains in this hpc/ feature folder.

set -euo pipefail

O2SD_DOCKER_HPC_ROOT="${O2SD_DOCKER_HPC_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
# shellcheck source=../util/o2_supply_demand_map_apptainer_runtime.sh
source "${O2SD_DOCKER_HPC_ROOT}/util/o2_supply_demand_map_apptainer_runtime.sh"

PROJECT_ROOT="${PROJECT_ROOT:-/share/lab_crd/lab_crd/taoli/Project/miningcloneid_soft_coupling}"
RESULT_ROOT="${RESULT_ROOT:-/share/lab_crd/lab_crd/taoli/Project/miningcloneid_soft_coupling/oxygen/results/fit_joint_multi_warmup_tsne_sigmaN0p1216_objmin_500seed_20260714_021540}"
OUTPUT_ROOT="${OUTPUT_ROOT:-${PROJECT_ROOT}/oxygen/results/analysis/joint_coupling/$(basename "${RESULT_ROOT%/}")}"
FIXED_O2_SOURCE_ROOT="${FIXED_O2_SOURCE_ROOT:-}"
R_MODULE="${R_MODULE:-R/4.4.2-gfbf-2024a}"
CLASS_THRESHOLD="${CLASS_THRESHOLD:-1.1}"
CLASS_LOWER_BOUND="${CLASS_LOWER_BOUND:-}"
CLASS_UPPER_BOUND="${CLASS_UPPER_BOUND:-}"
CLASS_BOUNDARY_RULE="${CLASS_BOUNDARY_RULE:-classb_inclusive}"
N_BOOT="${N_BOOT:-2000}"
PERMUTATIONS="${PERMUTATIONS:-999}"
TRAJECTORY_STEP="${TRAJECTORY_STEP:-10}"
MAX_PAIRS="${MAX_PAIRS:-}"
MAX_SEEDS="${MAX_SEEDS:-}"
QOS="${QOS:-xxlarge}"
PARTITION="${PARTITION:-}"
ACCOUNT="${ACCOUNT:-}"
CPUS="${CPUS:-4}"
MEMORY="${MEMORY:-64G}"
TIME_LIMIT="${TIME_LIMIT:-12:00:00}"
DRY_RUN="${DRY_RUN:-FALSE}"

for arg in "$@"; do
  case "${arg}" in
    --project_root=*) PROJECT_ROOT="${arg#*=}" ;;
    --result_root=*) RESULT_ROOT="${arg#*=}" ;;
    --output_root=*) OUTPUT_ROOT="${arg#*=}" ;;
    --fixed_o2_source_root=*) FIXED_O2_SOURCE_ROOT="${arg#*=}" ;;
    --r_module=*) R_MODULE="${arg#*=}" ;;
    --class_threshold=*) CLASS_THRESHOLD="${arg#*=}" ;;
    --class_lower_bound=*) CLASS_LOWER_BOUND="${arg#*=}" ;;
    --class_upper_bound=*) CLASS_UPPER_BOUND="${arg#*=}" ;;
    --class_boundary_rule=*) CLASS_BOUNDARY_RULE="${arg#*=}" ;;
    --n_boot=*) N_BOOT="${arg#*=}" ;;
    --permutations=*) PERMUTATIONS="${arg#*=}" ;;
    --trajectory_step=*) TRAJECTORY_STEP="${arg#*=}" ;;
    --max_pairs=*) MAX_PAIRS="${arg#*=}" ;;
    --max_seeds=*) MAX_SEEDS="${arg#*=}" ;;
    --qos=*) QOS="${arg#*=}" ;;
    --partition=*) PARTITION="${arg#*=}" ;;
    --account=*) ACCOUNT="${arg#*=}" ;;
    --cpus=*) CPUS="${arg#*=}" ;;
    --memory=*) MEMORY="${arg#*=}" ;;
    --time=*) TIME_LIMIT="${arg#*=}" ;;
    --dry_run=*) DRY_RUN="${arg#*=}" ;;
    --help|-h)
      echo "Usage: bash submit_joint_coupling_analysis.sh [--result_root=PATH] [--output_root=PATH] [--fixed_o2_source_root=PATH] [--class_threshold=N | --class_lower_bound=N --class_upper_bound=N --class_boundary_rule=RULE] [--max_pairs=N] [--max_seeds=N] [--dry_run=TRUE]"
      exit 0
      ;;
    *) echo "Unknown option: ${arg}" >&2; exit 2 ;;
  esac
done

worker="${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/Docker/hpc/joint_coupling_analysis/run_joint_coupling_analysis.sub"
if [[ "${OUTPUT_ROOT%/}" == "${RESULT_ROOT%/}" || "${OUTPUT_ROOT%/}" == "${RESULT_ROOT%/}"/* ]]; then
  echo "ERROR: OUTPUT_ROOT must be outside the read-only RESULT_ROOT" >&2
  exit 2
fi
log_dir="${OUTPUT_ROOT}/logs"
mkdir -p "${log_dir}"
export_args="ALL,PROJECT_ROOT=${PROJECT_ROOT},RESULT_ROOT=${RESULT_ROOT},OUTPUT_ROOT=${OUTPUT_ROOT},FIXED_O2_SOURCE_ROOT=${FIXED_O2_SOURCE_ROOT},R_MODULE=${R_MODULE},CLASS_THRESHOLD=${CLASS_THRESHOLD},CLASS_LOWER_BOUND=${CLASS_LOWER_BOUND},CLASS_UPPER_BOUND=${CLASS_UPPER_BOUND},CLASS_BOUNDARY_RULE=${CLASS_BOUNDARY_RULE},N_BOOT=${N_BOOT},PERMUTATIONS=${PERMUTATIONS},TRAJECTORY_STEP=${TRAJECTORY_STEP},MAX_PAIRS=${MAX_PAIRS},MAX_SEEDS=${MAX_SEEDS}"
command=(
  sbatch --job-name=o2_joint_coupling --cpus-per-task="${CPUS}" --mem="${MEMORY}" --time="${TIME_LIMIT}"
  --output="${log_dir}/%x_%j.out" --error="${log_dir}/%x_%j.err" --export="${export_args}"
)
[[ -n "${QOS}" ]] && command+=(--qos="${QOS}")
[[ -n "${PARTITION}" ]] && command+=(--partition="${PARTITION}")
[[ -n "${ACCOUNT}" ]] && command+=(--account="${ACCOUNT}")
command+=("${worker}")

if [[ "${DRY_RUN,,}" == "true" ]]; then
  printf 'DRY_RUN:'
  printf ' %q' "${command[@]}"
  printf '\n'
else
  "${command[@]}"
fi
