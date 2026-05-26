#!/bin/bash
# Submit small joint-only lambda sweeps for soft composite sharing.

set -euo pipefail

DEFAULT_PROJECT_ROOT="/share/lab_crd/lab_crd/taoli/Project/Rescomposite"
PROJECT_ROOT="${PROJECT_ROOT:-${DEFAULT_PROJECT_ROOT}}"
BASE_CONFIG="${BASE_CONFIG:-${PROJECT_ROOT}/oxygen/config/O2G_supply_demand.yaml}"
CONFIG_DIR="${CONFIG_DIR:-${PROJECT_ROOT}/oxygen/config/sweeps}"
OUT_ROOT="${OUT_ROOT:-${PROJECT_ROOT}/oxygen/results}"
SUB_SCRIPT="${SUB_SCRIPT:-${PROJECT_ROOT}/oxygen/code/O2G_supply_demand_MAP/hpc/submit_fit_seed_array_joint_buffering.sub}"
RUNNER_SCRIPT="${RUNNER_SCRIPT:-${PROJECT_ROOT}/oxygen/code/O2G_supply_demand_MAP/runner/run_fit_joint_model_O2G_supply_demand_MAP.sh}"
TOTAL_SEEDS="${TOTAL_SEEDS:-10}"
ARRAY_TASKS="${ARRAY_TASKS:-${TOTAL_SEEDS}}"
SEEDS_PER_TASK="${SEEDS_PER_TASK:-1}"
N_CORES="${N_CORES:-22}"
MEM="${MEM:-32G}"
QOS="${QOS:-xxlarge}"
TIME_LIMIT="${TIME_LIMIT:-12:00:00}"
AUTO_VIZ="${AUTO_VIZ:-TRUE}"
GLUCOSE="${GLUCOSE:-TRUE}"
R_MODULE="${R_MODULE:-R/4.4}"
DRY_RUN="${DRY_RUN:-FALSE}"

PROJECT_ROOT="$(cd "${PROJECT_ROOT}" && pwd)"
BASE_CONFIG="$(cd "$(dirname "${BASE_CONFIG}")" && pwd)/$(basename "${BASE_CONFIG}")"
SUB_SCRIPT="$(cd "$(dirname "${SUB_SCRIPT}")" && pwd)/$(basename "${SUB_SCRIPT}")"
RUNNER_SCRIPT="$(cd "$(dirname "${RUNNER_SCRIPT}")" && pwd)/$(basename "${RUNNER_SCRIPT}")"
mkdir -p "${CONFIG_DIR}" "${OUT_ROOT}"
CONFIG_DIR="$(cd "${CONFIG_DIR}" && pwd)"
OUT_ROOT="$(cd "${OUT_ROOT}" && pwd)"

if [[ ! -f "${BASE_CONFIG}" ]]; then
  echo "Missing base config: ${BASE_CONFIG}" >&2
  exit 1
fi
if [[ ! -f "${SUB_SCRIPT}" ]]; then
  echo "Missing joint submit script: ${SUB_SCRIPT}" >&2
  exit 1
fi
if [[ ! -f "${RUNNER_SCRIPT}" ]]; then
  echo "Missing joint runner script: ${RUNNER_SCRIPT}" >&2
  exit 1
fi
if (( ARRAY_TASKS * SEEDS_PER_TASK != TOTAL_SEEDS )); then
  echo "ARRAY_TASKS * SEEDS_PER_TASK must equal TOTAL_SEEDS." >&2
  exit 1
fi
if [[ "${DRY_RUN}" != "TRUE" && "${DRY_RUN}" != "true" && "${DRY_RUN}" != "1" ]] &&
   ! command -v sbatch >/dev/null 2>&1; then
  echo "sbatch not found; run this script on the HPC login node." >&2
  exit 1
fi

labels=(A B C D E)
lambda_pmis=(0.5 1 2 5 5)
lambda_death=(0.5 1 1 2 5)
lambda_loss=(1 2 5 5 10)

write_sweep_config() {
  local label="$1"
  local pmis="$2"
  local death="$3"
  local loss="$4"
  local run_prefix="$5"
  local config_path="$6"
  Rscript -e '
    args <- commandArgs(TRUE)
    if (!requireNamespace("yaml", quietly = TRUE)) stop("R package yaml is required")
    cfg <- yaml::read_yaml(args[[1]])
    base_dir <- dirname(normalizePath(args[[1]], mustWork = FALSE))
    resolve_path <- function(x) {
      if (is.null(x) || !length(x)) return(x)
      txt <- trimws(as.character(x[[1]]))
      if (!nzchar(txt)) return(x)
      if (grepl("^(/|~)", txt)) return(normalizePath(path.expand(txt), mustWork = FALSE))
      normalizePath(file.path(base_dir, txt), mustWork = FALSE)
    }
    path_keys <- c(
      "data_dir", "seeds_file", "parameter_table", "parameters", "init_params_tsv",
      "invitro_parameter_table", "parameter_table_invitro",
      "fit_objects_dir", "flow_density_path"
    )
    for (key in path_keys) {
      if (!is.null(cfg[[key]])) cfg[[key]] <- resolve_path(cfg[[key]])
    }
    cfg$run_prefix <- args[[2]]
    cfg$out_root <- args[[3]]
    cfg$append_run_prefix_timestamp <- FALSE
    cfg$joint_composite_penalty <- TRUE
    cfg$joint_composite_lambda_pmis <- as.numeric(args[[4]])
    cfg$joint_composite_lambda_death <- as.numeric(args[[5]])
    cfg$joint_composite_lambda_loss <- as.numeric(args[[6]])
    cfg$joint_composite_o2_grid <- list(0, 0.1, 0.2, 0.5, 1, 2)
    cfg$joint_composite_n_grid <- c(44, 66, 88)
    cfg$joint_composite_eps <- 1e-8
    cfg$joint_composite_log_eps <- 1e-8
    cfg$joint_composite_live_loss_mode <- "transition"
    cfg$joint_composite_low_o2_weight <- 2
    cfg$joint_composite_zero_o2_priority_weight <- 4
    cfg$joint_composite_pmis_enabled <- TRUE
    cfg$joint_composite_pmis_lambda <- as.numeric(args[[4]])
    yaml::write_yaml(cfg, args[[7]])
  ' "${BASE_CONFIG}" "${run_prefix}" "${OUT_ROOT}" "${pmis}" "${death}" "${loss}" "${config_path}"
}

echo "Submitting soft-composite joint lambda sweep"
echo "  project_root: ${PROJECT_ROOT}"
echo "  base_config: ${BASE_CONFIG}"
echo "  config_dir: ${CONFIG_DIR}"
echo "  out_root: ${OUT_ROOT}"
echo "  total_seeds: ${TOTAL_SEEDS}"
echo "  resources: qos=${QOS}, time=${TIME_LIMIT}, cpus=${N_CORES}, mem=${MEM}"
echo "  dry_run: ${DRY_RUN}"

for i in "${!labels[@]}"; do
  label="${labels[$i]}"
  pmis="${lambda_pmis[$i]}"
  death="${lambda_death[$i]}"
  loss="${lambda_loss[$i]}"
  run_prefix="fit_joint_softcomp_${label}_${TOTAL_SEEDS}seed"
  config_path="${CONFIG_DIR}/O2G_joint_softcomp_${label}.yaml"
  write_sweep_config "${label}" "${pmis}" "${death}" "${loss}" "${run_prefix}" "${config_path}"

  cmd=(
    sbatch
    --parsable
    --job-name="Ogj_sc${label}"
    --qos="${QOS}"
    --time="${TIME_LIMIT}"
    --cpus-per-task="${N_CORES}"
    --mem="${MEM}"
    --array="1-${ARRAY_TASKS}"
    --output="${OUT_ROOT}/Ogj_sc${label}_%A_%a.out"
    --error="${OUT_ROOT}/Ogj_sc${label}_%A_%a.err"
    --export="ALL,PROJECT_ROOT=${PROJECT_ROOT},RUNNER_SCRIPT=${RUNNER_SCRIPT},CONFIG_PATH=${config_path},OUT_ROOT=${OUT_ROOT},RUN_PREFIX=${run_prefix},TOTAL_SEEDS=${TOTAL_SEEDS},ARRAY_TASKS=${ARRAY_TASKS},SEEDS_PER_TASK=${SEEDS_PER_TASK},N_CORES=${N_CORES},AUTO_VIZ=${AUTO_VIZ},GLUCOSE=${GLUCOSE},R_MODULE=${R_MODULE}"
    "${SUB_SCRIPT}"
  )

  echo "Sweep ${label}: pmis=${pmis}, death=${death}, loss=${loss}, config=${config_path}"
  if [[ "${DRY_RUN}" == "TRUE" || "${DRY_RUN}" == "true" || "${DRY_RUN}" == "1" ]]; then
    printf "  command:"
    printf " %q" "${cmd[@]}"
    printf "\n"
  else
    job_id="$("${cmd[@]}")"
    echo "  submitted: ${job_id}"
  fi
done
