#!/bin/bash
# Submit fixed-O2 workflow jobs through one Slurm entry point.
#
# WORKFLOW_PART=simulation runs the fixed-O2 simulation array.
# WORKFLOW_PART=fixo2_analysis runs FixO2_invivo.R with the requested analysis
# parts after the simulation array dependency, when both are requested.

set -euo pipefail

O2SD_SHELL_UTILS="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../util" && pwd)/o2_supply_demand_map_shell_utils.sh"
# shellcheck source=../../util/o2_supply_demand_map_shell_utils.sh
source "${O2SD_SHELL_UTILS}"

ORIGINAL_SUBMIT_ARGS=("$@")

usage() {
  cat <<'EOF'
Usage:
  bash submit_fix_o2_simulation_array.sh \
    --run_parts=simulation \
    --fit_dir=oxygen/results/fit_invivo_O2_buffering_500seed \
    --simulation=invivo \
    --seeds=1:500 \
    --o2_values=0,0.1,0.5,1,2,5 \
    --initial_ploidy_values=2,4 \
    --time_days=1000 \
    --n_sim=3

  bash submit_fix_o2_simulation_array.sh \
    --run_parts=analysis \
    --run_dir=oxygen/results/fit_invivo_O2_buffering_500seed \
    --simulation_dir=oxygen/results/O2_fixed_simulation \
    --analysis_out_dir=oxygen/results/analysis/best_fit_parameter_feature/01_fixed_o2/FixO2_invivo_500seed \
    --analysis_cpus_per_task=62 \
    --analysis_mem=128G \
    --analysis_time_limit=12:00:00

Process selection:
  --run_parts=PARTS                simulation, analytical_simulation_agreement,
                                  attractors, counterfactual_trajectories,
                                  validation, analysis, or all.
                                  Here all means simulation array plus
                                  FixO2_invivo.R --run_parts=all; analysis
                                  means FixO2_invivo.R --run_parts=all only.

Required for simulation:
  --fit_dir=DIR or --run_dir=DIR   Parent result dir containing seedXX/best_params.tsv.
  --simulation=MODE                invivo, invitro, or joint.
  --o2_values=CSV or --o2=CSV      Fixed O2 values.
  --initial_ploidy_values=CSV      Initial ploidy values. --initial_ploidy=CSV also works.
  --time_days=N                    Simulation horizon in days.
  --n_sim=N                        Number of replicate simulations per seed/O2/ploidy.

Simulation options:
  --seeds=CSV/RANGE                Seed filter, e.g. 1:500, 25,322, or seed25,seed322.
  --seed_ids=CSV                   Alias that also sets agreement seed_ids.
  --out_dir=DIR                    Simulation output root; defaults PROJECT_ROOT/oxygen/results/O2_fixed_simulation.
  --simulation_dir=DIR             Alias for simulation output root and agreement simulation input dir.
  --dt=N
  --save_every_days=N
  --report_dt=N
  --initial_cells=N
  --joint_scope=shared_invivo|invitro_effective
  --Crowding=TRUE|FALSE
  --O2_growth=TRUE|FALSE
  --start_with=MODE
  --ploidy_O2_death=MODE
  --skip_existing=TRUE|FALSE       Worker skips tasks with all required outputs; defaults FALSE.

Agreement options:
  --agreement_out_dir=DIR          Analysis output root; defaults ~/oxygen/results/analysis/best_fit_parameter_feature/01_fixed_o2/FixO2_invivo_500seed.
  --agreement_analysis_dir=DIR     Source analysis dir for mode labels; defaults agreement_out_dir.
  --agreement_simulation_dir=DIR   Fixed-O2 simulation result dir; defaults simulation_dir.
  --agreement_fit_dir=DIR          Fitting result dir for final objectives; defaults run_dir/fit_dir.
  --agreement_run_dir=DIR          Fitting result dir used for analytical solutions; defaults agreement_fit_dir.
  --agreement_seed_ids=CSV         Optional seed filter. --seed_ids also works.
  --agreement_time_points=CSV      Defaults 25,50,100,200,300,500,700,1000.
  --agreement_o2_values=CSV        Defaults 0,0.1,0.5,1,2,5.
  --agreement_initial_ploidy_values=CSV
                                  Defaults 2,4.
  --agreement_simulation_ids=CSV   Defaults 1,2,3.
  --agreement_analytical_methods=CSV
                                  Defaults eigen,expm.
  --agreement_objective_transform=identity|log10
                                  Defaults identity.
  --agreement_cache_all_times=TRUE|FALSE
                                  Defaults TRUE.
  --agreement_recompute=TRUE|FALSE
                                  Defaults FALSE.
  --agreement_recompute_analytical=TRUE|FALSE
                                  Defaults FALSE.
  --agreement_progress_every=N     Defaults 200.
  --agreement_analytical_cache_table=FILE
  --agreement_all_time_simulation_metric_table=FILE
  --agreement_simulation_metric_table=FILE
  --agreement_simulation_summary_table=FILE
  --agreement_data_table=FILE

Analysis options:
  --analysis_out_dir=DIR           Alias for agreement_out_dir.
  --analysis_cpus_per_task=N       Alias for agreement_cpus_per_task.
  --analysis_n_workers=N           Alias for agreement_n_workers.
  --analysis_mem=SIZE              Alias for agreement_mem.
  --analysis_qos=NAME              Alias for agreement_qos.
  --analysis_time_limit=HH:MM:SS   Alias for agreement_time_limit.
  --analysis_job_name=NAME         Alias for agreement_job_name.
  --mode_reference_o2=VALUE        Forwarded to FixO2_invivo.R.
  --attractor_o2_grid=CSV          Forwarded to FixO2_invivo.R.
  --o2_grid=CSV                    Forwarded to FixO2_invivo.R for
                                  non-attractor fixed-O2 analyses.
  --simulation_ids=CSV             Forwarded to FixO2_invivo.R and agreement.
  --simulation_n_core=N            Cores used by FixO2_invivo.R if it needs to
                                  generate missing simulation files.
  --simulation_worker_threads=N    Threads per missing-simulation worker.
  --time_grid=CSV                  Forwarded to FixO2_invivo.R.
  --dt_grid=CSV                    Forwarded to FixO2_invivo.R.
  --plot_dt=N                      Forwarded to FixO2_invivo.R.
  --generate_missing_simulation=TRUE|FALSE
                                  Forwarded to FixO2_invivo.R.
  --include_simulation=TRUE|FALSE  Forwarded to FixO2_invivo.R.
  --generate_figures=TRUE|FALSE    Forwarded to FixO2_invivo.R.
  --max_seeds=N                    Forwarded to FixO2_invivo.R.
  --render_html_report=TRUE|FALSE  Forwarded to FixO2_invivo.R; defaults TRUE there.
  --html_report_dir=DIR            Forwarded to FixO2_invivo.R; defaults
                                  analysis_out_dir/html_report.
  --html_report_basename=NAME      Forwarded to FixO2_invivo.R; defaults index.
  --html_report_script=FILE        Optional override for render_fixo2_invivo_report.R.

HPC options:
  --project_root=DIR
  --simulation_script=FILE
  --fixo2_script=FILE
  --array_script=FILE              Defaults run_fix_o2_simulation_array.sub.
  --tasks_tsv=FILE                 Existing or generated simulation task table path.
  --refresh_task_list=TRUE|FALSE   Rebuild task table; defaults TRUE.
  --array_spec=SPEC                Slurm array spec; defaults to all task IDs.
  --array_max_concurrent=N         Adds %N to simulation array spec only.
  --cpus_per_task=N                Simulation CPUs per task. Defaults 1.
  --mem=SIZE                       Simulation memory. Defaults 2G.
  --qos=NAME                       Simulation QoS. Defaults small.
  --time_limit=HH:MM:SS            Simulation walltime. Defaults 4:00:00.
  --agreement_cpus_per_task=N      Agreement CPUs. Defaults cpus_per_task or 8.
  --agreement_n_workers=N          Agreement R workers. Defaults agreement_cpus_per_task.
  --agreement_mem=SIZE             Agreement memory. Defaults mem or 32G.
  --agreement_qos=NAME             Agreement QoS. Defaults qos or xxlarge.
  --agreement_time_limit=HH:MM:SS  Agreement walltime. Defaults time_limit or 04:00:00.
  --job_name=NAME                  Simulation job name. Defaults o2fix_sim.
  --agreement_job_name=NAME        Agreement job name. Defaults fixo2_agree.
  --log_root=DIR                   Defaults PROJECT_ROOT/oxygen/results/log.
  --r_module=MODULE                Defaults R/4.4.
  --dependency=SPEC                Extra dependency for agreement job.
  --dry_run=TRUE|FALSE             Print sbatch command(s) without submitting.
EOF
}

resolve_project_path() {
  local path="${1-}"
  if [[ -z "${path}" ]]; then
    return 0
  fi
  case "${path}" in
    "~")
      printf "%s" "${PROJECT_ROOT}"
      ;;
    "~/"*)
      printf "%s/%s" "${PROJECT_ROOT}" "${path:2}"
      ;;
    /*)
      printf "%s" "${path}"
      ;;
    *)
      printf "%s/%s" "${PROJECT_ROOT}" "${path}"
      ;;
  esac
}

write_env_line() {
  local file="$1"
  local name="$2"
  local value="${3-}"
  local quoted
  printf -v quoted "%q" "${value}"
  printf "%s=%s\n" "${name}" "${quoted}" >> "${file}"
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

append_fixo2_part() {
  local part="$1"
  if [[ "${part}" == "all" ]]; then
    FIXO2_RUN_PARTS="all"
    return 0
  fi
  if [[ "${FIXO2_RUN_PARTS}" == "all" ]]; then
    return 0
  fi
  if [[ -z "${FIXO2_RUN_PARTS}" ]]; then
    FIXO2_RUN_PARTS="${part}"
  elif [[ ",${FIXO2_RUN_PARTS}," != *",${part},"* ]]; then
    FIXO2_RUN_PARTS="${FIXO2_RUN_PARTS},${part}"
  fi
}

normalize_parts() {
  local input="${1:-simulation}"
  local item
  RUN_SIMULATION="FALSE"
  RUN_FIXO2_ANALYSIS="FALSE"
  FIXO2_RUN_PARTS=""
  IFS=',' read -ra items <<< "${input}"
  for item in "${items[@]}"; do
    item="$(echo "${item}" | tr '[:upper:]' '[:lower:]' | xargs)"
    case "${item}" in
      all)
        RUN_SIMULATION="TRUE"
        RUN_FIXO2_ANALYSIS="TRUE"
        append_fixo2_part "all"
        ;;
      simulation)
        RUN_SIMULATION="TRUE"
        ;;
      analysis|fixo2|fixo2_analysis|full_analysis)
        RUN_FIXO2_ANALYSIS="TRUE"
        append_fixo2_part "all"
        ;;
      attractor|attractors)
        RUN_FIXO2_ANALYSIS="TRUE"
        append_fixo2_part "attractors"
        ;;
      counterfactual|trajectories|counterfactual_trajectories)
        RUN_FIXO2_ANALYSIS="TRUE"
        append_fixo2_part "counterfactual_trajectories"
        ;;
      validation|simulation_validation|representative_simulation)
        RUN_FIXO2_ANALYSIS="TRUE"
        append_fixo2_part "simulation"
        ;;
      agreement|analytical_agreement|analytical_simulation|analytical_simulation_agreement|scatter|scatters)
        RUN_FIXO2_ANALYSIS="TRUE"
        append_fixo2_part "analytical_simulation_agreement"
        ;;
      "")
        ;;
      *)
        echo "Unknown run_parts value: ${item}" >&2
        exit 2
        ;;
    esac
  done
}

parse_args() {
  for arg in "$@"; do
    case "${arg}" in
      --help|-h) usage; exit 0 ;;
      --run_parts=*|--workflow_parts=*) RUN_PARTS="${arg#*=}" ;;
      --project_root=*) PROJECT_ROOT="${arg#*=}" ;;
      --simulation_script=*) SIM_SCRIPT="${arg#*=}" ;;
      --fixo2_script=*) FIXO2_SCRIPT="${arg#*=}" ;;
      --array_script=*) ARRAY_SCRIPT="${arg#*=}" ;;
      --tasks_tsv=*|--task_list=*) TASKS_TSV="${arg#*=}" ;;
      --refresh_task_list=*) REFRESH_TASK_LIST="${arg#*=}" ;;
      --fit_dir=*) FIT_DIR="${arg#*=}" ;;
      --run_dir=*) RUN_DIR="${arg#*=}" ;;
      --simulation=*) SIMULATION="${arg#*=}" ;;
      --seeds=*) SEEDS="${arg#*=}" ;;
      --seed_ids=*) SEEDS="${arg#*=}"; AGREEMENT_SEED_IDS="${arg#*=}" ;;
      --o2_values=*) O2_VALUES="${arg#*=}"; AGREEMENT_O2_VALUES="${arg#*=}" ;;
      --o2=*) O2_VALUES="${arg#*=}"; AGREEMENT_O2_VALUES="${arg#*=}" ;;
      --initial_ploidy_values=*) INITIAL_PLOIDY_VALUES="${arg#*=}"; AGREEMENT_INITIAL_PLOIDY_VALUES="${arg#*=}" ;;
      --initial_ploidy=*) INITIAL_PLOIDY_VALUES="${arg#*=}"; AGREEMENT_INITIAL_PLOIDY_VALUES="${arg#*=}" ;;
      --time_days=*) TIME_DAYS="${arg#*=}" ;;
      --n_sim=*) N_SIM="${arg#*=}" ;;
      --out_dir=*) OUT_ROOT="${arg#*=}" ;;
      --simulation_dir=*) OUT_ROOT="${arg#*=}"; AGREEMENT_SIMULATION_DIR="${arg#*=}" ;;
      --o2_grid=*) FIXO2_O2_GRID="${arg#*=}" ;;
      --mode_reference_o2=*) MODE_REFERENCE_O2="${arg#*=}" ;;
      --attractor_o2_grid=*) ATTRACTOR_O2_GRID="${arg#*=}" ;;
      --dt=*) DT="${arg#*=}" ;;
      --save_every_days=*) SAVE_EVERY_DAYS="${arg#*=}" ;;
      --report_dt=*) REPORT_DT="${arg#*=}" ;;
      --initial_cells=*) INITIAL_CELLS="${arg#*=}" ;;
      --joint_scope=*) JOINT_SCOPE="${arg#*=}" ;;
      --Crowding=*) CROWDING="${arg#*=}" ;;
      --O2_growth=*) O2_GROWTH="${arg#*=}" ;;
      --start_with=*) START_WITH="${arg#*=}" ;;
      --ploidy_O2_death=*) PLOIDY_O2_DEATH="${arg#*=}" ;;
      --agreement_out_dir=*|--analysis_out_dir=*) AGREEMENT_OUT_DIR="${arg#*=}" ;;
      --agreement_analysis_dir=*|--analysis_dir=*) AGREEMENT_ANALYSIS_DIR="${arg#*=}" ;;
      --agreement_simulation_dir=*) AGREEMENT_SIMULATION_DIR="${arg#*=}" ;;
      --agreement_fit_dir=*) AGREEMENT_FIT_DIR="${arg#*=}" ;;
      --agreement_run_dir=*) AGREEMENT_RUN_DIR="${arg#*=}" ;;
      --agreement_seed_ids=*) AGREEMENT_SEED_IDS="${arg#*=}" ;;
      --agreement_time_points=*|--time_points=*) AGREEMENT_TIME_POINTS="${arg#*=}" ;;
      --agreement_o2_values=*) AGREEMENT_O2_VALUES="${arg#*=}" ;;
      --agreement_initial_ploidy_values=*) AGREEMENT_INITIAL_PLOIDY_VALUES="${arg#*=}" ;;
      --agreement_simulation_ids=*|--simulation_ids=*) AGREEMENT_SIMULATION_IDS="${arg#*=}" ;;
      --agreement_progress_every=*|--progress_every=*) AGREEMENT_PROGRESS_EVERY="${arg#*=}" ;;
      --agreement_analytical_methods=*|--analytical_methods=*) AGREEMENT_ANALYTICAL_METHODS="${arg#*=}" ;;
      --agreement_objective_transform=*|--objective_transform=*) AGREEMENT_OBJECTIVE_TRANSFORM="${arg#*=}" ;;
      --agreement_recompute=*|--recompute=*) AGREEMENT_RECOMPUTE="${arg#*=}" ;;
      --agreement_recompute_analytical=*|--recompute_analytical=*) AGREEMENT_RECOMPUTE_ANALYTICAL="${arg#*=}" ;;
      --agreement_cache_all_times=*|--cache_all_times=*) AGREEMENT_CACHE_ALL_TIMES="${arg#*=}" ;;
      --agreement_analytical_cache_table=*|--analytical_cache_table=*) AGREEMENT_ANALYTICAL_CACHE_TABLE="${arg#*=}" ;;
      --agreement_all_time_simulation_metric_table=*|--all_time_simulation_metric_table=*) AGREEMENT_ALL_TIME_SIMULATION_METRIC_TABLE="${arg#*=}" ;;
      --agreement_simulation_metric_table=*|--simulation_metric_table=*) AGREEMENT_SIMULATION_METRIC_TABLE="${arg#*=}" ;;
      --agreement_simulation_summary_table=*|--simulation_summary_table=*) AGREEMENT_SIMULATION_SUMMARY_TABLE="${arg#*=}" ;;
      --agreement_data_table=*|--scatter_data_table=*) AGREEMENT_DATA_TABLE="${arg#*=}" ;;
      --simulation_n_core=*|--simulation_workers=*) FIXO2_SIMULATION_N_CORE="${arg#*=}" ;;
      --simulation_worker_threads=*|--simulation_threads_per_worker=*) FIXO2_SIMULATION_WORKER_THREADS="${arg#*=}" ;;
      --time_grid=*) FIXO2_TIME_GRID="${arg#*=}" ;;
      --dt_grid=*) FIXO2_DT_GRID="${arg#*=}" ;;
      --plot_dt=*) FIXO2_PLOT_DT="${arg#*=}" ;;
      --generate_missing_simulation=*) FIXO2_GENERATE_MISSING_SIMULATION="${arg#*=}" ;;
      --include_simulation=*) FIXO2_INCLUDE_SIMULATION="${arg#*=}" ;;
      --generate_figures=*) FIXO2_GENERATE_FIGURES="${arg#*=}" ;;
      --max_seeds=*) FIXO2_MAX_SEEDS="${arg#*=}" ;;
      --simulation_seed_selection_n_workers=*|--seed_selection_n_workers=*) FIXO2_SIMULATION_SEED_SELECTION_N_WORKERS="${arg#*=}" ;;
      --render_html_report=*|--generate_html_report=*) FIXO2_RENDER_HTML_REPORT="${arg#*=}" ;;
      --html_report_dir=*|--report_out_dir=*) FIXO2_HTML_REPORT_DIR="${arg#*=}" ;;
      --html_report_basename=*|--report_basename=*) FIXO2_HTML_REPORT_BASENAME="${arg#*=}" ;;
      --html_report_script=*|--report_script=*) FIXO2_HTML_REPORT_SCRIPT="${arg#*=}" ;;
      --array_spec=*|--array=*) ARRAY_SPEC="${arg#*=}" ;;
      --array_max_concurrent=*|--max_concurrent=*) ARRAY_MAX_CONCURRENT="${arg#*=}" ;;
      --cpus_per_task=*|--cpus-per-task=*) CPUS_PER_TASK="${arg#*=}" ;;
      --mem=*) MEM="${arg#*=}" ;;
      --qos=*) QOS="${arg#*=}" ;;
      --time_limit=*|--time=*) TIME_LIMIT="${arg#*=}" ;;
      --agreement_cpus_per_task=*|--agreement-cpus-per-task=*|--analysis_cpus_per_task=*|--analysis-cpus-per-task=*) AGREEMENT_CPUS_PER_TASK="${arg#*=}" ;;
      --agreement_n_workers=*|--agreement-n-workers=*|--analysis_n_workers=*|--analysis-n-workers=*) AGREEMENT_N_WORKERS="${arg#*=}" ;;
      --agreement_mem=*|--analysis_mem=*) AGREEMENT_MEM="${arg#*=}" ;;
      --agreement_qos=*|--analysis_qos=*) AGREEMENT_QOS="${arg#*=}" ;;
      --agreement_time_limit=*|--agreement_time=*|--analysis_time_limit=*|--analysis_time=*) AGREEMENT_TIME_LIMIT="${arg#*=}" ;;
      --job_name=*) JOB_NAME="${arg#*=}" ;;
      --agreement_job_name=*|--analysis_job_name=*) AGREEMENT_JOB_NAME="${arg#*=}" ;;
      --log_root=*|--log_dir=*) LOG_ROOT="${arg#*=}" ;;
      --r_module=*) R_MODULE="${arg#*=}" ;;
      --skip_existing=*) SKIP_EXISTING="${arg#*=}" ;;
      --dependency=*) AGREEMENT_DEPENDENCY="${arg#*=}" ;;
      --dry_run=*) DRY_RUN="${arg#*=}" ;;
      *) echo "Unknown argument: ${arg}" >&2; usage >&2; exit 2 ;;
    esac
  done
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKFLOW_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
DEFAULT_PROJECT_ROOT="$(cd "${WORKFLOW_ROOT}/../../.." && pwd)"

PROJECT_ROOT="${PROJECT_ROOT:-${DEFAULT_PROJECT_ROOT}}"
RUN_PARTS="${RUN_PARTS:-simulation}"
SIM_SCRIPT="${SIM_SCRIPT:-}"
FIXO2_SCRIPT="${FIXO2_SCRIPT:-}"
ARRAY_SCRIPT="${ARRAY_SCRIPT:-}"
TASKS_TSV="${TASKS_TSV:-}"
REFRESH_TASK_LIST="${REFRESH_TASK_LIST:-TRUE}"
FIT_DIR="${FIT_DIR:-}"
RUN_DIR="${RUN_DIR:-}"
SIMULATION="${SIMULATION:-}"
SEEDS="${SEEDS:-}"
O2_VALUES="${O2_VALUES:-}"
INITIAL_PLOIDY_VALUES="${INITIAL_PLOIDY_VALUES:-}"
TIME_DAYS="${TIME_DAYS:-}"
N_SIM="${N_SIM:-}"
OUT_ROOT="${OUT_ROOT:-}"
FIXO2_O2_GRID="${FIXO2_O2_GRID:-}"
MODE_REFERENCE_O2="${MODE_REFERENCE_O2:-}"
ATTRACTOR_O2_GRID="${ATTRACTOR_O2_GRID:-}"
DT="${DT:-}"
SAVE_EVERY_DAYS="${SAVE_EVERY_DAYS:-}"
REPORT_DT="${REPORT_DT:-}"
INITIAL_CELLS="${INITIAL_CELLS:-}"
JOINT_SCOPE="${JOINT_SCOPE:-}"
CROWDING="${CROWDING:-}"
O2_GROWTH="${O2_GROWTH:-}"
START_WITH="${START_WITH:-}"
PLOIDY_O2_DEATH="${PLOIDY_O2_DEATH:-}"
AGREEMENT_OUT_DIR="${AGREEMENT_OUT_DIR:-}"
AGREEMENT_ANALYSIS_DIR="${AGREEMENT_ANALYSIS_DIR:-}"
AGREEMENT_SIMULATION_DIR="${AGREEMENT_SIMULATION_DIR:-}"
AGREEMENT_FIT_DIR="${AGREEMENT_FIT_DIR:-}"
AGREEMENT_RUN_DIR="${AGREEMENT_RUN_DIR:-}"
AGREEMENT_SEED_IDS="${AGREEMENT_SEED_IDS:-}"
AGREEMENT_TIME_POINTS="${AGREEMENT_TIME_POINTS:-25,50,100,200,300,500,700,1000}"
AGREEMENT_O2_VALUES="${AGREEMENT_O2_VALUES:-0,0.1,0.5,1,2,5}"
AGREEMENT_INITIAL_PLOIDY_VALUES="${AGREEMENT_INITIAL_PLOIDY_VALUES:-2,4}"
AGREEMENT_SIMULATION_IDS="${AGREEMENT_SIMULATION_IDS:-1,2,3}"
AGREEMENT_PROGRESS_EVERY="${AGREEMENT_PROGRESS_EVERY:-200}"
AGREEMENT_ANALYTICAL_METHODS="${AGREEMENT_ANALYTICAL_METHODS:-eigen,expm}"
AGREEMENT_OBJECTIVE_TRANSFORM="${AGREEMENT_OBJECTIVE_TRANSFORM:-identity}"
AGREEMENT_RECOMPUTE="${AGREEMENT_RECOMPUTE:-FALSE}"
AGREEMENT_RECOMPUTE_ANALYTICAL="${AGREEMENT_RECOMPUTE_ANALYTICAL:-FALSE}"
AGREEMENT_CACHE_ALL_TIMES="${AGREEMENT_CACHE_ALL_TIMES:-TRUE}"
AGREEMENT_ANALYTICAL_CACHE_TABLE="${AGREEMENT_ANALYTICAL_CACHE_TABLE:-}"
AGREEMENT_ALL_TIME_SIMULATION_METRIC_TABLE="${AGREEMENT_ALL_TIME_SIMULATION_METRIC_TABLE:-}"
AGREEMENT_SIMULATION_METRIC_TABLE="${AGREEMENT_SIMULATION_METRIC_TABLE:-}"
AGREEMENT_SIMULATION_SUMMARY_TABLE="${AGREEMENT_SIMULATION_SUMMARY_TABLE:-}"
AGREEMENT_DATA_TABLE="${AGREEMENT_DATA_TABLE:-}"
FIXO2_SIMULATION_N_CORE="${FIXO2_SIMULATION_N_CORE:-}"
FIXO2_SIMULATION_WORKER_THREADS="${FIXO2_SIMULATION_WORKER_THREADS:-}"
FIXO2_TIME_GRID="${FIXO2_TIME_GRID:-}"
FIXO2_DT_GRID="${FIXO2_DT_GRID:-}"
FIXO2_PLOT_DT="${FIXO2_PLOT_DT:-}"
FIXO2_GENERATE_MISSING_SIMULATION="${FIXO2_GENERATE_MISSING_SIMULATION:-}"
FIXO2_INCLUDE_SIMULATION="${FIXO2_INCLUDE_SIMULATION:-}"
FIXO2_GENERATE_FIGURES="${FIXO2_GENERATE_FIGURES:-}"
FIXO2_MAX_SEEDS="${FIXO2_MAX_SEEDS:-}"
FIXO2_SIMULATION_SEED_SELECTION_N_WORKERS="${FIXO2_SIMULATION_SEED_SELECTION_N_WORKERS:-}"
FIXO2_RENDER_HTML_REPORT="${FIXO2_RENDER_HTML_REPORT:-}"
FIXO2_HTML_REPORT_DIR="${FIXO2_HTML_REPORT_DIR:-}"
FIXO2_HTML_REPORT_BASENAME="${FIXO2_HTML_REPORT_BASENAME:-}"
FIXO2_HTML_REPORT_SCRIPT="${FIXO2_HTML_REPORT_SCRIPT:-}"
ARRAY_SPEC="${ARRAY_SPEC:-}"
ARRAY_MAX_CONCURRENT="${ARRAY_MAX_CONCURRENT:-}"
CPUS_PER_TASK="${CPUS_PER_TASK:-}"
MEM="${MEM:-}"
QOS="${QOS:-}"
TIME_LIMIT="${TIME_LIMIT:-}"
AGREEMENT_CPUS_PER_TASK="${AGREEMENT_CPUS_PER_TASK:-}"
AGREEMENT_N_WORKERS="${AGREEMENT_N_WORKERS:-}"
AGREEMENT_MEM="${AGREEMENT_MEM:-}"
AGREEMENT_QOS="${AGREEMENT_QOS:-}"
AGREEMENT_TIME_LIMIT="${AGREEMENT_TIME_LIMIT:-}"
JOB_NAME="${JOB_NAME:-o2fix_sim}"
AGREEMENT_JOB_NAME="${AGREEMENT_JOB_NAME:-}"
LOG_ROOT="${LOG_ROOT:-}"
R_MODULE="${R_MODULE:-R/4.4}"
SKIP_EXISTING="${SKIP_EXISTING:-FALSE}"
AGREEMENT_DEPENDENCY="${AGREEMENT_DEPENDENCY:-}"
DRY_RUN="${DRY_RUN:-FALSE}"

parse_args "$@"
normalize_parts "${RUN_PARTS}"
if [[ -z "${AGREEMENT_JOB_NAME}" ]]; then
  if [[ "${FIXO2_RUN_PARTS}" == "analytical_simulation_agreement" ]]; then
    AGREEMENT_JOB_NAME="fixo2_agree"
  else
    AGREEMENT_JOB_NAME="fixo2_analysis"
  fi
fi

PROJECT_ROOT="$(cd "${PROJECT_ROOT}" && pwd)"
SIM_SCRIPT="${SIM_SCRIPT:-${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/simulation/o2/fixed_o2/run_fixed_o2_simulation.R}"
FIXO2_SCRIPT="${FIXO2_SCRIPT:-${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/analysis/best_fit_parameter_feature/01_fixed_o2/FixO2_invivo.R}"
ARRAY_SCRIPT="${ARRAY_SCRIPT:-${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/hpc/array_workers/run_fix_o2_simulation_array.sub}"
SIM_SCRIPT="$(resolve_project_path "${SIM_SCRIPT}")"
FIXO2_SCRIPT="$(resolve_project_path "${FIXO2_SCRIPT}")"
ARRAY_SCRIPT="$(resolve_project_path "${ARRAY_SCRIPT}")"
if [[ -n "${FIT_DIR}" ]]; then FIT_DIR="$(resolve_project_path "${FIT_DIR}")"; fi
if [[ -n "${RUN_DIR}" ]]; then RUN_DIR="$(resolve_project_path "${RUN_DIR}")"; fi
LOG_ROOT="${LOG_ROOT:-${PROJECT_ROOT}/oxygen/results/log}"
LOG_ROOT="$(resolve_project_path "${LOG_ROOT}")"
mkdir -p "${LOG_ROOT}"
LOG_ROOT="$(cd "${LOG_ROOT}" && pwd)"

SIM_SCRIPT="$(cd "$(dirname "${SIM_SCRIPT}")" && pwd)/$(basename "${SIM_SCRIPT}")"
FIXO2_SCRIPT="$(cd "$(dirname "${FIXO2_SCRIPT}")" && pwd)/$(basename "${FIXO2_SCRIPT}")"
ARRAY_SCRIPT="$(cd "$(dirname "${ARRAY_SCRIPT}")" && pwd)/$(basename "${ARRAY_SCRIPT}")"

if [[ ! -f "${ARRAY_SCRIPT}" ]]; then
  echo "Missing Slurm worker script: ${ARRAY_SCRIPT}" >&2
  exit 1
fi
if truthy "${RUN_SIMULATION}" && [[ ! -f "${SIM_SCRIPT}" ]]; then
  echo "Missing simulation script: ${SIM_SCRIPT}" >&2
  exit 1
fi
if truthy "${RUN_FIXO2_ANALYSIS}" && [[ ! -f "${FIXO2_SCRIPT}" ]]; then
  echo "Missing FixO2 workflow script: ${FIXO2_SCRIPT}" >&2
  exit 1
fi

if truthy "${RUN_SIMULATION}"; then
  if [[ -z "${RUN_DIR}" && -z "${FIT_DIR}" ]]; then
    echo "Simulation requires --fit_dir or --run_dir pointing to a parent directory with seedXX subdirectories." >&2
    exit 2
  fi
  for required_name in SIMULATION O2_VALUES INITIAL_PLOIDY_VALUES TIME_DAYS N_SIM; do
    if [[ -z "${!required_name}" ]]; then
      echo "Missing required simulation option: ${required_name}" >&2
      usage >&2
      exit 2
    fi
  done
fi

simulation_cpus="${CPUS_PER_TASK:-1}"
simulation_mem="${MEM:-2G}"
simulation_qos="${QOS:-small}"
simulation_time="${TIME_LIMIT:-4:00:00}"
agreement_cpus="${AGREEMENT_CPUS_PER_TASK:-${CPUS_PER_TASK:-8}}"
agreement_workers="${AGREEMENT_N_WORKERS:-${agreement_cpus}}"
agreement_mem="${AGREEMENT_MEM:-${MEM:-32G}}"
agreement_qos="${AGREEMENT_QOS:-${QOS:-xxlarge}}"
agreement_time="${AGREEMENT_TIME_LIMIT:-${TIME_LIMIT:-04:00:00}}"

for cpu_pair in "simulation:${simulation_cpus}" "agreement:${agreement_cpus}" "agreement_workers:${agreement_workers}"; do
  name="${cpu_pair%%:*}"
  value="${cpu_pair#*:}"
  if ! [[ "${value}" =~ ^[0-9]+$ ]] || (( value < 1 )); then
    echo "${name} CPU value must be a positive integer, got: ${value}" >&2
    exit 2
  fi
done

timestamp="$(date +%Y%m%d_%H%M%S)"
PROGRESS_LOG="${LOG_ROOT}/fixo2_submit_${timestamp}.log"
: > "${PROGRESS_LOG}"
JOBS_TSV="${LOG_ROOT}/fixo2_submit_jobs_${timestamp}.tsv"
printf "job_type\tjob_id\tarray_spec\tqos\twalltime\tmem\tcpus_per_task\ttask_table\tresult_dir\n" > "${JOBS_TSV}"

if ! truthy "${DRY_RUN}" && ! command -v sbatch >/dev/null 2>&1; then
  echo "sbatch not found; run on an HPC login node or use --dry_run=TRUE." >&2
  exit 1
fi

array_job_id=""
TASKS_TSV_FOR_REPORT=""
STATUS_DIR=""

if truthy "${RUN_SIMULATION}"; then
  if [[ -z "${OUT_ROOT}" ]]; then
    OUT_ROOT="${PROJECT_ROOT}/oxygen/results/O2_fixed_simulation"
  fi
  OUT_ROOT="$(resolve_project_path "${OUT_ROOT}")"
  mkdir -p "${OUT_ROOT}"
  OUT_ROOT="$(cd "${OUT_ROOT}" && pwd)"

  if [[ -z "${TASKS_TSV}" ]]; then
    TASKS_TSV="${OUT_ROOT}/${SIMULATION}/task_list.tsv"
  fi
  TASKS_TSV="$(resolve_project_path "${TASKS_TSV}")"
  mkdir -p "$(dirname "${TASKS_TSV}")"
  TASKS_TSV="$(cd "$(dirname "${TASKS_TSV}")" && pwd)/$(basename "${TASKS_TSV}")"
  TASKS_TSV_FOR_REPORT="${TASKS_TSV}"

  TASK_ROOT="$(cd "$(dirname "${TASKS_TSV}")" && pwd)"
  STATUS_DIR="${TASK_ROOT}/hpc_task_status"
  TASK_ENV_FILE="${TASK_ROOT}/hpc_submit_config.sh"
  JOBS_TSV="${TASK_ROOT}/hpc_submit_jobs.tsv"
  printf "job_type\tjob_id\tarray_spec\tqos\twalltime\tmem\tcpus_per_task\ttask_table\tresult_dir\n" > "${JOBS_TSV}"

  build_cmd=(
    Rscript "${SIM_SCRIPT}"
    "--build_task_list_only=TRUE"
    "--n_core=1"
    "--simulation=${SIMULATION}"
    "--o2_values=${O2_VALUES}"
    "--initial_ploidy_values=${INITIAL_PLOIDY_VALUES}"
    "--time_days=${TIME_DAYS}"
    "--n_sim=${N_SIM}"
    "--out_dir=${OUT_ROOT}"
  )
  if [[ -n "${FIT_DIR}" ]]; then build_cmd+=("--fit_dir=${FIT_DIR}"); fi
  if [[ -n "${RUN_DIR}" ]]; then build_cmd+=("--run_dir=${RUN_DIR}"); fi
  if [[ -n "${SEEDS}" ]]; then build_cmd+=("--seeds=${SEEDS}"); fi
  if [[ -n "${DT}" ]]; then build_cmd+=("--dt=${DT}"); fi
  if [[ -n "${SAVE_EVERY_DAYS}" ]]; then build_cmd+=("--save_every_days=${SAVE_EVERY_DAYS}"); fi
  if [[ -n "${REPORT_DT}" ]]; then build_cmd+=("--report_dt=${REPORT_DT}"); fi
  if [[ -n "${INITIAL_CELLS}" ]]; then build_cmd+=("--initial_cells=${INITIAL_CELLS}"); fi
  if [[ -n "${JOINT_SCOPE}" ]]; then build_cmd+=("--joint_scope=${JOINT_SCOPE}"); fi
  if [[ -n "${CROWDING}" ]]; then build_cmd+=("--Crowding=${CROWDING}"); fi
  if [[ -n "${O2_GROWTH}" ]]; then build_cmd+=("--O2_growth=${O2_GROWTH}"); fi
  if [[ -n "${START_WITH}" ]]; then build_cmd+=("--start_with=${START_WITH}"); fi
  if [[ -n "${PLOIDY_O2_DEATH}" ]]; then build_cmd+=("--ploidy_O2_death=${PLOIDY_O2_DEATH}"); fi

  if [[ ! -f "${TASKS_TSV}" ]] || truthy "${REFRESH_TASK_LIST}"; then
    {
      printf "Build fixed-O2 task list:"
      printf " %q" "${build_cmd[@]}"
      printf "\n"
    } | tee -a "${PROGRESS_LOG}" >&2
    load_r_module
    if ! command -v Rscript >/dev/null 2>&1; then
      echo "Rscript not found after loading ${R_MODULE}." >&2
      exit 1
    fi
    (cd "${PROJECT_ROOT}" && "${build_cmd[@]}") | tee -a "${PROGRESS_LOG}"
  fi

  if [[ ! -f "${TASKS_TSV}" ]]; then
    echo "Missing task table after build: ${TASKS_TSV}" >&2
    exit 1
  fi
  TASK_COUNT=$(( $(wc -l < "${TASKS_TSV}") - 1 ))
  if (( TASK_COUNT < 1 )); then
    echo "Task table has no data rows: ${TASKS_TSV}" >&2
    exit 1
  fi

  if [[ -z "${ARRAY_SPEC}" ]]; then
    ARRAY_SPEC="1-${TASK_COUNT}"
  fi
  if [[ -n "${ARRAY_MAX_CONCURRENT}" && "${ARRAY_SPEC}" != *%* ]]; then
    ARRAY_SPEC="${ARRAY_SPEC}%${ARRAY_MAX_CONCURRENT}"
  fi

  : > "${TASK_ENV_FILE}"
  write_env_line "${TASK_ENV_FILE}" "PROJECT_ROOT" "${PROJECT_ROOT}"
  write_env_line "${TASK_ENV_FILE}" "SIM_SCRIPT" "${SIM_SCRIPT}"
  write_env_line "${TASK_ENV_FILE}" "OUT_ROOT" "${OUT_ROOT}"
  write_env_line "${TASK_ENV_FILE}" "SIMULATION" "${SIMULATION}"
  write_env_line "${TASK_ENV_FILE}" "TIME_DAYS" "${TIME_DAYS}"
  write_env_line "${TASK_ENV_FILE}" "DT" "${DT}"
  write_env_line "${TASK_ENV_FILE}" "SAVE_EVERY_DAYS" "${SAVE_EVERY_DAYS}"
  write_env_line "${TASK_ENV_FILE}" "REPORT_DT" "${REPORT_DT}"
  write_env_line "${TASK_ENV_FILE}" "INITIAL_CELLS" "${INITIAL_CELLS}"
  write_env_line "${TASK_ENV_FILE}" "JOINT_SCOPE" "${JOINT_SCOPE}"
  write_env_line "${TASK_ENV_FILE}" "CROWDING" "${CROWDING}"
  write_env_line "${TASK_ENV_FILE}" "O2_GROWTH" "${O2_GROWTH}"
  write_env_line "${TASK_ENV_FILE}" "START_WITH" "${START_WITH}"
  write_env_line "${TASK_ENV_FILE}" "PLOIDY_O2_DEATH" "${PLOIDY_O2_DEATH}"

  SUBMIT_COMMAND_TEXT="$(shell_join bash "${BASH_SOURCE[0]}" "${ORIGINAL_SUBMIT_ARGS[@]}")"
  {
    printf "key\tvalue\n"
    printf "project_root\t%s\n" "${PROJECT_ROOT}"
    printf "simulation_script\t%s\n" "${SIM_SCRIPT}"
    printf "fixo2_script\t%s\n" "${FIXO2_SCRIPT}"
    printf "array_script\t%s\n" "${ARRAY_SCRIPT}"
    printf "task_table\t%s\n" "${TASKS_TSV}"
    printf "task_count\t%s\n" "${TASK_COUNT}"
    printf "out_root\t%s\n" "${OUT_ROOT}"
    printf "simulation\t%s\n" "${SIMULATION}"
    printf "o2_values\t%s\n" "${O2_VALUES}"
    printf "initial_ploidy_values\t%s\n" "${INITIAL_PLOIDY_VALUES}"
    printf "time_days\t%s\n" "${TIME_DAYS}"
    printf "n_sim\t%s\n" "${N_SIM}"
    printf "array_spec\t%s\n" "${ARRAY_SPEC}"
    printf "cpus_per_task\t%s\n" "${simulation_cpus}"
    printf "mem\t%s\n" "${simulation_mem}"
    printf "qos\t%s\n" "${simulation_qos}"
    printf "time_limit\t%s\n" "${simulation_time}"
    printf "r_module\t%s\n" "${R_MODULE}"
    printf "skip_existing\t%s\n" "${SKIP_EXISTING}"
    printf "run_parts\t%s\n" "${RUN_PARTS}"
    printf "submit_command\t%s\n" "${SUBMIT_COMMAND_TEXT}"
  } > "${TASK_ROOT}/hpc_submit_manifest.tsv"

  export_arg="ALL,WORKFLOW_PART=simulation,PROJECT_ROOT=${PROJECT_ROOT},TASKS_TSV=${TASKS_TSV},TASK_ENV_FILE=${TASK_ENV_FILE},TASK_LOOKUP_COLUMN=task_id,R_MODULE=${R_MODULE},SKIP_EXISTING=${SKIP_EXISTING},STATUS_DIR=${STATUS_DIR},DRY_RUN=FALSE"

  array_job_id="$(submit_or_print \
    "Submit fixed-O2 simulation array" \
    sbatch \
    "--job-name=${JOB_NAME}" \
    "--array=${ARRAY_SPEC}" \
    "--cpus-per-task=${simulation_cpus}" \
    "--mem=${simulation_mem}" \
    "--qos=${simulation_qos}" \
    "--time=${simulation_time}" \
    "--output=${LOG_ROOT}/o2fix_sim_%A_%a.out" \
    "--error=${LOG_ROOT}/o2fix_sim_%A_%a.err" \
    "--export=${export_arg}" \
    "${ARRAY_SCRIPT}")"

  printf "simulation_array\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "${array_job_id}" "${ARRAY_SPEC}" "${simulation_qos}" "${simulation_time}" "${simulation_mem}" "${simulation_cpus}" "${TASKS_TSV}" "${OUT_ROOT}" >> "${JOBS_TSV}"

  echo "Submitted fixed-O2 simulation array"
  echo "  job_id: ${array_job_id}"
  echo "  task_count: ${TASK_COUNT}"
  echo "  array_spec: ${ARRAY_SPEC}"
  echo "  task_table: ${TASKS_TSV}"
  echo "  status_dir: ${STATUS_DIR}"
  echo "  result_dir: ${OUT_ROOT}"
  echo "  logs: ${LOG_ROOT}/o2fix_sim_%A_%a.out"
fi

if truthy "${RUN_FIXO2_ANALYSIS}"; then
  AGREEMENT_OUT_DIR="${AGREEMENT_OUT_DIR:-~/oxygen/results/analysis/best_fit_parameter_feature/01_fixed_o2/FixO2_invivo_500seed}"
  AGREEMENT_ANALYSIS_DIR="${AGREEMENT_ANALYSIS_DIR:-${AGREEMENT_OUT_DIR}}"
  AGREEMENT_SIMULATION_DIR="${AGREEMENT_SIMULATION_DIR:-${OUT_ROOT:-~/oxygen/results/O2_fixed_simulation}}"
  AGREEMENT_FIT_DIR="${AGREEMENT_FIT_DIR:-${RUN_DIR:-${FIT_DIR:-~/oxygen/results/fit_invivo_O2_buffering_500seed}}}"
  AGREEMENT_RUN_DIR="${AGREEMENT_RUN_DIR:-${RUN_DIR:-${AGREEMENT_FIT_DIR}}}"
  FIXO2_O2_GRID="${FIXO2_O2_GRID:-${AGREEMENT_O2_VALUES}}"
  FIXO2_SIMULATION_N_CORE="${FIXO2_SIMULATION_N_CORE:-${agreement_workers}}"
  FIXO2_SIMULATION_WORKER_THREADS="${FIXO2_SIMULATION_WORKER_THREADS:-1}"
  if [[ -z "${FIXO2_RUN_PARTS}" ]]; then
    FIXO2_RUN_PARTS="all"
  fi

  if [[ "${FIXO2_RUN_PARTS}" == "analytical_simulation_agreement" ]]; then
    analysis_workflow_part="analytical_simulation_agreement"
    analysis_log_stem="fixo2_analytical_simulation_agreement"
    analysis_submit_label="Submit fixed-O2 analytical-simulation agreement"
    analysis_result_dir="${AGREEMENT_OUT_DIR}/simulation/analytical_simulation_agreement"
  else
    analysis_workflow_part="fixo2_analysis"
    analysis_log_stem="fixo2_analysis"
    analysis_submit_label="Submit fixed-O2 analysis"
    analysis_result_dir="${AGREEMENT_OUT_DIR}"
  fi

  ARG_ENV_FILE="${LOG_ROOT}/${analysis_log_stem}_${timestamp}.env"
  : > "${ARG_ENV_FILE}"
  write_env_line "${ARG_ENV_FILE}" "PROJECT_ROOT" "${PROJECT_ROOT}"
  write_env_line "${ARG_ENV_FILE}" "WORKFLOW_PART" "${analysis_workflow_part}"
  write_env_line "${ARG_ENV_FILE}" "FIXO2_RUN_PARTS" "${FIXO2_RUN_PARTS}"
  write_env_line "${ARG_ENV_FILE}" "FIXO2_SCRIPT" "${FIXO2_SCRIPT}"
  write_env_line "${ARG_ENV_FILE}" "R_MODULE" "${R_MODULE}"
  write_env_line "${ARG_ENV_FILE}" "AGREEMENT_SIMULATION_DIR" "${AGREEMENT_SIMULATION_DIR}"
  write_env_line "${ARG_ENV_FILE}" "AGREEMENT_ANALYSIS_DIR" "${AGREEMENT_ANALYSIS_DIR}"
  write_env_line "${ARG_ENV_FILE}" "AGREEMENT_FIT_DIR" "${AGREEMENT_FIT_DIR}"
  write_env_line "${ARG_ENV_FILE}" "AGREEMENT_RUN_DIR" "${AGREEMENT_RUN_DIR}"
  write_env_line "${ARG_ENV_FILE}" "AGREEMENT_OUT_DIR" "${AGREEMENT_OUT_DIR}"
  write_env_line "${ARG_ENV_FILE}" "SIMULATION_MODE" "${SIMULATION:-invivo}"
  write_env_line "${ARG_ENV_FILE}" "FIXO2_O2_GRID" "${FIXO2_O2_GRID}"
  write_env_line "${ARG_ENV_FILE}" "MODE_REFERENCE_O2" "${MODE_REFERENCE_O2}"
  write_env_line "${ARG_ENV_FILE}" "ATTRACTOR_O2_GRID" "${ATTRACTOR_O2_GRID}"
  write_env_line "${ARG_ENV_FILE}" "AGREEMENT_SIMULATION_IDS" "${AGREEMENT_SIMULATION_IDS}"
  write_env_line "${ARG_ENV_FILE}" "AGREEMENT_SEED_IDS" "${AGREEMENT_SEED_IDS}"
  write_env_line "${ARG_ENV_FILE}" "AGREEMENT_TIME_POINTS" "${AGREEMENT_TIME_POINTS}"
  write_env_line "${ARG_ENV_FILE}" "AGREEMENT_O2_VALUES" "${AGREEMENT_O2_VALUES}"
  write_env_line "${ARG_ENV_FILE}" "AGREEMENT_INITIAL_PLOIDY_VALUES" "${AGREEMENT_INITIAL_PLOIDY_VALUES}"
  write_env_line "${ARG_ENV_FILE}" "AGREEMENT_PROGRESS_EVERY" "${AGREEMENT_PROGRESS_EVERY}"
  write_env_line "${ARG_ENV_FILE}" "AGREEMENT_ANALYTICAL_METHODS" "${AGREEMENT_ANALYTICAL_METHODS}"
  write_env_line "${ARG_ENV_FILE}" "AGREEMENT_OBJECTIVE_TRANSFORM" "${AGREEMENT_OBJECTIVE_TRANSFORM}"
  write_env_line "${ARG_ENV_FILE}" "AGREEMENT_RECOMPUTE" "${AGREEMENT_RECOMPUTE}"
  write_env_line "${ARG_ENV_FILE}" "AGREEMENT_RECOMPUTE_ANALYTICAL" "${AGREEMENT_RECOMPUTE_ANALYTICAL}"
  write_env_line "${ARG_ENV_FILE}" "AGREEMENT_CACHE_ALL_TIMES" "${AGREEMENT_CACHE_ALL_TIMES}"
  write_env_line "${ARG_ENV_FILE}" "AGREEMENT_N_WORKERS" "${agreement_workers}"
  write_env_line "${ARG_ENV_FILE}" "AGREEMENT_ANALYTICAL_CACHE_TABLE" "${AGREEMENT_ANALYTICAL_CACHE_TABLE}"
  write_env_line "${ARG_ENV_FILE}" "AGREEMENT_ALL_TIME_SIMULATION_METRIC_TABLE" "${AGREEMENT_ALL_TIME_SIMULATION_METRIC_TABLE}"
  write_env_line "${ARG_ENV_FILE}" "AGREEMENT_SIMULATION_METRIC_TABLE" "${AGREEMENT_SIMULATION_METRIC_TABLE}"
  write_env_line "${ARG_ENV_FILE}" "AGREEMENT_SIMULATION_SUMMARY_TABLE" "${AGREEMENT_SIMULATION_SUMMARY_TABLE}"
  write_env_line "${ARG_ENV_FILE}" "AGREEMENT_DATA_TABLE" "${AGREEMENT_DATA_TABLE}"
  write_env_line "${ARG_ENV_FILE}" "FIXO2_SIMULATION_N_CORE" "${FIXO2_SIMULATION_N_CORE}"
  write_env_line "${ARG_ENV_FILE}" "FIXO2_SIMULATION_WORKER_THREADS" "${FIXO2_SIMULATION_WORKER_THREADS}"
  write_env_line "${ARG_ENV_FILE}" "FIXO2_TIME_GRID" "${FIXO2_TIME_GRID}"
  write_env_line "${ARG_ENV_FILE}" "FIXO2_DT_GRID" "${FIXO2_DT_GRID}"
  write_env_line "${ARG_ENV_FILE}" "FIXO2_PLOT_DT" "${FIXO2_PLOT_DT}"
  write_env_line "${ARG_ENV_FILE}" "FIXO2_GENERATE_MISSING_SIMULATION" "${FIXO2_GENERATE_MISSING_SIMULATION}"
  write_env_line "${ARG_ENV_FILE}" "FIXO2_INCLUDE_SIMULATION" "${FIXO2_INCLUDE_SIMULATION}"
  write_env_line "${ARG_ENV_FILE}" "FIXO2_GENERATE_FIGURES" "${FIXO2_GENERATE_FIGURES}"
  write_env_line "${ARG_ENV_FILE}" "FIXO2_MAX_SEEDS" "${FIXO2_MAX_SEEDS}"
  write_env_line "${ARG_ENV_FILE}" "FIXO2_SIMULATION_SEED_SELECTION_N_WORKERS" "${FIXO2_SIMULATION_SEED_SELECTION_N_WORKERS}"
  write_env_line "${ARG_ENV_FILE}" "FIXO2_RENDER_HTML_REPORT" "${FIXO2_RENDER_HTML_REPORT}"
  write_env_line "${ARG_ENV_FILE}" "FIXO2_HTML_REPORT_DIR" "${FIXO2_HTML_REPORT_DIR}"
  write_env_line "${ARG_ENV_FILE}" "FIXO2_HTML_REPORT_BASENAME" "${FIXO2_HTML_REPORT_BASENAME}"
  write_env_line "${ARG_ENV_FILE}" "FIXO2_HTML_REPORT_SCRIPT" "${FIXO2_HTML_REPORT_SCRIPT}"
  write_env_line "${ARG_ENV_FILE}" "TIME_DAYS" "${TIME_DAYS}"
  write_env_line "${ARG_ENV_FILE}" "DT" "${DT}"
  write_env_line "${ARG_ENV_FILE}" "SAVE_EVERY_DAYS" "${SAVE_EVERY_DAYS}"
  write_env_line "${ARG_ENV_FILE}" "REPORT_DT" "${REPORT_DT}"
  write_env_line "${ARG_ENV_FILE}" "INITIAL_CELLS" "${INITIAL_CELLS}"
  write_env_line "${ARG_ENV_FILE}" "JOINT_SCOPE" "${JOINT_SCOPE}"
  write_env_line "${ARG_ENV_FILE}" "CROWDING" "${CROWDING}"
  write_env_line "${ARG_ENV_FILE}" "O2_GROWTH" "${O2_GROWTH}"
  write_env_line "${ARG_ENV_FILE}" "START_WITH" "${START_WITH}"
  write_env_line "${ARG_ENV_FILE}" "PLOIDY_O2_DEATH" "${PLOIDY_O2_DEATH}"

  agreement_dependency="${AGREEMENT_DEPENDENCY}"
  if [[ -n "${array_job_id}" ]]; then
    if [[ -n "${agreement_dependency}" ]]; then
      agreement_dependency="${agreement_dependency},afterok:${array_job_id}"
    else
      agreement_dependency="afterok:${array_job_id}"
    fi
  fi

  agreement_cmd=(
    sbatch
    "--job-name=${AGREEMENT_JOB_NAME}"
    "--cpus-per-task=${agreement_cpus}"
    "--mem=${agreement_mem}"
    "--qos=${agreement_qos}"
    "--time=${agreement_time}"
    "--output=${LOG_ROOT}/${analysis_log_stem}_%j.out"
    "--error=${LOG_ROOT}/${analysis_log_stem}_%j.err"
    "--export=ALL,ARG_ENV_FILE=${ARG_ENV_FILE},WORKFLOW_PART=${analysis_workflow_part},DRY_RUN=FALSE"
  )
  if [[ -n "${agreement_dependency}" ]]; then
    agreement_cmd+=("--dependency=${agreement_dependency}")
  fi
  agreement_cmd+=("${ARRAY_SCRIPT}")

  agreement_job_id="$(submit_or_print "${analysis_submit_label}" "${agreement_cmd[@]}")"
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "${analysis_workflow_part}" "${agreement_job_id}" "NA" "${agreement_qos}" "${agreement_time}" "${agreement_mem}" "${agreement_cpus}" "${TASKS_TSV_FOR_REPORT:-NA}" "${analysis_result_dir}" >> "${JOBS_TSV}"

  echo "Submitted fixed-O2 analysis"
  echo "  job_id: ${agreement_job_id}"
  echo "  run_parts: ${FIXO2_RUN_PARTS}"
  echo "  dependency: ${agreement_dependency:-none}"
  echo "  arg_env_file: ${ARG_ENV_FILE}"
  echo "  result_dir: ${analysis_result_dir}"
  echo "  logs: ${LOG_ROOT}/${analysis_log_stem}_%j.out"
fi

echo "Submission manifest: ${JOBS_TSV}"
