#!/bin/bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKFLOW_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
PROJECT_ROOT="$(cd "${WORKFLOW_ROOT}/../../.." && pwd)"

CONFIG_PATH="${PROJECT_ROOT}/oxygen/config/O2G_supply_demand.yaml"
BASELINE_SEED_DIR="${PROJECT_ROOT}/oxygen/results/fit_invivo_o2_supply_demand_eq21_20260331_011709/seed2"
PROFILE_BOUNDS_TABLE="${BASELINE_SEED_DIR}/parameter_table_input.csv"
OUTPUT_ROOT="${PROJECT_ROOT}/oxygen/results/profile_likelihood_O2G_supply_demand_MAP_$(date +%Y%m%d_%H%M%S)"
MAX_STEPS_PER_DIRECTION="20"
SEEDS_PER_STEP="20"
N_CORES="62"
TARGET_DELTA_OBJECTIVE="0.2"
CI_DELTA_THRESHOLD="1.92"
STEP_FRACTION_INITIAL="0.10"
STEP_FRACTION_MIN="1e-6"
STEP_FRACTION_MAX="0.30"
BOUNDARY_START_TOLERANCE="1e-8"
MAX_ATTEMPTS_PER_STEP="5"
MIN_INTERIOR_POINTS_PER_DIRECTION="4"
CI_REFINE_TOLERANCE="0.05"
MAX_REFINE_STEPS_CI="6"
MAX_REFINE_STEPS_BOUNDARY="4"
BOUNDARY_REFINE_FRACTION_MIN="0.01"
MAX_STEP_GROWTH_FACTOR="1.5"
MAX_STEP_SHRINK_FACTOR="0.67"
USE_SOFT_PRIOR_FOR_PROFILE=""
LAMBDA_PRIOR_FOR_PROFILE=""
ARRAY_CONCURRENCY=""
MEMORY="128G"
WALLTIME="7-00:00:00"
JOB_NAME="o2sd_profile"
MAIL_USER=""
MAIL_TYPE="END"
R_MODULE="R/4.4"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config=*) CONFIG_PATH="${1#*=}" ;;
    --baseline_seed_dir=*) BASELINE_SEED_DIR="${1#*=}" ;;
    --profile_bounds_table=*) PROFILE_BOUNDS_TABLE="${1#*=}" ;;
    --output_root=*) OUTPUT_ROOT="${1#*=}" ;;
    --max_steps_per_direction=*) MAX_STEPS_PER_DIRECTION="${1#*=}" ;;
    --grid_points=*) MAX_STEPS_PER_DIRECTION="${1#*=}" ;;
    --seeds_per_step=*) SEEDS_PER_STEP="${1#*=}" ;;
    --seeds_per_point=*) SEEDS_PER_STEP="${1#*=}" ;;
    --n_cores=*) N_CORES="${1#*=}" ;;
    --target_delta_objective=*) TARGET_DELTA_OBJECTIVE="${1#*=}" ;;
    --ci_delta_threshold=*) CI_DELTA_THRESHOLD="${1#*=}" ;;
    --step_fraction_initial=*) STEP_FRACTION_INITIAL="${1#*=}" ;;
    --step_fraction_min=*) STEP_FRACTION_MIN="${1#*=}" ;;
    --step_fraction_max=*) STEP_FRACTION_MAX="${1#*=}" ;;
    --boundary_start_tolerance=*) BOUNDARY_START_TOLERANCE="${1#*=}" ;;
    --max_attempts_per_step=*) MAX_ATTEMPTS_PER_STEP="${1#*=}" ;;
    --min_interior_points_per_direction=*) MIN_INTERIOR_POINTS_PER_DIRECTION="${1#*=}" ;;
    --ci_refine_tolerance=*) CI_REFINE_TOLERANCE="${1#*=}" ;;
    --max_refine_steps_ci=*) MAX_REFINE_STEPS_CI="${1#*=}" ;;
    --max_refine_steps_boundary=*) MAX_REFINE_STEPS_BOUNDARY="${1#*=}" ;;
    --boundary_refine_fraction_min=*) BOUNDARY_REFINE_FRACTION_MIN="${1#*=}" ;;
    --max_step_growth_factor=*) MAX_STEP_GROWTH_FACTOR="${1#*=}" ;;
    --max_step_shrink_factor=*) MAX_STEP_SHRINK_FACTOR="${1#*=}" ;;
    --use_soft_prior_for_profile=*) USE_SOFT_PRIOR_FOR_PROFILE="${1#*=}" ;;
    --lambda_prior_for_profile=*) LAMBDA_PRIOR_FOR_PROFILE="${1#*=}" ;;
    --array_concurrency=*) ARRAY_CONCURRENCY="${1#*=}" ;;
    --mem=*) MEMORY="${1#*=}" ;;
    --time=*) WALLTIME="${1#*=}" ;;
    --job_name=*) JOB_NAME="${1#*=}" ;;
    --mail_user=*) MAIL_USER="${1#*=}" ;;
    --mail_type=*) MAIL_TYPE="${1#*=}" ;;
    --r_module=*) R_MODULE="${1#*=}" ;;
    --help)
      cat <<'EOF'
Usage:
  bash oxygen/code/O2G_supply_demand_MAP/hpc/submit_profile_likelihood_array.sh \
    --config=/abs/path/to/O2G_supply_demand.yaml \
    --baseline_seed_dir=/abs/path/to/baseline/seed2 \
    --profile_bounds_table=/abs/path/to/parameter_table_input.csv \
    --output_root=/abs/path/to/output_root \
    --max_steps_per_direction=20 \
    --seeds_per_step=20 \
    --n_cores=62

Optional:
  --target_delta_objective=0.2
  --ci_delta_threshold=1.92
  --step_fraction_initial=0.10
  --step_fraction_min=1e-6
  --step_fraction_max=0.30
  --boundary_start_tolerance=1e-8
  --max_attempts_per_step=5
  --min_interior_points_per_direction=4
  --ci_refine_tolerance=0.05
  --max_refine_steps_ci=6
  --max_refine_steps_boundary=4
  --boundary_refine_fraction_min=0.01
  --max_step_growth_factor=1.5
  --max_step_shrink_factor=0.67
  --use_soft_prior_for_profile=TRUE|FALSE
  --lambda_prior_for_profile=0.03
  --array_concurrency=4
  --mem=128G
  --time=7-00:00:00
  --job_name=o2sd_profile
  --mail_user=name@example.org
  --mail_type=END
  --r_module=R/4.4
EOF
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      exit 1
      ;;
  esac
  shift
done

if [[ ! -f "$CONFIG_PATH" ]]; then
  echo "Missing config file: $CONFIG_PATH" >&2
  exit 1
fi
if [[ ! -d "$BASELINE_SEED_DIR" ]]; then
  echo "Missing baseline seed directory: $BASELINE_SEED_DIR" >&2
  exit 1
fi
if [[ ! -f "$PROFILE_BOUNDS_TABLE" ]]; then
  echo "Missing profile bounds table: $PROFILE_BOUNDS_TABLE" >&2
  exit 1
fi

mkdir -p "$OUTPUT_ROOT"
mkdir -p "$OUTPUT_ROOT/slurm_logs"

CONFIG_PATH="$(cd "$(dirname "$CONFIG_PATH")" && pwd)/$(basename "$CONFIG_PATH")"
BASELINE_SEED_DIR="$(cd "$BASELINE_SEED_DIR" && pwd)"
PROFILE_BOUNDS_TABLE="$(cd "$(dirname "$PROFILE_BOUNDS_TABLE")" && pwd)/$(basename "$PROFILE_BOUNDS_TABLE")"
OUTPUT_ROOT="$(cd "$OUTPUT_ROOT" && pwd)"

cp "$PROFILE_BOUNDS_TABLE" "$OUTPUT_ROOT/profile_bounds_table.csv"
cp "$BASELINE_SEED_DIR/best_params.tsv" "$OUTPUT_ROOT/baseline_best_params.tsv"
cp "$BASELINE_SEED_DIR/fit_summary.tsv" "$OUTPUT_ROOT/baseline_fit_summary.tsv"

PARAM_COUNT="$(Rscript -e 'args <- commandArgs(TRUE); tab <- read.csv(args[1], stringsAsFactors = FALSE, check.names = FALSE); if (!("estimate" %in% names(tab)) || !("param_symbol" %in% names(tab))) stop("parameter table must contain estimate and param_symbol columns"); est <- tolower(trimws(as.character(tab$estimate))); keep <- est %in% c("true","t","1","yes","y"); cat(sum(keep), "\n", sep = "")' "$PROFILE_BOUNDS_TABLE")"

if [[ -z "$PARAM_COUNT" || "$PARAM_COUNT" -lt 1 ]]; then
  echo "No estimate=TRUE parameters were found in $PROFILE_BOUNDS_TABLE" >&2
  exit 1
fi

Rscript -e 'args <- commandArgs(TRUE); tab <- read.csv(args[1], stringsAsFactors = FALSE, check.names = FALSE); out <- args[2]; est <- tolower(trimws(as.character(tab$estimate))); keep <- est %in% c("true","t","1","yes","y"); target <- data.frame(param_index = seq_len(sum(keep)), param_symbol = as.character(tab$param_symbol[keep]), stringsAsFactors = FALSE); write.table(target, file = out, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)' "$PROFILE_BOUNDS_TABLE" "$OUTPUT_ROOT/parameter_targets.tsv"

ARRAY_SPEC="1-${PARAM_COUNT}"
if [[ -n "$ARRAY_CONCURRENCY" ]]; then
  ARRAY_SPEC="${ARRAY_SPEC}%${ARRAY_CONCURRENCY}"
fi

if [[ -z "$USE_SOFT_PRIOR_FOR_PROFILE" ]]; then
  USE_SOFT_PRIOR_FOR_PROFILE="$(Rscript -e 'args <- commandArgs(TRUE); tab <- read.delim(args[1], stringsAsFactors = FALSE, check.names = FALSE); hit <- which(tab$metric == "use_soft_prior"); if (length(hit) == 1L) cat(as.character(tab$value[[hit]]), "\n", sep = "") else cat("TRUE\n")' "$BASELINE_SEED_DIR/fit_summary.tsv")"
fi
if [[ -z "$LAMBDA_PRIOR_FOR_PROFILE" ]]; then
  LAMBDA_PRIOR_FOR_PROFILE="$(Rscript -e 'args <- commandArgs(TRUE); tab <- read.delim(args[1], stringsAsFactors = FALSE, check.names = FALSE); hit <- which(tab$metric == "lambda_prior"); if (length(hit) == 1L) cat(as.character(tab$value[[hit]]), "\n", sep = "") else cat("0.03\n")' "$BASELINE_SEED_DIR/fit_summary.tsv")"
fi

MANIFEST_PATH="$OUTPUT_ROOT/submission_manifest.tsv"
{
  printf "key\tvalue\n"
  printf "project_root\t%s\n" "$PROJECT_ROOT"
  printf "config_path\t%s\n" "$CONFIG_PATH"
  printf "baseline_seed_dir\t%s\n" "$BASELINE_SEED_DIR"
  printf "profile_bounds_table\t%s\n" "$PROFILE_BOUNDS_TABLE"
  printf "output_root\t%s\n" "$OUTPUT_ROOT"
  printf "max_steps_per_direction\t%s\n" "$MAX_STEPS_PER_DIRECTION"
  printf "seeds_per_step\t%s\n" "$SEEDS_PER_STEP"
  printf "n_cores\t%s\n" "$N_CORES"
  printf "target_delta_objective\t%s\n" "$TARGET_DELTA_OBJECTIVE"
  printf "ci_delta_threshold\t%s\n" "$CI_DELTA_THRESHOLD"
  printf "step_fraction_initial\t%s\n" "$STEP_FRACTION_INITIAL"
  printf "step_fraction_min\t%s\n" "$STEP_FRACTION_MIN"
  printf "step_fraction_max\t%s\n" "$STEP_FRACTION_MAX"
  printf "boundary_start_tolerance\t%s\n" "$BOUNDARY_START_TOLERANCE"
  printf "max_attempts_per_step\t%s\n" "$MAX_ATTEMPTS_PER_STEP"
  printf "min_interior_points_per_direction\t%s\n" "$MIN_INTERIOR_POINTS_PER_DIRECTION"
  printf "ci_refine_tolerance\t%s\n" "$CI_REFINE_TOLERANCE"
  printf "max_refine_steps_ci\t%s\n" "$MAX_REFINE_STEPS_CI"
  printf "max_refine_steps_boundary\t%s\n" "$MAX_REFINE_STEPS_BOUNDARY"
  printf "boundary_refine_fraction_min\t%s\n" "$BOUNDARY_REFINE_FRACTION_MIN"
  printf "max_step_growth_factor\t%s\n" "$MAX_STEP_GROWTH_FACTOR"
  printf "max_step_shrink_factor\t%s\n" "$MAX_STEP_SHRINK_FACTOR"
  printf "use_soft_prior_for_profile\t%s\n" "$USE_SOFT_PRIOR_FOR_PROFILE"
  printf "lambda_prior_for_profile\t%s\n" "$LAMBDA_PRIOR_FOR_PROFILE"
  printf "parameter_count\t%s\n" "$PARAM_COUNT"
  printf "array_spec\t%s\n" "$ARRAY_SPEC"
  printf "memory\t%s\n" "$MEMORY"
  printf "walltime\t%s\n" "$WALLTIME"
  printf "job_name\t%s\n" "$JOB_NAME"
  printf "mail_user\t%s\n" "$MAIL_USER"
  printf "mail_type\t%s\n" "$MAIL_TYPE"
  printf "r_module\t%s\n" "$R_MODULE"
} > "$MANIFEST_PATH"

SBATCH_ARGS=(
  --array="$ARRAY_SPEC"
  --cpus-per-task="$N_CORES"
  --mem="$MEMORY"
  --time="$WALLTIME"
  --job-name="$JOB_NAME"
  --output="$OUTPUT_ROOT/slurm_logs/%x_%A_%a.out"
  --error="$OUTPUT_ROOT/slurm_logs/%x_%A_%a.err"
)

if [[ -n "$MAIL_USER" ]]; then
  SBATCH_ARGS+=(--mail-user="$MAIL_USER" --mail-type="$MAIL_TYPE")
fi

EXPORTS="ALL,PROJECT_ROOT=${PROJECT_ROOT},CONFIG_PATH=${CONFIG_PATH},BASELINE_SEED_DIR=${BASELINE_SEED_DIR},PROFILE_BOUNDS_TABLE=${PROFILE_BOUNDS_TABLE},OUTPUT_ROOT=${OUTPUT_ROOT},MAX_STEPS_PER_DIRECTION=${MAX_STEPS_PER_DIRECTION},SEEDS_PER_STEP=${SEEDS_PER_STEP},N_CORES=${N_CORES},TARGET_DELTA_OBJECTIVE=${TARGET_DELTA_OBJECTIVE},CI_DELTA_THRESHOLD=${CI_DELTA_THRESHOLD},STEP_FRACTION_INITIAL=${STEP_FRACTION_INITIAL},STEP_FRACTION_MIN=${STEP_FRACTION_MIN},STEP_FRACTION_MAX=${STEP_FRACTION_MAX},BOUNDARY_START_TOLERANCE=${BOUNDARY_START_TOLERANCE},MAX_ATTEMPTS_PER_STEP=${MAX_ATTEMPTS_PER_STEP},MIN_INTERIOR_POINTS_PER_DIRECTION=${MIN_INTERIOR_POINTS_PER_DIRECTION},CI_REFINE_TOLERANCE=${CI_REFINE_TOLERANCE},MAX_REFINE_STEPS_CI=${MAX_REFINE_STEPS_CI},MAX_REFINE_STEPS_BOUNDARY=${MAX_REFINE_STEPS_BOUNDARY},BOUNDARY_REFINE_FRACTION_MIN=${BOUNDARY_REFINE_FRACTION_MIN},MAX_STEP_GROWTH_FACTOR=${MAX_STEP_GROWTH_FACTOR},MAX_STEP_SHRINK_FACTOR=${MAX_STEP_SHRINK_FACTOR},USE_SOFT_PRIOR_FOR_PROFILE=${USE_SOFT_PRIOR_FOR_PROFILE},LAMBDA_PRIOR_FOR_PROFILE=${LAMBDA_PRIOR_FOR_PROFILE},R_MODULE=${R_MODULE}"

SBATCH_ARGS+=(--export="$EXPORTS")

SUB_SCRIPT="${SCRIPT_DIR}/submit_profile_likelihood_array.sub"
if [[ ! -f "$SUB_SCRIPT" ]]; then
  echo "Missing submit script: $SUB_SCRIPT" >&2
  exit 1
fi

echo "Submitting supplement-style profile likelihood array:"
echo "  output_root: $OUTPUT_ROOT"
echo "  parameter_count: $PARAM_COUNT"
echo "  array_spec: $ARRAY_SPEC"
echo "  max_steps_per_direction: $MAX_STEPS_PER_DIRECTION"
echo "  seeds_per_step: $SEEDS_PER_STEP"
echo "  n_cores: $N_CORES"

sbatch "${SBATCH_ARGS[@]}" "$SUB_SCRIPT"
