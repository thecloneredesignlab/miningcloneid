#!/usr/bin/env bash
set -euo pipefail

# Batch launcher:
# 1) two-stage fit (growth -> ploidy)
# 2) equal-weight callback warm-started from stage result
# 3) repeat over K grid and seed grid

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
FIT_SCRIPT="${SCRIPT_DIR}/fit_invivo_ploidy_buffer.R"

if [[ ! -f "${FIT_SCRIPT}" ]]; then
  echo "ERROR: fit script not found: ${FIT_SCRIPT}" >&2
  exit 1
fi

# ----------------------------
# Defaults (override via env)
# ----------------------------
SEEDS_CSV="${SEEDS_CSV:-1,2,3}"
K_GRID="${K_GRID:-1e10 1e11 1e12}"

STAGE1_W_BURDEN="${STAGE1_W_BURDEN:-1}"
STAGE1_W_PLOIDY="${STAGE1_W_PLOIDY:-0}"
STAGE2_W_BURDEN="${STAGE2_W_BURDEN:-0.3}"
STAGE2_W_PLOIDY="${STAGE2_W_PLOIDY:-0.7}"

ITERMAX="${ITERMAX:-200}"
NP="${NP:-140}"
N_STARTS="${N_STARTS:-60}"
OPTIM_MAXIT="${OPTIM_MAXIT:-12000}"
N_CORES="${N_CORES:-}"
MAX_SCENARIOS="${MAX_SCENARIOS:-}"

USE_DEOPTIM="${USE_DEOPTIM:-FALSE}"
DEOPTIM_PARALLEL="${DEOPTIM_PARALLEL:-FALSE}"
DOSE_ZERO_ONLY="${DOSE_ZERO_ONLY:-TRUE}"
TRUNCATE_AT_TREATMENT="${TRUNCATE_AT_TREATMENT:-FALSE}"
PLOIDY_AT_HARVEST="${PLOIDY_AT_HARVEST:-TRUE}"

OUT_ROOT="${OUT_ROOT:-${PROJECT_ROOT}/oxygen/results}"
RUN_PREFIX="${RUN_PREFIX:-fit_two_stage_callback_kgrid}"

# Prevent BLAS/OpenMP oversubscription on HPC by default.
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export VECLIB_MAXIMUM_THREADS="${VECLIB_MAXIMUM_THREADS:-1}"

mkdir -p "${OUT_ROOT}"
IFS=',' read -r -a SEEDS <<< "${SEEDS_CSV}"

echo "Running two-stage + equal-weight callback over K grid..."
echo "  Seeds: ${SEEDS_CSV}"
echo "  K grid: ${K_GRID}"
echo "  Stage1 weights: w_burden=${STAGE1_W_BURDEN}, w_ploidy=${STAGE1_W_PLOIDY}"
echo "  Stage2 weights: w_burden=${STAGE2_W_BURDEN}, w_ploidy=${STAGE2_W_PLOIDY}"
echo "  Output root: ${OUT_ROOT}"
echo

for K in ${K_GRID}; do
  k_tag="$(echo "${K}" | sed 's/[.+-]/_/g')"

  for seed_raw in "${SEEDS[@]}"; do
    seed="$(echo "${seed_raw}" | xargs)"
    [[ -z "${seed}" ]] && continue

    run_base="${OUT_ROOT}/${RUN_PREFIX}_K${k_tag}_seed${seed}"
    stage_dir="${run_base}/two_stage"
    callback_dir="${run_base}/callback_equal"
    stage_log="${stage_dir}.log"
    callback_log="${callback_dir}.log"
    mkdir -p "${run_base}"

    echo "[$(date '+%F %T')] K=${K}, seed=${seed}: two-stage fit"
    cmd_stage=(
      Rscript "${FIT_SCRIPT}"
      "--two_stage=TRUE"
      "--stage1_w_burden=${STAGE1_W_BURDEN}"
      "--stage1_w_ploidy=${STAGE1_W_PLOIDY}"
      "--stage2_w_burden=${STAGE2_W_BURDEN}"
      "--stage2_w_ploidy=${STAGE2_W_PLOIDY}"
      "--w_burden=1"
      "--w_ploidy=1"
      "--K=${K}"
      "--use_deoptim=${USE_DEOPTIM}"
      "--deoptim_parallel=${DEOPTIM_PARALLEL}"
      "--truncate_at_treatment=${TRUNCATE_AT_TREATMENT}"
      "--dose_zero_only=${DOSE_ZERO_ONLY}"
      "--ploidy_at_harvest=${PLOIDY_AT_HARVEST}"
      "--itermax=${ITERMAX}"
      "--NP=${NP}"
      "--n_starts=${N_STARTS}"
      "--optim_maxit=${OPTIM_MAXIT}"
      "--seed=${seed}"
      "--out_dir=${stage_dir}"
    )
    if [[ -n "${N_CORES}" ]]; then
      cmd_stage+=("--n_cores=${N_CORES}")
    fi
    if [[ -n "${MAX_SCENARIOS}" ]]; then
      cmd_stage+=("--max_scenarios=${MAX_SCENARIOS}")
    fi
    {
      echo "Command:"
      printf '  %q' "${cmd_stage[@]}"
      echo
      echo
      "${cmd_stage[@]}"
    } 2>&1 | tee "${stage_log}"

    echo "[$(date '+%F %T')] K=${K}, seed=${seed}: callback equal-weight fit"
    cmd_callback=(
      Rscript "${FIT_SCRIPT}"
      "--two_stage=FALSE"
      "--w_burden=1"
      "--w_ploidy=1"
      "--init_params_tsv=${stage_dir}/fit_parameter_stages.tsv"
      "--K=${K}"
      "--use_deoptim=${USE_DEOPTIM}"
      "--deoptim_parallel=${DEOPTIM_PARALLEL}"
      "--truncate_at_treatment=${TRUNCATE_AT_TREATMENT}"
      "--dose_zero_only=${DOSE_ZERO_ONLY}"
      "--ploidy_at_harvest=${PLOIDY_AT_HARVEST}"
      "--itermax=${ITERMAX}"
      "--NP=${NP}"
      "--n_starts=${N_STARTS}"
      "--optim_maxit=${OPTIM_MAXIT}"
      "--seed=${seed}"
      "--out_dir=${callback_dir}"
    )
    if [[ -n "${N_CORES}" ]]; then
      cmd_callback+=("--n_cores=${N_CORES}")
    fi
    if [[ -n "${MAX_SCENARIOS}" ]]; then
      cmd_callback+=("--max_scenarios=${MAX_SCENARIOS}")
    fi
    {
      echo "Command:"
      printf '  %q' "${cmd_callback[@]}"
      echo
      echo
      "${cmd_callback[@]}"
    } 2>&1 | tee "${callback_log}"
  done
done

echo "All K x seed runs completed."
