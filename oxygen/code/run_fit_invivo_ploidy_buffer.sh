#!/usr/bin/env bash
set -euo pipefail

# Multi-seed launcher for fit_invivo_ploidy_buffer.R
# Default strategy:
# - single-stage joint fitting
# - progressive weight schedule (ploidy-prioritized -> balanced)
# - optim multi-start backend

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
FIT_SCRIPT="${SCRIPT_DIR}/fit_invivo_ploidy_buffer.R"

# ----------------------------
# Defaults (override via env)
# ----------------------------
SEEDS_CSV="${SEEDS_CSV:-1,2,3}"
W_BURDEN="${W_BURDEN:-0.05,0.2,1}"
W_PLOIDY="${W_PLOIDY:-8,3,1}"

ITERMAX="${ITERMAX:-200}"
NP="${NP:-140}"
N_STARTS="${N_STARTS:-60}"
OPTIM_MAXIT="${OPTIM_MAXIT:-12000}"
N_CORES="${N_CORES:-}"          # empty => auto detect in R script
MAX_SCENARIOS="${MAX_SCENARIOS:-}"  # optional

DOSE_ZERO_ONLY="${DOSE_ZERO_ONLY:-TRUE}"
TRUNCATE_AT_TREATMENT="${TRUNCATE_AT_TREATMENT:-FALSE}"
PLOIDY_AT_HARVEST="${PLOIDY_AT_HARVEST:-TRUE}"

USE_DEOPTIM="${USE_DEOPTIM:-FALSE}"
DEOPTIM_PARALLEL="${DEOPTIM_PARALLEL:-FALSE}"

OUT_ROOT="${OUT_ROOT:-${PROJECT_ROOT}/oxygen/results}"
RUN_PREFIX="${RUN_PREFIX:-fit_joint_sched_multi_seed}"

mkdir -p "${OUT_ROOT}"

IFS=',' read -r -a SEEDS <<< "${SEEDS_CSV}"

echo "Running multi-seed fit_invivo_ploidy_buffer jobs..."
echo "  Seeds: ${SEEDS_CSV}"
echo "  Weight schedule: w_burden=${W_BURDEN} | w_ploidy=${W_PLOIDY}"
echo "  Output root: ${OUT_ROOT}"
echo

for seed in "${SEEDS[@]}"; do
  seed="$(echo "${seed}" | xargs)"
  [[ -z "${seed}" ]] && continue

  out_base="${OUT_ROOT}/${RUN_PREFIX}_seed${seed}"
  log_file="${out_base}.log"

  echo "[$(date '+%F %T')] Start seed=${seed}"
  echo "  Log: ${log_file}"

  cmd=(
    Rscript "${FIT_SCRIPT}"
    "--two_stage=FALSE"
    "--use_deoptim=${USE_DEOPTIM}"
    "--deoptim_parallel=${DEOPTIM_PARALLEL}"
    "--truncate_at_treatment=${TRUNCATE_AT_TREATMENT}"
    "--dose_zero_only=${DOSE_ZERO_ONLY}"
    "--ploidy_at_harvest=${PLOIDY_AT_HARVEST}"
    "--w_burden=${W_BURDEN}"
    "--w_ploidy=${W_PLOIDY}"
    "--itermax=${ITERMAX}"
    "--NP=${NP}"
    "--n_starts=${N_STARTS}"
    "--optim_maxit=${OPTIM_MAXIT}"
    "--seed=${seed}"
    "--out_dir=${out_base}"
    "--append_timestamp_out_dir=TRUE"
  )

  if [[ -n "${N_CORES}" ]]; then
    cmd+=("--n_cores=${N_CORES}")
  fi
  if [[ -n "${MAX_SCENARIOS}" ]]; then
    cmd+=("--max_scenarios=${MAX_SCENARIOS}")
  fi

  {
    echo "Command:"
    printf '  %q' "${cmd[@]}"
    echo
    echo
    "${cmd[@]}"
  } 2>&1 | tee "${log_file}"

  echo "[$(date '+%F %T')] Done seed=${seed}"
  echo
done

echo "All seed jobs completed."
