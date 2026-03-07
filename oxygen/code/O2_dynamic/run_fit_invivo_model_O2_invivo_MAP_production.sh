#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUNNER="${SCRIPT_DIR}/run_fit_invivo_model_O2_invivo_MAP_it_stagewise_chain_callback.sh"

if [[ ! -f "${RUNNER}" ]]; then
  echo "ERROR: runner not found: ${RUNNER}" >&2
  exit 1
fi

OUT_ROOT_DEFAULT="$(cd "${SCRIPT_DIR}/../.." && pwd)/results"
DATA_DIR_DEFAULT="$(cd "${SCRIPT_DIR}/../../.." && pwd)/data/InVivoData_Gemcitabine"

RUN_PREFIX="${RUN_PREFIX:-fit_stagewise_chain_o2invivo_map_$(date +%Y%m%d_%H%M%S)}"
OUT_ROOT="${OUT_ROOT:-${OUT_ROOT_DEFAULT}}"
DATA_DIR="${DATA_DIR:-${DATA_DIR_DEFAULT}}"
SEEDS_FILE="${SEEDS_FILE:-${SCRIPT_DIR}/seeds.csv}"
N_CORES="${N_CORES:-32}"

bash "${RUNNER}" \
  --run_prefix="${RUN_PREFIX}" \
  --out_root="${OUT_ROOT}" \
  --data_dir="${DATA_DIR}" \
  --seeds_file="${SEEDS_FILE}" \
  --n_cores="${N_CORES}" \
  --use_deoptim=TRUE \
  --deoptim_parallel=TRUE \
  --fit_treatment=FALSE \
  --paired_only=TRUE \
  --auto_tune_iters=FALSE \
  --w_burden_chain=1,0.8,0.6,0.4,0.3,0.2,0.1,0 \
  --w_ploidy_chain=0,0.2,0.4,0.6,0.7,0.8,0.9,1 \
  --callback_w_burden=1 \
  --callback_w_ploidy=1 \
  --callback_init_pass=8 \
  --use_soft_prior=TRUE \
  --lambda_prior=0.03 \
  --burden_exclude_day0=TRUE \
  --de_init_mode=hybrid \
  --de_init_uniform_frac=0.5 \
  --de_init_sigma_frac=0.1 \
  --de_reltol=1e-3 \
  --de_steptol=25 \
  --o2_cache_bin_pct=0.01 \
  --o2_cache_hysteresis_pct=0.005 \
  --o2_cache_profile=FALSE \
  --tau_O2=3 \
  --pass_itermax=220 \
  --callback_itermax=480 \
  --np=256 \
  --sigma_burden=0.35 \
  --sigma_ploidy=0.15 \
  "$@"

