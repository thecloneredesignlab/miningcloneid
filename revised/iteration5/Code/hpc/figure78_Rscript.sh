#!/usr/bin/env bash
# Common interpreter launcher; no sbatch, preparation, or implicit simulation.
set -euo pipefail
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
source "${script_dir}/../config/figure78_hpc_paths.env"
if [[ -n "${FIGURE78_NATIVE_RSCRIPT:-}" ]]; then
  exec "${FIGURE78_NATIVE_RSCRIPT}" "$@"
fi
if [[ "${FIGURE78_WORKSPACE}" != /share/* && -n "$(command -v Rscript || true)" ]]; then
  exec Rscript "$@"
fi
container="${FIGURE78_CONTAINER:-/share/lab_crd/taoli/Docker/o2_supply_demand_map_r442_hpc_exact.sif}"
[[ -f "$container" ]] || { echo "Set FIGURE78_NATIVE_RSCRIPT or an approved FIGURE78_CONTAINER" >&2; exit 2; }
command -v apptainer >/dev/null || { echo "Apptainer is required on HPC" >&2; exit 2; }
ld_path=/app/eb/software/FlexiBLAS/3.4.4-GCC-13.3.0/lib64:/app/eb/software/OpenBLAS/0.3.27-GCC-13.3.0/lib:/app/eb/software/GCCcore/13.3.0/lib64:/app/eb/software/binutils/2.42-GCCcore-13.3.0/lib:/opt/rh/gcc-toolset-13/root/usr/lib64
binds=(--bind "${FIGURE78_WORKSPACE}:${FIGURE78_WORKSPACE}:rw"
       --bind "${FIGURE_MODEL_CODE_ROOT}:${FIGURE_MODEL_CODE_ROOT}:ro"
       --bind /app/eb:/app/eb:ro)
# Retained-data drawing needs no original fit files. Generation itself refuses
# missing original inputs, but a drawing-only container need not mount them.
if [[ -d "${FIGURE_JOINT_RESULT_ROOT}" ]]; then
  binds+=(--bind "${FIGURE_JOINT_RESULT_ROOT}:${FIGURE_JOINT_RESULT_ROOT}:ro")
fi
# Runtime-only Figure7 repair: keep the authenticated numerical code and its
# logical source/cache paths unchanged, but mount a distinct physical cache for
# each Slurm process. The locked model hard-codes its cache location; setting an
# unrecognized cache environment variable does not isolate its directory lock.
# Figure8 and local/native runs retain the original launch path.
figure7_cache_task_dir=
if [[ "${FIGURE78_FIGURE:-}" == 7 && -n "${SLURM_JOB_ID:-}" && -n "${FIGURE78_DATA_ROOT:-}" ]]; then
  [[ "$FIGURE78_DATA_ROOT" == /* ]] || { echo "Figure7 cache isolation requires an absolute data root" >&2; exit 2; }
  figure7_run="${FIGURE78_DATA_ROOT}/Figure7"
  [[ -f "${figure7_run}/prepared/analysis_bundle.rds" && -f "${figure7_run}/prepared/model_source/model/model_O2_supply_demand_MAP.cpp" ]] || {
    echo "Figure7 cache isolation requires the authenticated prepared model" >&2; exit 2;
  }
  figure7_model_cache="${figure7_run}/prepared/model_source/model/.rcpp_cache_o2_supply_demand_map"
  figure7_propagator_cache="${figure7_run}/.cache/propagators"
  mkdir -p "$figure7_model_cache" "$figure7_propagator_cache"
  figure7_cache_task_dir="$(mktemp -d /tmp/figure78-cache.XXXXXXXXXX)"
  cleanup_figure7_task_cache() {
    # Only this process's securely-created /tmp directory is disposable. Never
    # remove a shared cache, prepared model, result, or user-selected directory.
    python3 - "$figure7_cache_task_dir" <<'PY'
from pathlib import Path
import shutil
import sys
p=Path(sys.argv[1])
if p.parent!=Path('/tmp') or not p.name.startswith('figure78-cache.') or p.is_symlink():
    raise SystemExit('Refusing unsafe task-cache cleanup')
if p.is_dir():shutil.rmtree(str(p))
PY
  }
  trap cleanup_figure7_task_cache EXIT
  mkdir -p "${figure7_cache_task_dir}/model" "${figure7_cache_task_dir}/propagators"
  rsync -rt --exclude=.sourcecpp_lock "${figure7_model_cache}/" "${figure7_cache_task_dir}/model/"
  rsync -rt "${figure7_propagator_cache}/" "${figure7_cache_task_dir}/propagators/"
  binds+=(--bind "${figure7_cache_task_dir}/model:${figure7_model_cache}:rw"
          --bind "${figure7_cache_task_dir}/propagators:${figure7_propagator_cache}:rw")
  echo "FIGURE7_CACHE_ISOLATED job=${SLURM_JOB_ID} task=${SLURM_ARRAY_TASK_ID:-driver} physical=${figure7_cache_task_dir}"
fi
figure78_command=(apptainer exec --cleanenv --containall --pwd "$FIGURE78_WORKSPACE" \
  "${binds[@]}" \
  --env "R_HOME=/opt/R/4.4.2/lib64/R" --env "LD_LIBRARY_PATH=${ld_path}" \
  --env "FIGURE_MODEL_CODE_ROOT=${FIGURE_MODEL_CODE_ROOT}" \
  --env "FIGURE_JOINT_RESULT_ROOT=${FIGURE_JOINT_RESULT_ROOT}" \
  --env "OPENBLAS_NUM_THREADS=1" --env "OMP_NUM_THREADS=1" --env "MKL_NUM_THREADS=1" \
  --env "R_MAKEVARS_USER=${script_dir}/figure6_compiler_makevars" \
  "$container" Rscript --vanilla "$@")
if [[ -n "$figure7_cache_task_dir" ]]; then
  "${figure78_command[@]}"
else
  exec "${figure78_command[@]}"
fi
