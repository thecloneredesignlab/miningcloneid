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
exec apptainer exec --cleanenv --containall --pwd "$FIGURE78_WORKSPACE" \
  "${binds[@]}" \
  --env "R_HOME=/opt/R/4.4.2/lib64/R" --env "LD_LIBRARY_PATH=${ld_path}" \
  --env "FIGURE_MODEL_CODE_ROOT=${FIGURE_MODEL_CODE_ROOT}" \
  --env "FIGURE_JOINT_RESULT_ROOT=${FIGURE_JOINT_RESULT_ROOT}" \
  --env "OPENBLAS_NUM_THREADS=1" --env "OMP_NUM_THREADS=1" --env "MKL_NUM_THREADS=1" \
  --env "R_MAKEVARS_USER=${script_dir}/figure6_compiler_makevars" \
  "$container" Rscript --vanilla "$@"
