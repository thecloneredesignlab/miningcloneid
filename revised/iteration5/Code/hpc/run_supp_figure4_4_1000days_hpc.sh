#!/usr/bin/env bash
#SBATCH --job-name=supp4_4_1000days
#SBATCH --partition=red
#SBATCH --qos=xxlarge
#SBATCH --nodelist=hpctpa3pc0028
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
set -Eeuo pipefail
ITERATION_ROOT=/share/lab_crd/taoli/Project/HypoxiaLTEEFigures/revised/iteration5
FIT_ROOT=/share/lab_crd/taoli/Project/soft_couping_org/oxygen/results/fit_invivo_unified_np256_500seed_all_xxlarge_r442_exact_20260828_145253
SIF_IMAGE=/share/lab_crd/taoli/Docker/o2_supply_demand_map_r442_hpc_exact.sif
EXPECTED_NODE="${SUPP4_4_EXPECTED_NODE:-hpctpa3pc0028}"
[[ "$(hostname -s)" == "${EXPECTED_NODE}" ]] || { echo "Wrong compute node: expected ${EXPECTED_NODE}" >&2; exit 2; }
: "${SLURM_JOB_ID:?Submit through Slurm}"
: "${SUPP4_4_CODE_COMMIT:?Provide local source commit}"
: "${SUPP4_4_CODE_CHECKSUMS:?Provide eight-file SHA256 manifest}"
[[ "${SUPP4_4_CODE_CHECKSUMS}" == "${ITERATION_ROOT}/audit/"* ]] || exit 2
cd "${ITERATION_ROOT}"
RUN_ROOT="${ITERATION_ROOT}/audit/hpc_supp_figure4_4_1000days/${SLURM_JOB_ID}"
mkdir -p "${RUN_ROOT}/tmp/home" "${RUN_ROOT}/tmp/cache"
exec > >(tee -a "${RUN_ROOT}/run.log") 2>&1
stage=START
status() {
  printf 'job_id\tstatus\tstage\thost\tcode_commit\tupdated_at\n' > "${RUN_ROOT}/status.tsv"
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "${SLURM_JOB_ID}" "$1" "${stage}" "$(hostname -s)" "${SUPP4_4_CODE_COMMIT}" "$(TZ=America/New_York date -Iseconds)" >> "${RUN_ROOT}/status.tsv"
}
trap 'code=$?; status FAILED; exit "${code}"' ERR
stage=VERIFY_TASK_CODE; status RUNNING
[[ "$(wc -l < "${SUPP4_4_CODE_CHECKSUMS}")" -eq 8 ]] || exit 2
sha256sum -c "${SUPP4_4_CODE_CHECKSUMS}" | tee "${RUN_ROOT}/code_checks.txt"
cp "${SUPP4_4_CODE_CHECKSUMS}" "${RUN_ROOT}/code_sha256.txt"
if command -v apptainer >/dev/null 2>&1; then runtime=$(command -v apptainer); else runtime=$(command -v singularity); fi
RED_LIB=/app/eb/software
LD_PATH="${RED_LIB}/FlexiBLAS/3.4.4-GCC-13.3.0/lib64:${RED_LIB}/OpenBLAS/0.3.27-GCC-13.3.0/lib:${RED_LIB}/GCCcore/13.3.0/lib64:${RED_LIB}/binutils/2.42-GCCcore-13.3.0/lib:${RED_LIB}/libiconv/1.17-GCCcore-13.3.0/lib:${RED_LIB}/ICU/75.1-GCCcore-13.3.0/lib:/opt/rh/gcc-toolset-13/root/usr/lib64"
args=(exec --cleanenv --containall --pwd "${ITERATION_ROOT}" --home "${RUN_ROOT}/tmp/home"
  --env "R_HOME=/opt/R/4.4.2/lib64/R" --env "LD_LIBRARY_PATH=${LD_PATH}"
  --env "LD_PRELOAD=${RED_LIB}/libiconv/1.17-GCCcore-13.3.0/lib/libiconv.so.2"
  --env TMPDIR=/tmp --env TMP=/tmp --env TEMP=/tmp --env "XDG_CACHE_HOME=${RUN_ROOT}/tmp/cache"
  --env "R_PROFILE_USER=${ITERATION_ROOT}/Code/config/supp_figure4_4.Rprofile"
  --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 --env MKL_NUM_THREADS=1
  --env RCPP_PARALLEL_NUM_THREADS=1 --env KMP_USE_SHM=0 --env "SUPP4_4_FIT_ROOT=${FIT_ROOT}"
  --bind /app/eb:/app/eb:ro --bind "${ITERATION_ROOT}:${ITERATION_ROOT}:rw"
  --bind "${FIT_ROOT}:${FIT_ROOT}:ro" --bind "${RUN_ROOT}/tmp:/tmp:rw" --bind "${RUN_ROOT}/tmp:/var/tmp:rw")
container() { "${runtime}" "${args[@]}" "${SIF_IMAGE}" "$@"; }
find "${ITERATION_ROOT}/Figures" -maxdepth 1 -type f \( -name 'Figure4B*' -o -name 'assembled_fig4.*' -o -name 'Supp_Figure4_3.*' -o -name 'supp_fig4-1*' -o -name 'supp_fig4-2*' -o -name 'Supp_Figure4_4.*' \) -print0 |
  sort -z | xargs -0 sha256sum > "${RUN_ROOT}/protected_before.sha256"
find "${ITERATION_ROOT}/data/Figures/Supp_Figure4_4" -maxdepth 1 -type f -print0 |
  sort -z | xargs -0 sha256sum >> "${RUN_ROOT}/protected_before.sha256"
stage=UNIT_TESTS; status RUNNING
container Rscript "${ITERATION_ROOT}/Code/Figures/test_Supp_Figure4_4.R"
if [[ "${SUPP4_4_RENDER_ONLY:-0}" == 1 ]]; then
  stage=REUSE_VALIDATED_DATA; status RUNNING
  : "${SUPP4_4_DATA_CHECKSUMS:?Provide checksums of validated analysis tables}"
  [[ "${SUPP4_4_DATA_CHECKSUMS}" == "${ITERATION_ROOT}/audit/"* ]] || exit 2
  [[ "$(wc -l < "${SUPP4_4_DATA_CHECKSUMS}")" -eq 9 ]] || exit 2
  sha256sum -c "${SUPP4_4_DATA_CHECKSUMS}" | tee "${RUN_ROOT}/reused_data_checks.txt"
  cp "${SUPP4_4_DATA_CHECKSUMS}" "${RUN_ROOT}/reused_data_sha256.txt"
else
  stage=VALIDATE_1000DAY_CACHES_AND_OVERLAP; status RUNNING
  container Rscript "${ITERATION_ROOT}/Code/Figures/data_Supp_Figure4_4_1000days.R"
fi
stage=RENDER_1000DAY_VERSION; status RUNNING
container Rscript "${ITERATION_ROOT}/Code/Figures/draw_Supp_Figure4_4_1000days.R"
if [[ "${SUPP4_4_RENDER_ONLY:-0}" == 1 ]]; then
  sha256sum -c "${SUPP4_4_DATA_CHECKSUMS}" > "${RUN_ROOT}/reused_data_after_render.txt"
fi
stage=VERIFY_ORIGINAL_FIGURES_AND_DATA; status RUNNING
sha256sum -c "${RUN_ROOT}/protected_before.sha256" > "${RUN_ROOT}/protected_check.txt"
stage=OUTPUT_CHECKSUMS; status RUNNING
find "${ITERATION_ROOT}/data/Figures/Supp_Figure4_4_1000days" -maxdepth 1 -type f -print0 |
  sort -z | xargs -0 sha256sum > "${RUN_ROOT}/output_sha256.txt"
sha256sum "${ITERATION_ROOT}/Figures/Supp_Figure4_4_1000days.pdf" "${ITERATION_ROOT}/Figures/Supp_Figure4_4_1000days.png" "${ITERATION_ROOT}/Figures/Supp_Figure4_4_1000days.svg" >> "${RUN_ROOT}/output_sha256.txt"
stage=COMPLETE; status COMPLETE
echo supp_figure4_4_1000days_complete
