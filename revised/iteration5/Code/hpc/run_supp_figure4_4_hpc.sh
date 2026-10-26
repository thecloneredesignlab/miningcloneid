#!/usr/bin/env bash
#SBATCH --job-name=supp4_4_500seeds
#SBATCH --partition=red
#SBATCH --qos=xxlarge
#SBATCH --nodelist=hpctpa3pc0028
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
set -Eeuo pipefail
ITERATION_ROOT="/share/lab_crd/taoli/Project/HypoxiaLTEEFigures/revised/iteration5"
FIT_ROOT="/share/lab_crd/taoli/Project/soft_couping_org/oxygen/results/fit_invivo_unified_np256_500seed_all_xxlarge_r442_exact_20260828_145253"
SIF_IMAGE="/share/lab_crd/taoli/Docker/o2_supply_demand_map_r442_hpc_exact.sif"
[[ "$(hostname -s)" == hpctpa3pc0028 ]] || { echo "Wrong compute node" >&2; exit 2; }
: "${SLURM_JOB_ID:?Submit this script through Slurm}"
: "${SUPP4_4_CODE_COMMIT:?Provide verified source commit}"
cd "${ITERATION_ROOT}"
RUN_ROOT="${ITERATION_ROOT}/audit/hpc_supp_figure4_4/${SLURM_JOB_ID}"
mkdir -p "${RUN_ROOT}/tmp/home" "${RUN_ROOT}/tmp/cache"
exec > >(tee -a "${RUN_ROOT}/run.log") 2>&1
stage=START
status() {
  printf 'job_id\tstatus\tstage\thost\tcode_commit\tupdated_at\n' > "${RUN_ROOT}/status.tsv"
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "${SLURM_JOB_ID}" "$1" "${stage}" "$(hostname -s)" "${SUPP4_4_CODE_COMMIT}" "$(TZ=America/New_York date -Iseconds)" >> "${RUN_ROOT}/status.tsv"
}
trap 'code=$?; status FAILED; exit "${code}"' ERR
status RUNNING
stage=VERIFY_TASK_CODE; status RUNNING
for path in \
  Code/Figures/util/analysis/supp_figure4_4.R \
  Code/Figures/test_Supp_Figure4_4.R \
  Code/Figures/data_Supp_Figure4_4.R \
  Code/Figures/draw_Supp_Figure4_4.R \
  Code/Figures/validate_Supp_Figure4_4.py \
  Code/Figures/Supp_Figure4_4_notes.md \
  Code/config/supp_figure4_4.Rprofile \
  Code/hpc/run_supp_figure4_4_hpc.sh; do
  expected=$(git rev-parse "${SUPP4_4_CODE_COMMIT}:revised/iteration5/${path}")
  actual=$(git hash-object "${path}")
  [[ "${actual}" == "${expected}" ]] || { echo "Code mismatch: ${path}" >&2; exit 2; }
  printf '%s\t%s\n' "${actual}" "${path}" >> "${RUN_ROOT}/code_blobs.tsv"
done
if command -v apptainer >/dev/null 2>&1; then runtime=$(command -v apptainer); else runtime=$(command -v singularity); fi
RED_LIB=/app/eb/software
LD_PATH="${RED_LIB}/FlexiBLAS/3.4.4-GCC-13.3.0/lib64:${RED_LIB}/OpenBLAS/0.3.27-GCC-13.3.0/lib:${RED_LIB}/GCCcore/13.3.0/lib64:${RED_LIB}/binutils/2.42-GCCcore-13.3.0/lib:${RED_LIB}/libiconv/1.17-GCCcore-13.3.0/lib:${RED_LIB}/ICU/75.1-GCCcore-13.3.0/lib:/opt/rh/gcc-toolset-13/root/usr/lib64"
args=(exec --cleanenv --containall --pwd "${ITERATION_ROOT}" --home "${RUN_ROOT}/tmp/home"
  --env "R_HOME=/opt/R/4.4.2/lib64/R" --env "LD_LIBRARY_PATH=${LD_PATH}"
  --env "LD_PRELOAD=${RED_LIB}/libiconv/1.17-GCCcore-13.3.0/lib/libiconv.so.2"
  --env TMPDIR=/tmp --env TMP=/tmp --env TEMP=/tmp
  --env "XDG_CACHE_HOME=${RUN_ROOT}/tmp/cache"
  --env "R_PROFILE_USER=${ITERATION_ROOT}/Code/config/supp_figure4_4.Rprofile"
  --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 --env MKL_NUM_THREADS=1
  --env RCPP_PARALLEL_NUM_THREADS=1 --env KMP_USE_SHM=0
  --env "SUPP4_4_FIT_ROOT=${FIT_ROOT}"
  --bind /app/eb:/app/eb:ro --bind "${ITERATION_ROOT}:${ITERATION_ROOT}:rw"
  --bind "${FIT_ROOT}:${FIT_ROOT}:ro" --bind "${RUN_ROOT}/tmp:/tmp:rw" --bind "${RUN_ROOT}/tmp:/var/tmp:rw")
container() { "${runtime}" "${args[@]}" "${SIF_IMAGE}" "$@"; }
find "${ITERATION_ROOT}/Figures" -maxdepth 1 -type f \( -name 'Figure4B*' -o -name 'assembled_fig4.*' -o -name 'Supp_Figure4_3.*' -o -name 'supp_fig4-1*' -o -name 'supp_fig4-2*' \) -print0 |
  sort -z | xargs -0 sha256sum > "${RUN_ROOT}/unchanged_figure4_before.sha256"
stage=UNIT_TESTS; status RUNNING
container Rscript "${ITERATION_ROOT}/Code/Figures/test_Supp_Figure4_4.R"
stage=VALIDATE_CACHES_AND_AGGREGATE; status RUNNING
container Rscript "${ITERATION_ROOT}/Code/Figures/data_Supp_Figure4_4.R"
stage=RENDER; status RUNNING
container Rscript "${ITERATION_ROOT}/Code/Figures/draw_Supp_Figure4_4.R"
stage=VERIFY_EXISTING_FIGURES_UNCHANGED; status RUNNING
sha256sum -c "${RUN_ROOT}/unchanged_figure4_before.sha256" > "${RUN_ROOT}/unchanged_figure4_check.txt"
stage=OUTPUT_CHECKSUMS; status RUNNING
find "${ITERATION_ROOT}/data/Figures/Supp_Figure4_4" -maxdepth 1 -type f -print0 |
  sort -z | xargs -0 sha256sum > "${RUN_ROOT}/output_sha256.txt"
sha256sum "${ITERATION_ROOT}/Figures/Supp_Figure4_4.pdf" "${ITERATION_ROOT}/Figures/Supp_Figure4_4.png" "${ITERATION_ROOT}/Figures/Supp_Figure4_4.svg" >> "${RUN_ROOT}/output_sha256.txt"
stage=COMPLETE; status COMPLETE
echo "supp_figure4_4_complete"
