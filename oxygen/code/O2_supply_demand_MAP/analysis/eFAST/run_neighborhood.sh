#!/usr/bin/env bash
set -euo pipefail
mode=${1:?Use audit, smoke, convergence, full, or pipeline}
case "$mode" in audit|smoke|convergence|full|pipeline) ;; *) exit 2 ;; esac
repo=${EFAST_REPO_ROOT:-/share/lab_crd/taoli/Project/soft_couping_eFAST}
fit=${EFAST_FIT_ROOT:-/share/lab_crd/taoli/Project/soft_couping_org/oxygen/results/fit_invivo_unified_np256_500seed_all_xxlarge_r442_exact_20260828_145253}
figure=${EFAST_FIGURE4_DIR:-/share/lab_crd/taoli/Project/HypoxiaLTEEFigures/revised/iteration4/data/Figures/Figure4}
sif=${EFAST_SIF:-/share/lab_crd/taoli/Docker/o2_supply_demand_map_r442_hpc_exact_salib152_20260923.sif}
workers=${EFAST_WORKERS:-16}
pilot_seeds=${EFAST_PILOT_SEEDS:-10}
code="$repo/oxygen/code/O2_supply_demand_MAP/analysis/eFAST"
out="$repo/oxygen/results/eFAST/neighborhood10pct"
[[ $(hostname -s) == hpctpa3pc0028 ]] || { echo "Run on hpctpa3pc0028" >&2; exit 1; }
[[ "$workers" =~ ^[1-9][0-9]*$ ]] || exit 2
[[ "$pilot_seeds" =~ ^([3-9]|10)$ ]] || exit 2
mkdir -p "$out"
# One controller per analysis root. The kernel releases this lock on exit.
exec 9>"$out/controller.lock"
flock -n 9 || { echo "Another neighborhood controller is running" >&2; exit 1; }
sif_hash=$(sha256sum "$sif" | awk '{print $1}')
[[ "$sif_hash" == 0b60f6cdab8a91f6bbeea1ab6cf02dd79660a295f24bc6e8ad6c5fed04d982ad ]] || exit 1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export APPTAINERENV_OMP_NUM_THREADS=1 APPTAINERENV_OPENBLAS_NUM_THREADS=1 APPTAINERENV_MKL_NUM_THREADS=1
container() {
  apptainer exec --cleanenv "$sif" env R_HOME= R_ENVIRON_USER=/dev/null \
    R_LIBS_USER=/opt/R/4.4.2/lib64/R/library OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 "$@"
}
stamp=$(date -u +%Y%m%dT%H%M%SZ)
{
  printf 'field\tvalue\n'
  printf 'git_commit\t%s\n' "$(git -C "$repo" rev-parse HEAD)"
  printf 'host\t%s\n' "$(hostname -s)"
  printf 'sif\t%s\n' "$sif"
  printf 'sif_sha256\t%s\n' "$sif_hash"
  printf 'workers\t%s\n' "$workers"
  printf 'convergence_policy\tdiagnostic_only\n'
  printf 'phase_repeats\t5\nmode\t%s\nstarted_utc\t%s\n' "$mode" "$stamp"
} > "$out/environment_${stamp}.tsv"
audit() { container python3 "$code/neighborhood.py" audit --fit-root "$fit" --figure4-dir "$figure" --out-dir "$out"; }
run() { container python3 "$code/neighborhood.py" run --fit-root "$fit" --figure4-dir "$figure" \
        --out-dir "$out" --workers "$workers" --pilot-seeds "$pilot_seeds" --stage "$1"; }
if [[ "$mode" == audit ]]; then
  audit
elif [[ "$mode" == pipeline ]]; then
  audit
  run smoke
  # Extend the global design to five phases, preserving verified R1/R2 designs.
  EFAST_WORKERS="$workers" bash "$code/run_efast.sh" full
  run convergence
  # Convergence is diagnostic only: always continue to all 500 audited endpoints.
  run full
else
  [[ -s "$out/input_manifest.json" ]] || audit
  run "$mode"
fi
echo "Neighborhood $mode complete at $(date -Is)"
