#!/usr/bin/env bash
set -euo pipefail

mode=${1:?Use pilot or full}
case "$mode" in pilot|full) ;; *) echo "Use pilot or full" >&2; exit 2 ;; esac

repo_root=${EFAST_REPO_ROOT:-/share/lab_crd/taoli/Project/soft_couping_eFAST}
fit_root=${EFAST_FIT_ROOT:-/share/lab_crd/taoli/Project/soft_couping_org/oxygen/results/fit_invivo_unified_np256_500seed_all_xxlarge_r442_exact_20260828_145253}
figure4_dir=${EFAST_FIGURE4_DIR:-/share/lab_crd/taoli/Project/HypoxiaLTEEFigures/revised/iteration4/data/Figures/Figure4}
sif=${EFAST_SIF:-/share/lab_crd/taoli/Docker/o2_supply_demand_map_r442_hpc_exact_salib152_20260923.sif}
workers=${EFAST_WORKERS:-4}
code_dir="$repo_root/oxygen/code/O2_supply_demand_MAP/analysis/eFAST"
out_dir="$repo_root/oxygen/results/eFAST"
if [[ "$mode" == pilot ]]; then out_dir="$out_dir/pilot"; fi
mkdir -p "$out_dir"

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
echo "host=$(hostname) mode=$mode workers=$workers repo=$repo_root sif=$sif"
if [[ $(hostname -s) != hpctpa3pc0028 ]]; then
  echo "Refusing to run outside hpctpa3pc0028" >&2
  exit 1
fi
sif_sha256=$(sha256sum "$sif" | awk '{print $1}')
if [[ "$sif_sha256" != 0b60f6cdab8a91f6bbeea1ab6cf02dd79660a295f24bc6e8ad6c5fed04d982ad ]]; then
  echo "SIF SHA-256 differs from the validated SALib 1.5.2 image" >&2
  exit 1
fi

{
  printf 'field\tvalue\n'
  printf 'git_commit\t%s\n' "$(git -C "$repo_root" rev-parse HEAD)"
  printf 'host\t%s\n' "$(hostname -s)"
  printf 'sif\t%s\n' "$sif"
  printf 'sif_sha256\t%s\n' "$sif_sha256"
  printf 'fit_root\t%s\n' "$fit_root"
  printf 'figure4_dir\t%s\n' "$figure4_dir"
  printf 'workers\t%s\n' "$workers"
  printf 'run_started\t%s\n' "$(date -Is)"
} > "$out_dir/environment.tsv"

container() {
  apptainer exec --cleanenv "$sif" env \
    R_HOME= R_ENVIRON_USER=/dev/null \
    R_LIBS_USER=/opt/R/4.4.2/lib64/R/library "$@"
}
if [[ "$mode" == pilot ]]; then
  resolutions=(129)
  replicates=1
  grid_flag=()
else
  resolutions=(129 257)
  replicates=2
  grid_flag=(--full-grid)
fi

for n in "${resolutions[@]}"; do
  container python3 "$code_dir/efast.py" prepare \
    --fit-root "$fit_root" --figure4-dir "$figure4_dir" --out-dir "$out_dir" \
    --n "$n" --replicates "$replicates" "${grid_flag[@]}"
done

for n in "${resolutions[@]}"; do
  for ((rep=1; rep<=replicates; rep++)); do
    run_dir="$out_dir/runs/N${n}_R${rep}"
    if [[ ! -s "$run_dir/outputs.tsv.gz" ]]; then
      container Rscript --vanilla "$code_dir/evaluate_fixed_o2.R" \
        --samples="$run_dir/samples.tsv.gz" --metadata="$run_dir/metadata.json" \
        --fit_root="$fit_root" --figure4_dir="$figure4_dir" \
        --out="$run_dir/outputs.tsv.gz" --workers="$workers" --validate=TRUE
    fi
  done
done

container python3 "$code_dir/efast.py" summarize --out-dir "$out_dir"
if [[ "$mode" == full ]]; then
  container python3 "$code_dir/efast.py" plot --out-dir "$out_dir" --figure4-dir "$figure4_dir"
fi
echo "eFAST $mode complete at $(date -Is)"
