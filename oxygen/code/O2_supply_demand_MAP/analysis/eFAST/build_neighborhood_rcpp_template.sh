#!/usr/bin/env bash
set -euo pipefail
repo=${EFAST_REPO_ROOT:-/share/lab_crd/taoli/Project/soft_couping_eFAST}
sif=${EFAST_SIF:-/share/lab_crd/taoli/Docker/o2_supply_demand_map_r442_hpc_exact_salib152_20260923.sif}
out=${EFAST_OUT_DIR:-$repo/oxygen/results/eFAST/neighborhood10pct}
model="$repo/oxygen/code/O2_supply_demand_MAP/model"
mkdir -p "$out/slurm"
[[ ! -e "$out/slurm/rcpp_template.tar.gz" ]] || { echo "Rcpp template already exists" >&2; exit 1; }
build=$(mktemp -d "$out/slurm/.rcpp_build.XXXXXX")
trap 'rm -rf "$build"' EXIT
mkdir "$build/cache"
apptainer exec --cleanenv --bind "$build/cache:$model/.rcpp_cache_o2_supply_demand_map" "$sif" \
  env R_HOME= R_ENVIRON_USER=/dev/null R_LIBS_USER=/opt/R/4.4.2/lib64/R/library \
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 MININGCLONEID_RCPP_REBUILD=TRUE \
  Rscript --vanilla -e 'args <- commandArgs(TRUE); source(args[[1]]); stopifnot(isTRUE(.USE_CPP_O2SIMPS_BACKEND)); cat("Validated C++ backend initialized\n")' \
  "$model/model_O2_supply_demand_MAP.R"
tar -czf "$out/slurm/rcpp_template.tar.gz" -C "$build/cache" .
apptainer exec --cleanenv "$sif" env OPENBLAS_NUM_THREADS=1 python3 -c '
import datetime,hashlib,json,pathlib,sys
def sha(p):
 h=hashlib.sha256()
 with open(p,"rb") as f:
  for b in iter(lambda:f.read(1048576),b""):h.update(b)
 return h.hexdigest()
out,model,sif=map(pathlib.Path,sys.argv[1:])
archive=out/"slurm/rcpp_template.tar.gz"
data=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),sif_sha256=sha(sif),model_cpp_sha256=sha(model/"model_O2_supply_demand_MAP.cpp"),model_r_sha256=sha(model/"model_O2_supply_demand_MAP.R"),archive_sha256=sha(archive),archive_bytes=archive.stat().st_size,archive_mtime_ns=archive.stat().st_mtime_ns,logical_bind_target=str(model/".rcpp_cache_o2_supply_demand_map"),runtime_rebuild=False,cache_scope="unique node-local directory per Slurm task")
(out/"slurm/rcpp_template_manifest.json").write_text(json.dumps(data,indent=2)+"\n")
print(json.dumps(data,indent=2))' "$out" "$model" "$sif"
