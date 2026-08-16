# RED HPC-exact SIF workflow

This directory contains the build and validation path used when numerical
reproduction requires the exact RED EasyBuild R/4.4.2 runtime rather than a
source reconstruction with matching package versions.

The accepted host reference is the in-vitro seed10 run produced with:

- R `/app/eb/software/R/4.4.2-gfbf-2024a`
- FlexiBLAS 3.4.4 with OpenBLAS 0.3.27 and LAPACK 3.12.0
- GCC 13.3.0 and binutils 2.42
- 22 process workers and one BLAS/OpenMP thread per worker

## Stages

1. `run_seed10_bridge_validation.sh` runs inside the existing SIF namespace
   but bind-mounts the exact RED EasyBuild tree and R package library. It uses
   an isolated empty `sourceCpp` cache so the model DLL is rebuilt by the RED
   compiler instead of reusing a pre-existing object.
2. `compare_seed10_results.sh` requires byte-identical raw/transformed
   parameter tables and exact DE/local/total objective strings.
3. `build_hpc_exact_sif.sh` copies every prefix in the captured module stack
   and the frozen R package library into a new standalone image. It never
   overwrites an existing sandbox or SIF.
4. `run_seed10_standalone_validation.sh` runs the new image with
   `--cleanenv --containall`, binds only the project and an isolated C++ cache,
   regenerates the HTML report, and repeats the exact comparison.

`o2_hpc_exact_rscript` is the standalone image entrypoint. It restores the
captured RED `PATH` and `LD_LIBRARY_PATH`, selects the internal R package
library, and invokes the EasyBuild `Rscript` binary stored in the image.

## Acceptance contract

The standalone run is accepted only when all of the following hold:

- `optimizer_deoptim_objective` is exactly `4.22415977007025`.
- `optimizer_local_objective` and `objective_total` are exactly
  `3.85253526260594`.
- `best_params.tsv` and `best_params_transformed.tsv` are byte-identical to
  the canonical seed10 reference.
- the fit, post-fit pipeline, and non-empty `fit_report_seed10.html` complete.
- `.libPaths()` begins with `/opt/o2-host-r-library/4.4`; the host user R
  library is not mounted into the standalone run.

The serialized `fit_result.rds` and HTML SHA-256 are not required to match the
historical artifacts because current code records additional Death metadata,
paths, and report-generation details. Numerical outputs and parameter tables
are the exact-reproduction boundary.
