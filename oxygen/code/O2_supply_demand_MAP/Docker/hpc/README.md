# Image-backed HPC layer

This directory is the one-to-one Apptainer counterpart of `../../hpc/`.
Submission remains on the Slurm host, with the same job resources, arrays,
dependencies, log paths, called workflow entrypoints, and command-line
interfaces. Only the language runtime changes: R and Python execute inside the
prebuilt container image instead of loading host modules.

## Runtime

The default immutable HPC artifact is:

```text
/share/lab_crd/lab_crd/taoli/Docker/o2_supply_demand_map_r44.sif
```

It was produced from Docker image `zafiro/o2_supply_demand_map:r44`. Override
the artifact only when intentionally testing a different image:

```bash
export O2SD_CONTAINER_IMAGE=/share/path/to/another.sif
```

The SIF SHA-256 verified on RED on 2026-07-25 is
`ee9fe0f5ab6d3fb0689b23c02d330de677afa4d0cc044df9f2eeeae42f8e5e68`.
The matching OCI index and amd64 manifest digests are recorded in
`../image_runtime_lock.tsv`.

Additional comma-separated Apptainer bind specifications can be supplied with
`O2SD_CONTAINER_BINDS`. `/share` and the repository root are bound
automatically. The Linux Rcpp compilation cache is overlaid at
`/tmp/o2sd-rcpp-cache-<uid>` so the protected repository cache is not modified.
The runtime uses `apptainer exec --cleanenv`, disables user R profiles,
redirects Python's user-package location to an empty container-local home, and
explicitly forwards Slurm, project-prefixed, and thread-count environment
variables.

Legacy `--r_module` and `R_MODULE` inputs are retained for call compatibility
but are ignored by this layer.

## Entry points

Use the same relative script as under `hpc/`. For example:

```bash
bash oxygen/code/O2_supply_demand_MAP/Docker/hpc/submit/submit_o2_fit.sh \
  --fitting_mode=invivo \
  --config_path=/absolute/path/to/O2_supply_demand.yaml \
  --out_root=/absolute/path/to/results \
  --dry_run=TRUE
```

Remove `--dry_run=TRUE` only when ready to submit. The main directories retain
their original roles:

| Folder | Responsibility |
|---|---|
| `submit/` | Fit submission and multi-warmup orchestration. |
| `array_workers/` | Slurm array task implementations. |
| `postprocess/` | Dependent postprocessing and report rerender jobs. |
| `best_fit_parameter_feature/` | Best-fit feature workflow submission. |
| `dense_grid_monotonicity_classification/` | Dense-grid monotonicity submission. |
| `fix_o2_simulation/` | Fixed-O2 simulation submission. |
| `fixo2_eigen_attractor/` | Fixed-O2 eigen-attractor tasks. |
| `joint_coupling_analysis/` | Joint coupling analysis submission. |
| `parameter_landscape/` | Parameter-landscape submission. |
| `warm_up_joint_fitting_results_extra/` | Joint warm-up result processing. |
| `bin/`, `util/` | Apptainer language shims and shared runtime. |

## One-to-one and resource checks

`hpc_script_mapping.tsv` records all 28 original files and their counterparts.
The verifier checks complete inventory coverage, shell syntax, exact Slurm
resource/dependency/log option parity, and absence of active host module
commands:

```bash
bash oxygen/code/O2_supply_demand_MAP/Docker/hpc/verify_hpc_parity.sh
```

This check is read-only and never calls `sbatch`.
