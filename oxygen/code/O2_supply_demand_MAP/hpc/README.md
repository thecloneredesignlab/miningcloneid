# RED HPC Usage Guide (Apptainer SIF)

This document describes how to run the O2 supply-demand model on the Moffitt
RED HPC. Slurm orchestration runs on the host, but **all R and Python commands
must run through the Apptainer SIF image**.

## 1. Fixed paths and execution architecture

| Item | Current convention |
|---|---|
| Login | `4482173@red.moffitt.org` |
| Established project checkout | `/share/lab_crd/lab_crd/taoli/Project/miningcloneid_soft_coupling` |
| Default SIF | `/share/lab_crd/lab_crd/taoli/Docker/o2_supply_demand_map_r44.sif` |
| Recorded SIF SHA-256 | `08d52a9661d8ab21d083a488c994bebc4144c40fd38c11aeb1acb948abe237ad` |
| SIF-backed entry points | `oxygen/code/O2_supply_demand_MAP/Docker/hpc/` |
| Host-runtime reference layer | `oxygen/code/O2_supply_demand_MAP/hpc/` |

The execution chain is:

```text
local machine -> SSH login node -> host sbatch -> Slurm compute node
                                           -> Docker/hpc runtime wrapper
                                           -> apptainer exec <SIF>
                                           -> R/Python inside the SIF
```

Core rules:

- `sbatch`, `squeue`, `sacct`, and dependency orchestration remain on the HPC
  host.
- Do not run `sbatch` from inside the SIF.
- When SIF execution is required, submit the corresponding script under
  `Docker/hpc/`. Do not use the host-R version under `hpc/`.
- `Docker/hpc/` mirrors the Slurm resources, arrays, dependencies, log paths,
  workflow entry points, and command-line interfaces under `hpc/`. Only the
  R/Python runtime changes.
- The runtime uses `apptainer exec --cleanenv`, disables user `.Rprofile`,
  `.Renviron`, and Python user packages, and avoids contamination from personal
  environments.

## 2. Log in to RED

Interactive login:

```bash
ssh 4482173@red.moffitt.org
```

When running one remote command from a local machine, explicitly start a login
shell. A non-login shell may not expose `sbatch`, `squeue`, `sacct`, or other
Slurm commands:

```bash
ssh 4482173@red.moffitt.org \
  "bash -lc 'squeue -u 4482173'"
```

If already logged in but Slurm commands are unavailable, start a fresh login
shell:

```bash
exec bash -l
```

## 3. Verify the checkout, branch, and commit

Enter the checkout assigned to the run. The established `soft_coupling`
checkout is:

```bash
HPC_PROJECT_ROOT=/share/lab_crd/lab_crd/taoli/Project/miningcloneid_soft_coupling
cd "${HPC_PROJECT_ROOT}"

git status --short --branch
git branch --show-current
git rev-parse HEAD
git remote -v
```

Only synchronize after confirming that the branch is correct and that no
existing work will be overwritten:

```bash
git pull --ff-only
git rev-parse HEAD
```

Use a separate checkout or worktree for a new branch or isolated experiment.
Do not switch a shared checkout that may be used by another active job. Record
the following before submission:

- absolute HPC project path;
- branch name;
- full `git rev-parse HEAD` commit;
- `git status --short` output;
- result root for this run.

## 4. SIF preflight

### 4.1 Verify Apptainer, image readability, and checksum

```bash
SIF_IMAGE=/share/lab_crd/lab_crd/taoli/Docker/o2_supply_demand_map_r44.sif
EXPECTED_SIF_SHA256=08d52a9661d8ab21d083a488c994bebc4144c40fd38c11aeb1acb948abe237ad

command -v apptainer
test -r "${SIF_IMAGE}"
printf '%s  %s\n' "${EXPECTED_SIF_SHA256}" "${SIF_IMAGE}" | sha256sum -c -
```

Stop if the checksum does not match. Confirm whether the image was
intentionally replaced; do not silently fall back to host R.

The checksum above was verified on RED on 2026-07-30. This SIF predates the
current target-lock refresh for `svglite 2.2.2` and `textshaping 1.0.5`. If a
run must validate the refreshed lock, rebuild and republish the image, then
verify and record the replacement SIF before using it. Do not claim that the
recorded older SIF matches the refreshed lock.

### 4.2 Verify host/SIF script parity

```bash
cd "${HPC_PROJECT_ROOT}"
bash oxygen/code/O2_supply_demand_MAP/Docker/hpc/verify_hpc_parity.sh
```

This command checks inventory coverage, shell syntax, Slurm resource parity,
and host-module leakage. It does not call `sbatch`. The expected final line is:

```text
HPC parity verification passed: 28 mapped files; Slurm resources unchanged.
```

### 4.3 Verify R and Python inside the SIF

```bash
(
  export PROJECT_ROOT="${HPC_PROJECT_ROOT}"
  export O2SD_CONTAINER_IMAGE="${SIF_IMAGE}"
  source "${HPC_PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/Docker/hpc/util/o2_supply_demand_map_apptainer_runtime.sh"

  o2sd_container_prepare
  command -v Rscript
  command -v python3
  Rscript -e 'cat(R.version.string, "\n"); print(.libPaths())'
  python3 --version
)
```

`Rscript` and `python3` should resolve under `Docker/hpc/bin/`. Those shims
invoke the interpreters inside the SIF. SIF execution does not require a host R
module; the legacy `--r_module` option is retained for compatibility but is
ignored by this layer.

## 5. Use `xxlarge` correctly

`xxlarge` is a Slurm QoS, not an automatic CPU, memory, or GPU allocation.
Every submission must still specify or verify:

- walltime;
- CPUs per task;
- memory per task;
- array size and concurrency;
- whether a GPU is required. The model-fitting workflow is CPU-only by default.

The current unified submitter defaults are:

| Mode | QoS | Walltime | CPUs/task | Memory/task |
|---|---:|---:|---:|---:|
| in-vivo | `xlarge` | `12:00:00` | 22 | 16G |
| in-vitro | `xxlarge` | `12:00:00` | 22 | 16G |
| joint | `xxlarge` | `12:00:00` | 22 | 16G |

This repository treats 12 hours as the `xxlarge` maximum. Cluster policy may
change, so query the live definition before a production submission:

```bash
sacctmgr show qos xxlarge \
  format=Name,MaxWall,MaxTRES,MaxJobsPU,MaxSubmitPU -P
```

Do not use `xxlarge` as a partition name, and do not omit CPU, memory, or time
requests merely because `xxlarge` is selected. The unified submitter uses:

```text
--invitro_qos=xxlarge
--invitro_time_limit=12:00:00
--invitro_n_cores=22
--invitro_mem=16G
```

## 6. Use the unified SIF-backed submitter

The primary production entry point is:

```text
oxygen/code/O2_supply_demand_MAP/Docker/hpc/submit/submit_o2_fit.sh
```

Inspect all options supported by the exact checkout being submitted:

```bash
bash "${HPC_PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/Docker/hpc/submit/submit_o2_fit.sh" \
  --help
```

### 6.1 Run a dry-run first

The example below validates a one-seed in-vitro submission. `--dry_run=TRUE`
prints the exact `sbatch` command without submitting it:

```bash
RUN_TAG="fit_invitro_sif_smoke_$(date +%Y%m%d_%H%M%S)"
RESULT_ROOT="${HPC_PROJECT_ROOT}/oxygen/results"
DEATH_DATA="${HPC_PROJECT_ROOT}/data/InVitroData/sum159_dead_cell_endpoint_likelihood_ready_20260731.tsv"
export O2SD_CONTAINER_IMAGE="${SIF_IMAGE}"

bash "${HPC_PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/Docker/hpc/submit/submit_o2_fit.sh" \
  --fitting_mode=invitro \
  --project_root="${HPC_PROJECT_ROOT}" \
  --config_path="${HPC_PROJECT_ROOT}/oxygen/config/O2_supply_demand.yaml" \
  --out_root="${RESULT_ROOT}" \
  --invitro_run_prefix="${RUN_TAG}" \
  --invitro_total_seeds=1 \
  --invitro_array_tasks=1 \
  --invitro_seeds_per_task=1 \
  --invitro_qos=xxlarge \
  --invitro_time_limit=12:00:00 \
  --invitro_n_cores=22 \
  --invitro_mem=16G \
  --death_data_path="${DEATH_DATA}" \
  --dry_run=TRUE
```

Review the dry-run output and confirm that:

- every path points to the intended HPC checkout or `/share`;
- the worker path contains `/Docker/hpc/`;
- QoS, walltime, CPUs, memory, and array size are correct;
- config, parameter table, fit objects, and data files exist;
- result and log directories are not used by another production run;
- no unintended GPU resource is requested.

### 6.2 Submit a smoke test

After approving the dry-run, repeat the same command with the last option set
to:

```text
--dry_run=FALSE
```

The submitter prints the Slurm job ID and records the submission command,
resources, and script paths in `run_provenance.tsv` under the run directory.
Save the job ID and immediately verify the submitted resources with `squeue`
or `scontrol`.

### 6.3 Production-size example

Increase scale only after the smoke test completes successfully and its output
is validated. For 500 seeds with one seed per array element, use:

```text
--invitro_total_seeds=500
--invitro_array_tasks=500
--invitro_seeds_per_task=1
```

The seed plan must satisfy:

```text
array_tasks * seeds_per_task = total_seeds
```

The current single-side in-vivo/in-vitro submitter does not use the
multi-warmup `--array_max_concurrent` option. Do not assume that this option
limits ordinary fitting arrays. If a concurrency cap is needed, inspect the
actual `sbatch --array` command generated by the unified submitter and confirm
the current QoS policy before submission.

### 6.4 Supported fitting modes

| Purpose | Key options |
|---|---|
| in-vivo only | `--fitting_mode=invivo` |
| in-vitro only | `--fitting_mode=invitro` |
| direct joint fit | `--fitting_mode=joint --joint_fitting_mode=DIRECT` |
| single-side fits followed by joint | `--fitting_mode=joint --joint_fitting_mode=JOINT` |
| multi-warmup joint fit | `--fitting_mode=joint --joint_fitting_mode=MULTI_WARMUP` |

A JOINT or MULTI_WARMUP run may submit preparation, selection, array,
collection, and reporting jobs connected by dependencies. Save every job ID
and generated task manifest, not only the final job ID.

## 7. Put a custom Slurm worker inside the SIF

If a new workflow does not yet have a `Docker/hpc/` wrapper, use the following
structure. Submit `sbatch` on the host; inside the worker, source the shared
runtime before invoking `Rscript`:

```bash
#!/usr/bin/env bash
#SBATCH --job-name=o2_sif_job
#SBATCH --qos=xxlarge
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/share/lab_crd/lab_crd/taoli/Project/miningcloneid_soft_coupling/oxygen/results/log/o2_sif_%j.out
#SBATCH --error=/share/lab_crd/lab_crd/taoli/Project/miningcloneid_soft_coupling/oxygen/results/log/o2_sif_%j.err

set -euo pipefail

PROJECT_ROOT=/share/lab_crd/lab_crd/taoli/Project/miningcloneid_soft_coupling
export PROJECT_ROOT
export O2SD_CONTAINER_IMAGE=/share/lab_crd/lab_crd/taoli/Docker/o2_supply_demand_map_r44.sif

source "${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/Docker/hpc/util/o2_supply_demand_map_apptainer_runtime.sh"
o2sd_container_ignore_r_module ""

Rscript "${PROJECT_ROOT}/absolute/or/repository/relative/runner.R" \
  --out_dir="${PROJECT_ROOT}/oxygen/results/example"
```

Submit the worker from the host:

```bash
sbatch /absolute/path/to/worker.sub
```

A permanent repository worker must be added to both the `hpc/` and
`Docker/hpc/` one-to-one layers and must pass `verify_hpc_parity.sh`. Do not
maintain an unverified container-only copy.

## 8. Binds and filesystem paths

The runtime automatically binds:

- `/share:/share`;
- `PROJECT_ROOT` at the same absolute path inside the container;
- a temporary Rcpp cache at the model's Rcpp cache path.

Inputs or outputs outside `/share` and the project checkout require an explicit
bind. Separate multiple bind specifications with commas:

```bash
export O2SD_CONTAINER_BINDS="/external/input:/external/input:ro,/external/output:/external/output"
```

Then run the `Docker/hpc/` submitter. Do not bind an entire user home, SSH
directory, credential directory, or another unnecessary sensitive path.

## 9. Monitor, diagnose, and cancel jobs

### 9.1 Pending or running jobs

```bash
squeue -u 4482173
squeue -j JOB_ID
scontrol show job JOB_ID
```

Inspect a specific array element with:

```bash
scontrol show job JOB_ID_TASK_ID
```

### 9.2 Finished jobs

After a job leaves `squeue`, inspect its terminal state with `sacct`:

```bash
sacct -j JOB_ID \
  --format=JobID,JobName,State,ExitCode,Elapsed,Timelimit,AllocCPUS,ReqMem,MaxRSS
```

Absence from `squeue` does not mean success. Treat a run as successful only
when all of the following are true:

1. the relevant `sacct` entries are `COMPLETED` with `ExitCode=0:0`;
2. the `.err` log contains no unhandled error;
3. the `.out` tail matches the expected workflow completion;
4. expected result files exist, are readable, and are complete;
5. the commit, scripts, resources, and command in `run_provenance.tsv` match the
   intended run;
6. the SIF path and SHA-256 were recorded and verified separately.

Inspect logs with:

```bash
tail -n 100 /absolute/path/to/job.out
tail -n 100 /absolute/path/to/job.err
```

### 9.3 Cancel a job

```bash
scancel JOB_ID
```

Cancel only the explicitly named job or jobs submitted for the current task.
Do not cancel unrelated jobs owned by the same account.

## 10. Common problems

### `sbatch: command not found`

The command is running in a non-login shell. Log in again, or use `bash -lc` or
`exec bash -l`.

### `Container SIF is not readable`

Check `O2SD_CONTAINER_IMAGE`, file permissions, and `/share` availability. Do
not fall back automatically to host R.

### SIF checksum mismatch

Stop the submission and determine whether the image was intentionally updated.
A replacement image must have a recorded path, SHA-256, platform, and
verification date.

### The log shows a host R module being loaded

The submission likely used `hpc/` instead of `Docker/hpc/`. Check the absolute
submitter and worker paths. The SIF layer retains `R_MODULE` only for interface
compatibility and must not load host R.

### An input file is missing inside the container

Prefer inputs under the project checkout or `/share`. Bind any other path with
`O2SD_CONTAINER_BINDS`, and verify that the compute node can access it.

### `PENDING (ReqNodeNotAvail)`

Do not assume that an invalid node name was requested. Inspect
`scontrol show job`, requested walltime, QoS limits, partition state, and
maintenance reservations. A long walltime may make the job overlap a
maintenance window.

### The job is absent from `squeue`

Use `sacct` to inspect the terminal state, then verify logs and result files.
The job may have completed, failed, timed out, or been cancelled.

## 11. Final pre-submission checklist

- [ ] The session is logged in to `4482173@red.moffitt.org`, and Slurm commands
      come from a login shell.
- [ ] The run uses the assigned isolated HPC checkout, correct branch, and
      correct commit.
- [ ] `git status --short` was reviewed and no other task's work is overwritten.
- [ ] The submitter or worker comes from `Docker/hpc/`.
- [ ] The SIF is readable and its SHA-256 matches the approved artifact.
- [ ] `verify_hpc_parity.sh` passes.
- [ ] The container R/Python sanity check passes.
- [ ] `--dry_run=TRUE` was run and the generated `sbatch` command was reviewed.
- [ ] QoS, walltime, CPUs, memory, array size, and concurrency match the plan.
- [ ] All input, output, and log paths are absolute and visible to compute nodes.
- [ ] The result directory uses a new, traceable run prefix.
- [ ] All job IDs, dependencies, commit, SIF SHA, and result roots are recorded.
- [ ] Final acceptance uses `sacct`, logs, and result artifacts together.

## 12. HPC directory responsibilities

| Folder | Responsibility |
|---|---|
| `submit/` | Unified fitting submission and multi-warmup orchestration. |
| `array_workers/` | In-vivo, in-vitro, joint, fixed-O2, landscape, and warm-up array workers. |
| `postprocess/` | Dependent postprocessing and report rendering. |
| `best_fit_parameter_feature/` | Best-fit feature workflow submission. |
| `dense_grid_monotonicity_classification/` | Dense-grid monotonicity submission. |
| `fix_o2_simulation/` | Fixed-O2 simulation array submission. |
| `fixo2_eigen_attractor/` | Fixed-O2 eigen-attractor tasks. |
| `joint_coupling_analysis/` | Joint coupling analysis submission. |
| `parameter_landscape/` | Parameter-landscape workflow submission. |
| `warm_up_joint_fitting_results_extra/` | Joint warm-up result processing. |
| `util/` | Shared HPC shell and provenance helpers. |

Additional image implementation and lock information is available in:

- `oxygen/code/O2_supply_demand_MAP/Docker/hpc/README.md`
- `oxygen/code/O2_supply_demand_MAP/Docker/image_runtime_lock.tsv`
- `oxygen/code/O2_supply_demand_MAP/Docker/hpc/util/o2_supply_demand_map_apptainer_runtime.sh`
