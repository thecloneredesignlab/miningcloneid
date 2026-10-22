# Figure7 and Figure8: complete current forward pipelines

This is the execution guide for the current code. It covers parameter preparation,
simulation, seed-level data generation, aggregation, validation, and plotting.
It is not limited to retained-data rendering. Chinese: [中文说明](Figure78_pipeline_ZH.md).

## Entry points and scope

- `run_Figure7.R`: instantaneous and historical Figure7, supplementary figures1–7.
- `run_Figure8.R`: Figure8 and supplementary figures1–9.
- `data_Figure7.R` / `data_Figure8.R`: the same real generation pipeline, not validator stubs.
- `draw_Figure7.R` / `draw_Figure8.R`: drawing only; no implicit simulation.

Each figure has one current configuration (`Code/config/figure7.json` or
`figure8.json`) and one data directory (`data/Figures/Figure7` or `Figure8`).
Shared modules and numerical kernels under `util/analysis` are part of those
pipelines, not alternative versions. `_ver2` publication filenames identify the
historical-lineage scientific analysis; they are not separate code/data versions.

## Upstream inputs

`Code/config/figure78_inputs.json` identifies the original final joint fits,
C01/C02 family IDs, the approved Figure7 winner seeds94/80, and the private locked
`Code/vendor/figure78_model` source closure. Preparation
requires all input files; there is no fallback to a Figure6/7/8 prediction bundle.
Original joint fitting is upstream and is not rerun by these forward pipelines.

For each seed, in vivo reads all26 `shared_invivo` parameters and `fit_config.rds`.
In vitro reads all33 `invitro_effective` parameters and that seed's
`fit_result.rds$ctx$invitro_cfg`. Every native parameter is checked on passage
through preparation and scenario construction. Only imposed oxygen and explicitly
fixed native `p_misseg` may change. `p_mis_base` and `k_o_mis` remain seed-specific.
No fitted endpoint is deduplicated; production Figure8 retains500 seeds per family.

The model closure is locked by `Code/config/manifests/figure78_model_sources.tsv`.
Byte-checked current sources are staged under `prepared/model_source`, allowing
compilation without writing into the upstream model checkout. These are current
model sources, not frozen retired Figure7/8 workers. Both local and HPC execution
require the workspace's private copy. Shared/external model overrides are refused.
`FIGURE_JOINT_RESULT_ROOT` or CLI `--fit-root` may relocate the original fits;
`--model-root` / `FIGURE_MODEL_CODE_ROOT`, if supplied, must resolve to that private copy.
Defaults are workspace-relative and portable. See [entry-point guide](Figure78_entrypoints_EN.md).

Figure8 empirical inputs are read directly from `data/share_bundle_figures_3_5`:
the two Figure3 sample/line CSVs and the original19-row `organ_lut.csv`.
The saved raw experimental adapters in C01 seed1 supply six passage lineages;
their equal-lineage geometric aggregation is recomputed, giving569799 seeded
cells and6180783 threshold cells. No old canonical-rule TSV is required.

## Run stages

From the iteration5 directory, with R and the required packages available:

```sh
Rscript Code/Figures/run_Figure7.R --stage=plan
Rscript Code/Figures/run_Figure8.R --stage=plan
```

`plan` checks source availability/model hashes and reports the real task count;
it does not simulate or submit jobs. Production has320 Figure7 tasks (three beta
conditions inside each task) and9160 Figure8 seed-block tasks. The task manifest,
not a guessed array size, is authoritative.

To generate a new full dataset in a separate directory:

```sh
Rscript Code/Figures/run_Figure7.R --stage=prepare --data-root=/absolute/new/data
Rscript Code/Figures/run_Figure7.R --stage=simulate --data-root=/absolute/new/data
Rscript Code/Figures/run_Figure7.R --stage=aggregate --data-root=/absolute/new/data
Rscript Code/Figures/run_Figure7.R --stage=validate --data-root=/absolute/new/data
Rscript Code/Figures/run_Figure7.R --stage=plot --data-root=/absolute/new/data --output-root=/absolute/new/rendered
```

Use `run_Figure8.R` for Figure8. `--stage=all` performs the same steps; it is a real
potentially expensive computation. `simulate --task-index=N` runs one manifest
task; `--task-ids=N,M` runs an explicit subset. Without a task filter it runs all
tasks sequentially. Completed compatible tasks are reused. Figure8 checkpoints
after each seed inside a block. Aggregation refuses missing tasks.

Existing published/retained data are protected. To deliberately replace a current
data directory, `prepare --rebuild=yes` first moves it to a recoverable audit
backup; it does not delete it or create a new version-numbered active directory.
Do not use this option merely to inspect or redraw existing results.

## Scientific steps preserved

Figure7 creates its own total trajectories, then replays passive history and
non-inherited direct-origin labels against those same trajectories and passage
samples. Instantaneous A is `100 F_HL/(F_HL+F_LL)`; instantaneous B is
`100 sum(J_LL D_L)/sum(J_LL X_L)`. D contains only currently surviving direct H→L
daughters before their first division. Their daughters do not inherit D. History A
is the current low-cell stock fraction of the permanently tagged lineage; history B
is that lineage's current gross low-cell birth contribution. No growth share
subtracts death. Figure7 intentionally uses pure2N/4N states in both contexts.

Figure8 uses finite-time normalized states, not dominant-eigenvector predictions.
In vivo uses stable finite matrix exponentials. In vitro uses each seed's fitted
Gaussian cohorts, its own DT/density configuration, and the threshold-triggered
integer-day stochastic passage protocol. The live-state Euler step, randomized
integerization, hypergeometric without-replacement sample, and separate label RNG
are preserved. Passage is not a fixed5-day dilution and is not continuous growth.

Culture repeats are averaged within each seed first. Seed-level values are then
pooled500+500 or grouped500 per family, with type-7 median/Q25/Q75. Mixed cohorts
are combined as cell states/counts, not endpoint ploidy averages. All ten time
snapshots are taken from the same continuing1000-day trajectory. A proven
no-threshold-crossing bound may reduce200 identical repeats to one Euler run;
both executed and nominal repeat counts are retained.

## Data layout and provenance

```text
Figure7/ or Figure8/
  prepared/     full seed/context inputs, task manifest, bundle, current model closure
  tasks/        current task outputs and compatible resumable checkpoints
  seed_data/    Figure8 per-seed endpoints and time snapshots
  tables/       all main/supplementary plotting tables
  provenance/   original-input hashes, parameter audit, current-code hashes,
                raw passage derivation, runtime, numerical gates, plot-input hashes
  .cache/       generated compilation cache, not scientific prediction input
```

New result provenance is generated from current inputs and code. It is not re-signed
as if produced by an old worker. A changed numerical code/config/input fingerprint
rejects resume. Plotting verifies the new full-dataset receipt without needing any
retired worker or archive. The separate retained-data compatibility path remains
available for the unchanged published results. All26 main/supplementary product
entries use the same accepted plot builders and publication filenames.

## Representative cold-start acceptance

```sh
Rscript Code/Figures/validate_Figure78_cold_start.R --data-root=/absolute/empty/test/data
Rscript Code/Figures/validate_Figure78_pipeline.R --data-root=/absolute/empty/test/data
```

The first command uses exactly the production workers and1000-day horizon but two
passage repeats and, for Figure8, two seeds per family. R-level source/read guards
reject archives, frozen retired workers, and retained Figure7/8 scientific results.
Before the guard starts, published files are byte-hashed only for protection;
no reference values enter the calculation. The cold-start entry writes its own
protection baseline, so the two commands do not require a manually prepared audit
fixture. Compiler caches may be reused; no previous scientific bundle/task/summary
may be present. The second command performs independent metric, passage, quantile, and
retained-reference regression checks after the guarded calculation has completed.
Representative outputs are explicitly refused as500+500 publication data.

## HPC boundary

`Code/hpc/run_Figure78_task.sh` runs the same current worker using a prepared task
manifest and `SLURM_ARRAY_TASK_ID`. It automatically loads the scoped paths and
checks the current preparation/task manifest. It does not submit an array. A later approved
submission must use `xxlarge`; its index range must cover the manifest, without a
`%N` concurrency cap. Full runs, cluster environment setup, and Slurm execution are
not performed or claimed by local representative acceptance.
