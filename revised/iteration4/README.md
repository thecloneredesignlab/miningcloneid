# Reproducible figure and manuscript workspace

This directory is a relocatable, self-bootstrapping figure and manuscript
workspace. Its name and depth inside the repository are not part of the
runtime contract.

The minimum portable bundle contains only:

```text
manager.sh
Code/
```

Read the command-line usage from any current working directory:

```bash
bash /absolute/path/to/workspace/manager.sh --help
```

Run the complete workflow with explicit scientific-input roots:

```bash
bash /absolute/path/to/workspace/manager.sh \
  --invitro-result-root=/path/to/fit_invitro_unified_500seed_r442_exact_20260825_032031 \
  --invivo-result-root=/path/to/fit_invivo_unified_500seed_r442_exact_20260825_032031 \
  --joint-result-root=/path/to/fit_joint_invivo_clusters_global_invitro_best_500seed_r442_exact_20260826_033633 \
  --gemcitabine-data-root=/path/to/data/InVivoData_Gemcitabine \
  --ltee-data-root=/path/to/data/InVitroData_LTEE
```

Container wrappers are stored in `Code/Docker/`. Each wrapper accepts the same
scientific-input and optional control arguments as `manager.sh`. Display their
complete embedded usage with either the standard `--help` spelling or the
supported `--hlep` compatibility alias:

```bash
bash Code/Docker/manager_docker.sh --help
bash Code/Docker/manager_hpc.sh --help
```

The local Docker wrapper additionally requires `--docker-image=IMAGE` and
enforces the approved `o2_supply_demand_map:r44` image ID. The HPC wrapper
accepts `--sif-image=PATH` and defaults to
`/share/lab_crd/lab_crd/taoli/Docker/o2_supply_demand_map_r44.sif`. The HPC
wrapper directly invokes Apptainer or Singularity; it does not submit a Slurm
job or choose scheduler resources.

## Current local Docker full-recomputation guide

The command below is specific to the current local checkout. It forces all
three persistent scientific-analysis cache families to be rebuilt:

- fixed-O2 analysis;
- Figure 4 in-vivo t-SNE landscape;
- Figure 6 objective-filtered multi-seed response ensemble.

It also configures headless R graphics with `bitmapType = "cairo"`, mounts the
current Git worktree metadata read-only, mounts all external scientific inputs
read-only, and leaves the figure workspace writable for generated outputs.

```bash
docker run --rm --init \
  --platform linux/amd64 \
  --user "$(id -u):$(id -g)" \
  --network none \
  --workdir /Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration4 \
  --env HOME=/tmp \
  --env TMPDIR=/tmp \
  --env XDG_CACHE_HOME=/tmp/.cache \
  --mount "type=bind,source=/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures,target=/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures" \
  --mount "type=bind,source=/Users/4482173/Documents/GitHub/miningcloneid/.git,target=/Users/4482173/Documents/GitHub/miningcloneid/.git,readonly" \
  --mount "type=bind,source=/Users/4482173/Documents/GitHub/soft_couping_org/oxygen/code/O2_supply_demand_MAP,target=/Users/4482173/Documents/GitHub/soft_couping_org/oxygen/code/O2_supply_demand_MAP,readonly" \
  --mount "type=bind,source=/Users/4482173/Documents/GitHub/soft_couping_org/oxygen/results/fit_invitro_unified_500seed_r442_exact_20260825_032031,target=/Users/4482173/Documents/GitHub/soft_couping_org/oxygen/results/fit_invitro_unified_500seed_r442_exact_20260825_032031,readonly" \
  --mount "type=bind,source=/Users/4482173/Documents/GitHub/soft_couping_org/oxygen/results/fit_invivo_unified_500seed_r442_exact_20260825_032031,target=/Users/4482173/Documents/GitHub/soft_couping_org/oxygen/results/fit_invivo_unified_500seed_r442_exact_20260825_032031,readonly" \
  --mount "type=bind,source=/Users/4482173/Documents/GitHub/soft_couping_org/oxygen/results/fit_joint_invivo_clusters_global_invitro_best_500seed_r442_exact_20260826_033633,target=/Users/4482173/Documents/GitHub/soft_couping_org/oxygen/results/fit_joint_invivo_clusters_global_invitro_best_500seed_r442_exact_20260826_033633,readonly" \
  --mount "type=bind,source=/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/data/InVivoData_Gemcitabine,target=/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/data/InVivoData_Gemcitabine,readonly" \
  o2_supply_demand_map:r44 \
  bash -lc '
    set -euo pipefail
    printf "%s\n" "options(bitmapType = \"cairo\")" > /tmp/figure-workspace-Rprofile
    export R_PROFILE_USER=/tmp/figure-workspace-Rprofile
    exec bash /Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration4/manager.sh \
      --invitro-result-root=/Users/4482173/Documents/GitHub/soft_couping_org/oxygen/results/fit_invitro_unified_500seed_r442_exact_20260825_032031 \
      --invivo-result-root=/Users/4482173/Documents/GitHub/soft_couping_org/oxygen/results/fit_invivo_unified_500seed_r442_exact_20260825_032031 \
      --joint-result-root=/Users/4482173/Documents/GitHub/soft_couping_org/oxygen/results/fit_joint_invivo_clusters_global_invitro_best_500seed_r442_exact_20260826_033633 \
      --gemcitabine-data-root=/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/data/InVivoData_Gemcitabine \
      --ltee-data-root=/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/data/InVitroData_LTEE \
      --n-core=8 \
      --recompute-fixed-o2 \
      --recompute-invivo-tsne \
      --rebuild-figure6-grid
  '
```

If `manuscript/` is absent, the manager copies the repository-root
`manuscript/` tree before analysis. Existing workspace manuscripts are
preserved and are never silently replaced.

Optional controls are `--n-core=N`, `--recompute-fixed-o2`,
`--recompute-invivo-tsne`, and `--rebuild-figure6-grid`. Every run executes
all per-figure analysis and drawing entries, publishes figures, validates the
scientific inputs before analysis and after report generation, hashes
generated intermediates, verifies published-figure MD5 identity, and
rebuilds the model-dependent manuscript tables before generating the embedded
HTML report. The manager does not compile the TeX
manuscript or create, replace, or remove its PDF.

## Targeted Figure 5D DE-initial/optimizer rerun

The active Figure 5D compares the actual differential-evolution (DE) initial
population with the empirical distribution of optimizer endpoints. It does not
use generalized-posterior, MCMC, Welsch-weighted reference, or HPC outputs.
C01--C06 each contribute one selected pair from the six primary warm-start
regions; no secondary-cluster layer is used. The common *in vitro* anchor is
seed228.

Regenerate and audit the active products locally with:

```bash
Rscript Code/Figures/prepare_Figure5F_de_initial_population.R
Rscript Code/Figures/audit_Figure5F_prior_optimizer_inputs.R
Rscript Code/Figures/build_Figure5F_prior_optimizer_products.R
Rscript Code/Figures/build_Figure5F_supplementary_table.R
Rscript Code/Figures/draw_Figure5.R \
  --invitro-result-root=/path/to/fit_invitro_unified_500seed_r442_exact_20260825_032031 \
  --invivo-result-root=/path/to/fit_invivo_unified_500seed_r442_exact_20260825_032031 \
  --joint-result-root=/path/to/fit_joint_invivo_clusters_global_invitro_best_500seed_r442_exact_20260826_033633 \
  --gemcitabine-data-root=/path/to/data/InVivoData_Gemcitabine \
  --ltee-data-root=/path/to/data/InVitroData_LTEE
```

The DE initialization is replayed for joint seeds 1--500 with 400 members per
run, yielding 200,000 paired starting members per family and parameter. The
active panel shows the two reconstructed natural-scale marginals separately:
the upper blue half is *in vivo* and the lower pink half is *in vitro*.
Pair-specific DE initial populations have gray dashed outlines, whereas the
500 optimizer endpoints have C01--C06 solid outlines. Density peaks are
normalized separately within each displayed distribution. The 5th--95th
percentile spans are descriptive numerical-search summaries—not posterior or
confidence intervals, biological replicates, or proof of structural
identifiability. Paired `log2(in vivo / in vitro)` products are retained only
as a secondary directional/search audit and are not the Figure 5D coordinate.
The full data and interpretation contract is in
`data/Figures/Figure5/figure5f_chart_contract.md`.

## Targeted Figure 6 supplementary rerun

Figure 6 and its parent-indexed Supplementary Figures 6-1 through 6-4 use the
analytical inputs packaged under this iteration4 workspace and load every
model calculation from the required read-only implementation at
`/Users/4482173/Documents/GitHub/soft_couping_org/oxygen/code/O2_supply_demand_MAP`.
They do not use the model copy in the HypoxiaLTEEFigures checkout. The exact
external model-source paths and MD5 values used in a run are
written to `data/Figures/Figure6/invitro_model_code_provenance.tsv`; this MD5
record, rather than a user-reported reference commit alone, is the
reproduction identity.

After endpoint caches are available, regenerate the analysis and figures with:

```bash
cd /absolute/path/to/workspace
Rscript Code/Figures/data_Figure6.R --n-core=8 --rebuild=FALSE --n-resample=100
Rscript Code/Figures/draw_Figure6.R
Rscript Code/Figures/data_Supp_Figure6_1.R --n-core=8 --rebuild=FALSE
Rscript Code/Figures/draw_Supp_Figure6_1.R
Rscript Code/Figures/data_Supp_Figure6_2.R
Rscript Code/Figures/draw_Supp_Figure6_2.R
Rscript Code/Figures/data_Supp_Figure6_3.R --n-core=8 --rebuild=FALSE
Rscript Code/Figures/draw_Supp_Figure6_3.R
Rscript Code/Figures/data_Supp_Figure6_4.R --n-core=8 --rebuild=FALSE
Rscript Code/Figures/draw_Supp_Figure6_4.R
Rscript Code/Figures/build_figure6_manuscript_tables.R
```

On macOS, long *in vitro* endpoint grids can be precomputed with independent R
processes to avoid inheriting an initialized OpenMP runtime through `fork`.
Four or eight workers may be launched with distinct `--worker-id` values and
the same `--worker-count`, first for `--stage=q20` and then for
`--stage=dense`. Workers receive disjoint manifest rows and write the same
atomic caches consumed by `data_Figure6.R`; process isolation changes only
task scheduling. For example, one of four workers is:

```bash
Rscript Code/Figures/compute_Figure6_invitro_endpoint.R \
  --stage=q20 --worker-id=1 --worker-count=4 --rebuild=false
```

The analysis ranks the 500 endpoints separately within each of six joint
warm-start pairs, uses the lowest 10% as the primary equal-weight ensemble,
and evaluates nested 5% and 20% sensitivity sets. It also reports exact
14-parameter endpoint multiplicity and repeats the qualitative-claim audit
after weighting each unique parameter endpoint once. Per-seed and
unique-endpoint caches make the calculation restartable. The joint *in vivo*
and *in vitro* panels use the two parameter vectors from each same retained
joint endpoint; the separate-fit *in vitro* best seed is not substituted into
the joint analysis.

## RED HPC Figure 6 data-only rebuild

`Code/hpc/run_figure6_data_hpc.sh` is the data-only entry point for the
allocated RED node `hpctpa3pc0028`. It uses 63 independent R workers inside
`/share/lab_crd/taoli/Docker/o2_supply_demand_map_r442_hpc_exact.sif`, loads
the external model from
`/share/lab_crd/taoli/Project/soft_couping_org/oxygen/code/O2_supply_demand_MAP`,
and reads the three current fit-result trees under the corresponding HPC
`oxygen/results` directory. The runner archives only prior Figure 6 data
products and then rebuilds the main and Supplementary Figure 6-1 through 6-4
data with fresh caches. It never calls a `draw_*.R` script or manuscript
renderer.

The current fresh Figure 3/4/5 staged dependencies listed by the runner must
be synchronized into the HPC iteration4 workspace before launch. Validate the
complete environment on the allocated node with:

```bash
bash Code/hpc/run_figure6_data_hpc.sh --preflight-only
```

Then start the full data rebuild from the existing Slurm allocation. Drawing
and manuscript-table generation are performed locally only after the validated
HPC data directories have been copied back.

Directory contract:

- `Code/Figures/data_FigureX.R` and `data_Supp_FigureX_Y.R`: analysis entry points.
- `Code/Figures/draw_FigureX.R` and `draw_Supp_FigureX_Y.R`: drawing entry points.
- `Code/Figures/util`: reusable runtime, analysis, graphics, model, and
  workflow functions organized by function.
- `Code/config`: bootstrap configuration and static analysis parameters.
- `data/Figures/FigureX` and `data/Figures/Supp_FigureX_Y`: all staged scientific
  inputs, runtime configurations, and generated intermediates.
- `Figures`: final assembled figure outputs only.
- `manuscript/Figures`: published copies referenced by the TeX manuscript.
- `audit/md5/scientific_input_md5.tsv`: observed input sizes and MD5 values,
  compared with the fixed baseline.
- `audit/logs`: complete manager logs.
- `audit/parameters`: resolved pipeline parameters.
- `audit/manifests`: allowed input roots and figure entry points.
- `audit/reports`: validation results.

Scientific inputs are restricted to the roots published in
`audit/manifests/allowed_scientific_inputs.txt`. Model calculations must call
the exact external `soft_couping_org/oxygen/code/O2_supply_demand_MAP` root
shown above; the local vendor/model copies are not calculation inputs. The
resolved paths are recorded in
`audit/parameters/pipeline_parameters.tsv`.

The immutable validation baseline is
`Code/config/manifests/expected_scientific_input_md5.tsv`. It identifies each
input by a functional root ID and a relative path, so it does not depend on the
workspace name or location. `audit/reports/input_md5_validation.tsv` records
the latest validation summary. No MD5 validation is performed for copied
inputs, generated intermediates, code, figures, manuscript bootstrap files, or
published figure copies.
