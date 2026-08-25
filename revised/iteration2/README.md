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
  --invitro-result-root=/path/to/fit_invitro_O2_buffering_500seed \
  --invivo-result-root=/path/to/fit_invivo_O2_buffering_500seed \
  --joint-result-root=/path/to/fit_joint_multi_warmup_results \
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
  --workdir /Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration2 \
  --env HOME=/tmp \
  --env TMPDIR=/tmp \
  --env XDG_CACHE_HOME=/tmp/.cache \
  --mount "type=bind,source=/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures,target=/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures" \
  --mount "type=bind,source=/Users/4482173/Documents/GitHub/miningcloneid/.git,target=/Users/4482173/Documents/GitHub/miningcloneid/.git,readonly" \
  --mount "type=bind,source=/Users/4482173/Documents/GitHub/soft_coupling/oxygen/results/fit_invitro_O2_buffering_500seed,target=/Users/4482173/Documents/GitHub/soft_coupling/oxygen/results/fit_invitro_O2_buffering_500seed,readonly" \
  --mount "type=bind,source=/Users/4482173/Documents/GitHub/soft_coupling/oxygen/results/fit_invivo_O2_buffering_500seed,target=/Users/4482173/Documents/GitHub/soft_coupling/oxygen/results/fit_invivo_O2_buffering_500seed,readonly" \
  --mount "type=bind,source=/Users/4482173/Documents/GitHub/soft_coupling/oxygen/results/fit_joint_multi_warmup_tsne_sigmaN0p1216_objmin_500seed_20260714_021540,target=/Users/4482173/Documents/GitHub/soft_coupling/oxygen/results/fit_joint_multi_warmup_tsne_sigmaN0p1216_objmin_500seed_20260714_021540,readonly" \
  --mount "type=bind,source=/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/data/InVivoData_Gemcitabine,target=/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/data/InVivoData_Gemcitabine,readonly" \
  o2_supply_demand_map:r44 \
  bash -lc '
    set -euo pipefail
    printf "%s\n" "options(bitmapType = \"cairo\")" > /tmp/figure-workspace-Rprofile
    export R_PROFILE_USER=/tmp/figure-workspace-Rprofile
    exec bash /Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration2/manager.sh \
      --invitro-result-root=/Users/4482173/Documents/GitHub/soft_coupling/oxygen/results/fit_invitro_O2_buffering_500seed \
      --invivo-result-root=/Users/4482173/Documents/GitHub/soft_coupling/oxygen/results/fit_invivo_O2_buffering_500seed \
      --joint-result-root=/Users/4482173/Documents/GitHub/soft_coupling/oxygen/results/fit_joint_multi_warmup_tsne_sigmaN0p1216_objmin_500seed_20260714_021540 \
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
generates the embedded HTML report. The manager does not compile the TeX
manuscript or create, replace, or remove its PDF.

## Targeted Figure 5D DE-initial/optimizer rerun

The active Figure 5D compares the actual differential-evolution (DE) initial
population with the empirical distribution of optimizer endpoints. It does not
use generalized-posterior, MCMC, Welsch-weighted reference, or HPC outputs.
C01, C02, and C03 each contribute one selected pair: within each family, the
retained pair is the one whose *in vivo* member has the lower original
separate-*in vivo* objective. The common *in vitro* anchor is seed10.

Regenerate and audit the active products locally with:

```bash
Rscript Code/Figures/prepare_Figure5F_de_initial_population.R
Rscript Code/Figures/audit_Figure5F_prior_optimizer_inputs.R
Rscript Code/Figures/build_Figure5F_prior_optimizer_products.R
Rscript Code/Figures/build_Figure5F_supplementary_table.R
Rscript Code/Figures/draw_Figure5.R \
  --invitro-result-root=/path/to/fit_invitro_O2_buffering_500seed \
  --invivo-result-root=/path/to/fit_invivo_O2_buffering_500seed \
  --joint-result-root=/path/to/fit_joint_multi_warmup_results \
  --gemcitabine-data-root=/path/to/data/InVivoData_Gemcitabine \
  --ltee-data-root=/path/to/data/InVitroData_LTEE
```

The DE initialization is replayed for joint seeds 1--500 with 400 members per
run, yielding 200,000 paired starting members per family and parameter. The
active panel shows the two reconstructed natural-scale marginals separately:
the upper blue half is *in vivo* and the lower pink half is *in vitro*.
Pair-specific DE initial populations have gray dashed outlines, whereas the
500 optimizer endpoints have C01/C02/C03 solid outlines. Density peaks are
normalized separately within each displayed distribution. The 5th--95th
percentile spans are descriptive numerical-search summaries—not posterior or
confidence intervals, biological replicates, or proof of structural
identifiability. Paired `log2(in vivo / in vitro)` products are retained only
as a secondary directional/search audit and are not the Figure 5D coordinate.
The full data and interpretation contract is in
`data/Figures/Figure5/figure5f_chart_contract.md`.

## Targeted Figure 6 supplementary rerun

Figure 6 and its parent-indexed Supplementary Figures 6-1 through 6-3 can be regenerated entirely from
the analytical inputs packaged under this iteration2 workspace:

```bash
cd /absolute/path/to/workspace
Rscript Code/Figures/data_Figure6.R --n-core=8 --rebuild=FALSE --n-resample=100
Rscript Code/Figures/draw_Figure6.R
Rscript Code/Figures/data_Supp_Figure6_1.R
Rscript Code/Figures/draw_Supp_Figure6_1.R
Rscript Code/Figures/data_Supp_Figure6_2.R
Rscript Code/Figures/draw_Supp_Figure6_2.R
Rscript Code/Figures/data_Supp_Figure6_3.R --n-core=8 --rebuild=FALSE
Rscript Code/Figures/draw_Supp_Figure6_3.R
```

The analysis ranks the 500 endpoints separately within each of six joint
warm-start pairs, uses the lowest 10% as the primary equal-weight ensemble,
and evaluates nested 5% and 20% sensitivity sets. It also reports exact
14-parameter endpoint multiplicity and repeats the qualitative-claim audit
after weighting each unique parameter endpoint once. Per-seed caches make the
calculation restartable. These entrypoints do not source scripts from another
checkout.

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
`audit/manifests/allowed_scientific_inputs.txt`. Model calculations may call
the repository's `oxygen/code/O2_supply_demand_MAP`; no other external code is
sourced at runtime. The resolved paths are recorded in
`audit/parameters/pipeline_parameters.tsv`.

The immutable validation baseline is
`Code/config/manifests/expected_scientific_input_md5.tsv`. It identifies each
input by a functional root ID and a relative path, so it does not depend on the
workspace name or location. `audit/reports/input_md5_validation.tsv` records
the latest validation summary. No MD5 validation is performed for copied
inputs, generated intermediates, code, figures, manuscript bootstrap files, or
published figure copies.
