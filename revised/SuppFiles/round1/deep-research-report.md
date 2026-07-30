# Reprocessing Audit of the LTEE Joint-Fitting Top10 Archive

## Executive summary

I inspected the uploaded `top10.zip`, the public `HypoxiaLTEEFigures` branch of your repository, and the O2 workflow documentation/code paths referenced there. The repository clearly defines an R-based oxygen supply-demand fitting workflow with a unified optimizer, a joint-fit wrapper, a YAML configuration file, and documented per-seed outputs for in vivo, in vitro, and joint fits. The joint workflow specifically expects warm-start soft-coupling tables and produces per-seed diagnostic tables such as `joint_soft_coupling.tsv`, `joint_components.tsv`, `joint_soft_coupling_initial_values.tsv`, `joint_warmup_initial_values.tsv`, plus configuration/result RDS files. citeturn2view0turn3view0turn7view0turn8view0turn9view0turn12view0turn13view0

From direct inspection of the archive, `top10.zip` is largely internally consistent as a curated top-results bundle: it contains 80 ranked entries, and every objective recorded in `top10_index.tsv` exactly matches the corresponding seed-level `fit_summary.tsv`. The most important archive problems are packaging and portability problems rather than ranking errors: the bundle contains a large macOS metadata payload, all six joint pairings are missing the run-level `joint_soft_coupling_tables/joint_soft_coupling_parameters_table__<label>.csv` files that the resolved configs point to, the top-level index uses machine-specific absolute paths, and all 60 joint seed directories are missing `fit_result.rds`. The first three issues were repairable; the last one was not repairable without rerunning the joint fitter. 

Because `joint_pre.zip` was not present in the uploaded materials and the current runtime does not have R installed, I could not perform a faithful numerical rerun of the repository’s R workflow here. The repo documentation and optimizer code make clear that the workflow is R-based and depends on specific in vivo and in vitro inputs and R packages, so a true refit requires that environment and the original preprocessed inputs. citeturn4view1turn7view0turn12view0

I produced a corrected archival deliverable that fixes the repairable issues and documents the unresolved ones:

- [Corrected archive](sandbox:/mnt/data/top10_corrected.zip)
- [Corrections manifest](sandbox:/mnt/data/top10_corrected/top10/corrections_manifest.tsv)
- [Reconstructed joint-table manifest](sandbox:/mnt/data/top10_corrected/top10/reconstructed_joint_soft_coupling_tables.tsv)
- [Execution log](sandbox:/mnt/data/top10_corrected/top10/execution_log.tsv)
- [Command log](sandbox:/mnt/data/top10_corrected/top10/commands_run.txt)

## Repository-grounded workflow and expected outputs

On the `HypoxiaLTEEFigures` branch, the repository root contains both `manuscript/` and `oxygen/`; the `oxygen/` subtree is the relevant workflow for the uploaded fitting results. The repository’s own README identifies the main paths for this workflow: the unified optimizer is `code/O2_supply_demand_MAP/optimizer/fit_model_O2_supply_demand_MAP.R`, the local runner is `code/O2_supply_demand_MAP/runner/run_fit_model_O2_supply_demand_MAP.sh`, the joint-only wrapper is `code/O2_supply_demand_MAP/runner/run_fit_joint_model_O2_supply_demand_MAP.sh`, and the default config is `config/O2_supply_demand.yaml`. citeturn2view0turn7view0

The joint wrapper script is intentionally thin: it exports provenance variables, checks whether `--mode` was supplied, and then dispatches to the unified runner with `--fit_joint`, defaulting to `--mode=run` when no mode was provided. That means the expected behavior of your archived joint runs is governed by the same unified runner/optimizer/config stack as the separate in vivo and in vitro runs. citeturn3view0

The repo README describes the expected run-directory layouts. For in vivo multi-seed runs, the runner creates `<out_root>/<run_prefix>/seed<seed>/`. For in vitro runs, the main outputs include `fit_summary.tsv`, `best_params.tsv`, `best_params_transformed.tsv`, and `fit_result.rds`. For joint runs, the runner resolves a seed plan and, for each seed, performs joint fitting, in vivo visualization, in vitro visualization, and joint report rendering. The joint objective is the weighted sum of in vivo and in vitro objectives plus soft-coupling and constraint penalties. citeturn8view0turn9view0

The config file further defines the expected computational behavior. In the repo default, joint soft coupling is enabled, the soft-coupled parameters are explicitly listed, `joint_soft_coupling_sigma_default` is `0.65`, `joint_soft_coupling_welsch_c` is `0.4`, warm-starting is enabled, and `joint_warmup_sigmaN` defaults to `0.1`. The same config also controls report density (`viz_report_dt`) and the number of top-ranked seeds highlighted in auto-generated reports (`viz_top_n`). citeturn12view0

The README also explains the warm-start/start-table mechanism that matters most for your joint archive. A labeled joint soft-coupling start table is generated by `make_joint_soft_coupling_parameters_table.R`; by default, it is written as `oxygen/data/O2_supply_demand/joint_soft_coupling_parameters_table__<seed_label>.csv`, and it contains `param_name`, `value`, `scale`, `seed_label`, `invivo_seed_label`, and `invitro_seed_label`. During joint fitting, per-seed outputs can include `joint_soft_coupling.tsv`, `joint_components.tsv`, `fit_summary.tsv`, `joint_soft_coupling_initial_values.tsv`, `joint_warmup_initial_values.tsv`, `fit_config.rds`, and `fit_result.rds`. The report renderer reads `fit_summary.tsv`, `best_params.tsv`, parameter-table snapshots, and optional `joint_soft_coupling.tsv` to annotate and render the report. citeturn9view0turn13view0

```mermaid
flowchart LR
    A[Repo oxygen README and config] --> B[Infer expected runners, inputs, and outputs]
    B --> C[Inspect top10.zip contents and seed directories]
    C --> D[Check index-objective consistency]
    D --> E[Compare expected files to actual files]
    E --> F[Repair portable/archive-level issues]
    F --> G[Emit corrected archive and manifests]
```

## Inspection method and archived run parameters

I unpacked `top10.zip`, enumerated all result groups, parsed `top10_index.tsv`, inspected representative `config.input.yaml`, `config.resolved.yaml`, `run_provenance.tsv`, `run_command.txt`, `fit_summary.tsv`, and status logs, and compared observed files against the output schema described in the repo documentation.

The archive contains three result groups:

| Result group | Contents | Count |
|---|---:|---:|
| `fit_invitro_O2_buffering_500seed` | top 10 in vitro seeds | 10 |
| `fit_invivo_O2_buffering_500seed` | top 10 in vivo seeds | 10 |
| `fit_joint_multi_warmup_tsne_sigmaN0p1216_objmin_500seed_20260714_021540` | six warm-start pairings × top 10 joint seeds each | 60 |

Within the joint group, the six pairings are:

| Pair label | Best objective | Worst objective | Spread |
|---|---:|---:|---:|
| `fit_joint_tsne_vi_seed138_C03Sc01_vt_seed10` | 19.978236 | 20.019243 | 0.041007 |
| `fit_joint_tsne_vi_seed25_C02Sc01_vt_seed10` | 18.970464 | 18.972960 | 0.002496 |
| `fit_joint_tsne_vi_seed290_C01Sc02_vt_seed10` | 19.791314 | 19.805523 | 0.014209 |
| `fit_joint_tsne_vi_seed311_C03Sc02_vt_seed10` | 19.414487 | 19.434210 | 0.019723 |
| `fit_joint_tsne_vi_seed322_C02Sc02_vt_seed10` | 18.890060 | 18.925572 | 0.035513 |
| `fit_joint_tsne_vi_seed366_C01Sc01_vt_seed10` | 18.852286 | 18.866626 | 0.014341 |

The archived provenance shows these joint runs were produced from a Git working tree named `miningcloneid_soft_coupling` on branch `soft_coupling` at commit `3d4985997c54f34c3b1a30cc72462a282cff4817`. The resolved joint configs show that the archived runs overrode several defaults from the repo config: they used `seeds_csv=500`, `n_cores=22`, `NP=80`, `itermax=500`, `append_run_prefix_timestamp=no`, and `joint_warmup_sigmaN=0.1216`, while keeping soft coupling enabled and retaining the same global Welsch penalty settings documented in the repo config. The repo default config confirms the baseline settings that were overridden in the archived run objects. citeturn12view0

Just as importantly, every row in `top10_index.tsv` matched the corresponding objective recorded in the seed’s `fit_summary.tsv`. That means the ranking logic in the exported archive is correct; the archive’s problems are about portability, completeness of referenced support files, and reproducibility, not about the top-10 selection itself.

## Discrepancy analysis

The expected-vs-actual comparison is clearest when separated into structural integrity, metadata portability, and missing computational state. The repo documentation establishes what a joint run should point to and which warm-start/joint diagnostics should exist; the archive inspection shows which of those expectations were satisfied and which were not. citeturn8view0turn9view0turn13view0

| Path or class | Expected from repo/code | Actual in uploaded archive | Corrected status |
|---|---|---|---|
| `__MACOSX/` subtree | Not part of repo outputs | Present in large quantity inside zip | Removed |
| `.DS_Store` files | Not part of repo outputs | Present throughout archive | Removed |
| `top10/top10_index.tsv` | Portable manifest if reused externally | `source_path` points to host-specific absolute paths under `/Users/.../GitHub/soft_coupling/...` | Rewritten to archive-relative paths; original preserved as `top10_index.original.tsv` |
| `fit_joint.../joint_soft_coupling_tables/joint_soft_coupling_parameters_table__<label>.csv` | Warm-start labeled CSV expected by joint config/start-table workflow | Missing for all six joint pairings, even though each resolved config points to one | Reconstructed for all six pairings |
| Joint seed `fit_summary.tsv`, `joint_soft_coupling.tsv`, `joint_components.tsv`, `joint_soft_coupling_initial_values.tsv`, `joint_warmup_initial_values.tsv`, `fit_config.rds` | Expected joint diagnostics and summaries | Present | No change needed |
| Joint seed `fit_result.rds` | Repo documents it as a possible joint per-seed output | Missing for all 60 joint seeds | Unresolved without rerun |

The missing run-level joint soft-coupling tables were the most important repairable defect. The repo documents that the start-table CSV is generated once per warm-start label and then fed into the joint run. In your archive, every seed directory within a given joint pair contains a `joint_soft_coupling_parameters_table_input.csv`, and all 10 copies within each pair are byte-identical. That made it safe to reconstruct the six missing run-level CSVs under `joint_soft_coupling_tables/` using the already-exported per-seed copies. This aligns with the repo’s documented column contract for those CSVs and with the resolved-config paths embedded in the joint run directories. citeturn9view0

The unresolved defect is the absence of `fit_result.rds` from all 60 joint seed directories. The repo documentation says a joint seed directory can contain both `fit_config.rds` and `fit_result.rds`; your archive contains the former but not the latter. I did not fabricate replacements, because any synthetic placeholder would be misleading and non-reproducible. A faithful repair would require rerunning the joint fitter in the repo’s R environment with access to the original inputs. citeturn9view0

A smaller but still real issue is portability of the top-level index. The repo README uses explicit repo-root-relative command patterns and output-root conventions; absolute local filesystem paths in an exported archive are therefore less useful than archive-relative paths for downstream analysis and sharing. I rewrote `source_path` accordingly while preserving the original manifest. citeturn7view0turn8view0

The status logs did not indicate numerical failure in the exported top-10 fits. The recurring warning was report-generation related: `pandoc unavailable`, which led to HTML-only report fallbacks. That is a presentation limitation, not a fitting failure.

## Corrected outputs and file mapping

The corrected deliverable is a repaired top-results archive, not a newly refit numerical result set. Its purpose is to make the exported subset portable, internally self-consistent, and more aligned with the repo’s documented warm-start file conventions without inventing missing optimizer state.

The file mapping is:

| Original state | Corrected state |
|---|---|
| `top10/top10_index.tsv` with absolute `source_path` values | `top10/top10_index.tsv` rewritten to archive-relative paths |
| no preserved original manifest | `top10/top10_index.original.tsv` added |
| missing `top10/fit_joint.../joint_soft_coupling_tables/joint_soft_coupling_parameters_table__<label>.csv` | six reconstructed CSVs added under the expected run-level directory |
| scattered `.DS_Store` files | removed |
| entire `__MACOSX/` metadata subtree | removed |
| absent joint `fit_result.rds` | still absent; listed as unresolved in manifest |

I also added archive-level documentation files:

| Added file | Purpose |
|---|---|
| `top10/README_corrections.md` | Human-readable summary of repairs and unresolved issues |
| `top10/corrections_manifest.tsv` | Machine-readable list of fixed and unresolved discrepancies |
| `top10/reconstructed_joint_soft_coupling_tables.tsv` | Exact provenance of each reconstructed joint start table |
| `top10/execution_log.tsv` | Environment snapshot |
| `top10/commands_run.txt` | Command-level log of the repair steps |

The reconstructed start tables are especially well-grounded because the repo defines the expected columns for these CSVs, and the report/fit logic is designed around the same `param_name` / `value` / `scale` convention. citeturn9view0

Downloadable repaired artifacts:

- [Corrected archive](sandbox:/mnt/data/top10_corrected.zip)
- [Corrections manifest](sandbox:/mnt/data/top10_corrected/top10/corrections_manifest.tsv)
- [Reconstructed joint-table manifest](sandbox:/mnt/data/top10_corrected/top10/reconstructed_joint_soft_coupling_tables.tsv)

## Execution log and reproducibility fixes

This inspection-and-repair pass was performed in a Linux container with Python `3.13.5`. R was **not** installed in the runtime, which was the main reason a true repo-native rerun could not be executed here. That matters because the repo explicitly defines this workflow as R-based, names R packages that must be installed, and routes both fitting and reporting through R scripts. citeturn7view0turn12view0

The practical command flow I executed was:

```bash
unzip -q /mnt/data/top10.zip -d /mnt/data/top10_unz
# inspect file tree, indexes, configs, logs, and seed outputs
# verify objective consistency between top10_index.tsv and fit_summary.tsv
# reconstruct missing run-level joint_soft_coupling_tables/*.csv from identical per-seed inputs
# rewrite top10_index.tsv to archive-relative paths and preserve original
zip -qr /mnt/data/top10_corrected.zip /mnt/data/top10_corrected
```

The repo-native command structure you would need for a full rerun, once the original inputs and an R environment are available, follows the repository documentation closely:

```bash
# from repository root

# generate a labeled joint soft-coupling start table
Rscript oxygen/code/O2_supply_demand_MAP/analysis/warm_start/make_joint_soft_coupling_parameters_table.R \
  --invivo-seed-dir oxygen/results/fit_invivo_O2_buffering_500seed/seed<best_invivo_seed> \
  --invitro-seed-dir oxygen/results/fit_invitro_O2_buffering_500seed/seed<best_invitro_seed> \
  --seed-label <label>

# launch the joint run
bash oxygen/code/O2_supply_demand_MAP/runner/run_fit_joint_model_O2_supply_demand_MAP.sh \
  --mode=run \
  --config=oxygen/config/O2_supply_demand.yaml

# aggregate postprocessed results
Rscript oxygen/code/O2_supply_demand_MAP/analysis/fit_results/extra_results.R \
  --run_dir=oxygen/results/<joint_run_prefix>
```

Those command patterns, as well as the existence and meaning of the start-table CSV and extra-results workflow, are documented in the repo README. The optimizer code also validates that the in vivo dataset exists under `data_dir` and specifically checks for the burden workbook `dt_Gem_VT_20260209_v5.xlsx`, which is another reason the missing `joint_pre.zip` / preprocessed-input bundle is a real blocker for full rerun reproducibility here. citeturn4view1turn8view0turn9view0

The most important reproducibility fixes are straightforward:

1. Export the actual preprocessed input bundle used for the run alongside the top-results archive.
2. Include the run-level `joint_soft_coupling_tables/` directory whenever exporting joint warm-start runs.
3. Keep `top10_index.tsv` portable by storing archive-relative or repo-relative paths, not machine-specific absolute paths.
4. Decide whether `fit_result.rds` is required for downstream reproducibility; if yes, export it consistently for joint seeds.
5. Capture `sessionInfo()` and package versions in each run root, because the workflow is R-based and report generation already has environment-sensitive behavior such as pandoc-dependent output mode. citeturn7view0turn9view0

Overall, the archive is scientifically usable as a ranked-result subset, but before it can serve as a fully reproducible joint-fitting deliverable, it still needs one real computational repair: regeneration of the missing joint `fit_result.rds` files in the proper R environment with the original preprocessed inputs.