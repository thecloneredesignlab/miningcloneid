# Reanalysis tables and data dictionary

These files were generated from the supplied `top10.zip`, `joint_pre.zip`, and the output semantics defined by the Repo branch `HypoxiaLTEEFigures`.

## Main summary

- `analysis_summary.json`: machine-readable headline ranges and counts.
- `reanalysis_stdout.json`: script stdout capture.

## Separate in vitro

- `invitro_top10_robustness.csv`: objective, 4N decline, dead-compartment fractions, and late 2N prediction per selected seed.

## Separate in vivo

- `invivo_top10_fit_quality.csv`: burden RMSE, terminal mean/distribution distances, necrosis MAE and reconstructed necrosis objectives.
- `independent_parameter_values_long.csv`: natural-scale parameter values across selected separate fits.
- `independent_parameter_boundary_hits.csv`: bound-hit audit.
- `provenance_audit.csv`: available commit/worktree provenance.

## Joint selected solutions

- `joint_selected60_fit_summary.csv`: objective, optimizer and total soft-penalty diagnostics.
- `joint_pair_summary.csv`: pair-level objective ranges.
- `joint_soft_coupling_all60.csv`: all soft-coupled parameter rows.
- `joint_soft_coupling_parameter_summary.csv`: parameter direction, ratio and saturation summary.
- `joint_survival_function_all60.csv`: derived s(N) values.
- `joint_survival_function_summary_all60.csv`: derived survival-function summary.
- `joint_missegregation_function_all60.csv`: derived p_mis(N,O2) values.
- `joint_missegregation_function_comparison_all60.csv`: context ratios at reference O2/N.
- `joint_warm_start_movement_all60.csv`, `joint_warm_start_movement_by_seed.csv`: optimizer movement from warm starts.

## joint_pre landscape and fixed-O2

- `landscape_primary_cluster_summary.csv`
- `landscape_primary_cluster_silhouette.csv`
- `landscape_subcluster_summary.csv`
- `landscape_subcluster_silhouette.csv`
- `joint_warm_start_fixed_o2_crosswalk.csv`
- `fixed_o2_by_seed_500.csv`
- `fixed_o2_class_counts.csv`
- `fixed_o2_reliability_counts.csv`
- `fixed_o2_final_interpretation_counts.csv`
- `fixed_o2_class_by_reliability.csv`
- `fixed_o2_cluster_class_association.csv`
- `fixed_o2_landscape_cluster_crosswalk.csv`
- `fixed_o2_objective_definition_audit.csv`
- `fixed_o2_objective_spearman_correlations.csv`
- `fixed_o2_run_arguments.csv`

## File inventories

- `joint_pre_file_inventory.txt`: all files in the preprocessing archive.
- `joint_pre_AUC_file_inventory.txt`: search result for files whose names indicate AUC inputs/results; no definitive Figure 4 seed-level AUC table was found.

## Reproduction

Extract `top10.zip` and `joint_pre.zip`, then run from the delivery-package root:

```bash
python scripts/reanalyse_results.py \
  --top10-root /path/to/extracted/top10 \
  --joint-pre-root /path/to/extracted/joint_pre \
  --out-dir analysis
```

`--top10-root` must contain `top10_index.tsv`. `--joint-pre-root` must contain `cross_validation/` and `landscape_subcluster/`. The two input paths can alternatively be supplied through the `TOP10_ROOT` and `JOINT_PRE_ROOT` environment variables.
