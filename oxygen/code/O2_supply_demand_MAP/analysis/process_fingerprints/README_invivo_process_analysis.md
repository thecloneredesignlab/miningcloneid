# In vivo O2 process fingerprint analysis

This directory contains a reproducible analysis pipeline for clustering in vivo
single-seed O2_supply_demand_MAP fits by biological process fingerprints rather
than raw optimizer parameters alone.

Main entry point:

```bash
Rscript oxygen/code/O2_supply_demand_MAP/analysis/process_fingerprints/run_invivo_process_analysis.R \
  --run_dir oxygen/results/fit_invivo_O2_buffering_500seed \
  --out_dir oxygen/results/analysis \
  --objective_source auto \
  --objective_deltas 2,5,10 \
  --main_objective_delta 10 \
  --n_bootstrap 200 \
  --random_seed 20260623 \
  --workers 8 \
  --generate_viz false
```

The code does not hard-code local or HPC absolute paths. The input fit directory
is always supplied through `--run_dir`; all generated tables, figures, caches,
logs, and reports are written below `--out_dir`.

## Output layout

```text
<out_dir>/
  tables/
  figures/
  cache/
  logs/
  report/
```

Important tables include:

- `seed_manifest.tsv`, `seed_qc_summary.tsv`, `seed_exclusion_log.tsv`
- `parameter_values_long.tsv`, `parameter_matrix_raw.tsv`,
  `parameter_matrix_transformed.tsv`, `parameter_transform_metadata.tsv`
- `process_fingerprint_static_full_long.tsv`
- `process_fingerprint_static_18only_long.tsv`
- `process_fingerprint_dynamic_long.tsv`
- `process_fingerprint_static_scaled.tsv`
- `process_fingerprint_static_18only_scaled.tsv`
- `process_fingerprint_dynamic_scaled.tsv`
- `distance_parameter_matrix.tsv`, `distance_static_full_matrix.tsv`,
  `distance_static_18only_matrix.tsv`, `distance_dynamic_matrix.tsv`
- `pairwise_regime_diagnostics.tsv`
- `cluster_k_diagnostics.tsv`, `cluster_membership_static_full.tsv`,
  `cluster_medoids.tsv`, `cluster_qc_summary.tsv`
- `unavailable_dynamic_features.tsv`
- `report/analysis_summary.md`

## Model coupling

Static fingerprints source the current repository model file:

```text
oxygen/code/O2_supply_demand_MAP/model/model_O2_supply_demand_MAP.R
```

The pipeline calls the model's existing helpers for hypoxia stress,
proliferation, death, missegregation, buffering, WGD, and O2 supply-demand
responses. It does not reimplement a separate biological model.

The actual O2 scale is percent O2. For example, `1` means 1% O2, not a 0-1
fraction.

## Parameter aliases and transforms

The target 18 parameters are read from `best_params.tsv` when seed directories
are present. If seed directories are unavailable, the pipeline falls back to
`extra_results/seed_summary.tsv` `value__*` columns or a sibling
`*_param_matrix_with_seed.tsv` file.

Aliases are explicit in `process_fingerprint_utils.R`. Missing target
parameters are recorded and cause that seed to be marked invalid; missing values
are not silently replaced.

Transforms follow the current optimizer scale when known. In this model most
positive rates/scales and fitted probabilities are optimized as `log10_*`,
while `buffer_smax`, `gamma_growth`, `gamma_mu`, and `n_O` are identity scale.

## Dynamic fingerprints

Dynamic fingerprints are extracted only from existing prediction tables, such as:

```text
extra_results/predict1000_burden_total_seed_day_mean.tsv
extra_results/predict1000_ploidy_seed_day_mean.tsv
```

If O2 trajectories, entropy, or high-ploidy lineage fractions are not present in
the output schema, they are listed in `unavailable_dynamic_features.tsv`. They
are not imputed as zero.

## HPC commands

Small validation from the HPC repository root:

```bash
Rscript oxygen/code/O2_supply_demand_MAP/analysis/process_fingerprints/run_invivo_process_analysis.R \
  --run_dir /share/lab_crd/lab_crd/taoli/Project/miningcloneid_soft_coupling/oxygen/results/fit_invivo_O2_buffering_500seed \
  --out_dir /share/lab_crd/lab_crd/taoli/Project/miningcloneid_soft_coupling/oxygen/results/analysis \
  --objective_source auto \
  --objective_deltas 2,5,10 \
  --main_objective_delta 10 \
  --n_bootstrap 20 \
  --random_seed 20260623 \
  --workers 4 \
  --generate_viz false \
  --max_seeds 10 \
  --force true
```

Full run from the HPC repository root:

```bash
Rscript oxygen/code/O2_supply_demand_MAP/analysis/process_fingerprints/run_invivo_process_analysis.R \
  --run_dir /share/lab_crd/lab_crd/taoli/Project/miningcloneid_soft_coupling/oxygen/results/fit_invivo_O2_buffering_500seed \
  --out_dir /share/lab_crd/lab_crd/taoli/Project/miningcloneid_soft_coupling/oxygen/results/analysis \
  --objective_source auto \
  --objective_deltas 2,5,10 \
  --main_objective_delta 10 \
  --n_bootstrap 500 \
  --random_seed 20260623 \
  --workers 16 \
  --generate_viz false \
  --force true
```

## Interpretation boundary

Process clustering represents data-compatible mechanistic regimes; it does not
automatically equal true biological heterogeneity. Only differences that exceed
fitting uncertainty and predict independent data should be interpreted as real
biological mechanisms.

