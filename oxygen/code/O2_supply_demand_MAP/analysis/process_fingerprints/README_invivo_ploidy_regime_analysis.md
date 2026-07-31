# In vivo ploidy regime analysis

This pipeline labels long-horizon chromosome-number trajectories and compares those labels against model-derived process fingerprints. It uses the current repository model helpers and C++ backend; it does not hard-code HPC result paths and does not modify fit result directories.

## Main entry

```bash
Rscript oxygen/code/O2_supply_demand_MAP/analysis/process_fingerprints/run_invivo_ploidy_regime_analysis.R \
  --run_dir oxygen/results/ver1/fit_invivo_O2_buffering_500seed__20260612_143000 \
  --out_dir oxygen/results/analysis \
  --objective_source auto \
  --objective_deltas 2,5,10 \
  --main_objective_delta 10 \
  --trajectory_end_day 1000 \
  --mid_window 200,700 \
  --late_window 900,1000 \
  --o2_grid 0.05,0.1,0.25,0.5,1,2,5 \
  --n_bootstrap 100 \
  --n_permutations 1000 \
  --workers 1 \
  --analysis_level core \
  --generate_viz false \
  --random_seed 20260623
```

Outputs are written under:

- `tables/`
- `figures/`
- `cache/`
- `logs/`
- `report/analysis_summary.md`

## Important constraints

- Process clustering matrices are built only from process fingerprints, not from actual 1000-day trajectory features.
- Trajectory features are used for phenotype labels, validation, and concordance only.
- Optimizer seeds are not posterior samples; seed proportions must not be interpreted as biological probabilities.
- Objective-compatible counterfactual interpolation and alternative-boundary sensitivity simulations require safe objective/simulation APIs. The core pipeline records unavailable diagnostics explicitly rather than fabricating those outputs.

## HPC examples

Smoke:

```bash
Rscript oxygen/code/O2_supply_demand_MAP/analysis/process_fingerprints/run_invivo_ploidy_regime_analysis.R \
  --run_dir /share/lab_crd/lab_crd/taoli/Project/miningcloneid_soft_coupling/oxygen/results/fit_invivo_O2_buffering_500seed \
  --out_dir /share/lab_crd/lab_crd/taoli/Project/miningcloneid_soft_coupling/oxygen/results/analysis/ploidy_regime_smoke \
  --objective_source auto \
  --objective_deltas 2,5,10 \
  --main_objective_delta 10 \
  --trajectory_end_day 1000 \
  --n_bootstrap 20 \
  --n_permutations 100 \
  --workers 1 \
  --max_seeds 10 \
  --analysis_level core \
  --generate_viz false \
  --random_seed 20260623
```

Full:

```bash
Rscript oxygen/code/O2_supply_demand_MAP/analysis/process_fingerprints/run_invivo_ploidy_regime_analysis.R \
  --run_dir /share/lab_crd/lab_crd/taoli/Project/miningcloneid_soft_coupling/oxygen/results/fit_invivo_O2_buffering_500seed \
  --out_dir /share/lab_crd/lab_crd/taoli/Project/miningcloneid_soft_coupling/oxygen/results/analysis/ploidy_regime_full \
  --objective_source auto \
  --objective_deltas 2,5,10 \
  --main_objective_delta 10 \
  --trajectory_end_day 1000 \
  --mid_window 200,700 \
  --late_window 900,1000 \
  --o2_grid 0.05,0.1,0.25,0.5,1,2,5 \
  --n_bootstrap 500 \
  --n_permutations 2000 \
  --workers 1 \
  --analysis_level full \
  --module_factorial false \
  --n_local_perturbations 0 \
  --generate_viz false \
  --random_seed 20260623
```

`--workers` is accepted for command compatibility; the current implementation runs the model probes in one R process. Full-mode counterfactual swaps, factorial swaps, and local perturbations are recorded as unavailable unless a safe objective/trajectory simulation API is wired.
