# Profile-likelihood analysis

This folder contains analysis-only consumers. Numerical state derived from
fitted parameters is created under `simulation/`; figures are created under
`vis/`.

| File | Responsibility | Inputs | Outputs |
|---|---|---|---|
| `collect_profile_likelihood_results.R` | Collect and summarize completed profile-likelihood fit results. | Existing profile-likelihood result files. | Profile-likelihood summary tables. |
| `compare_sigma_burden_live_effective_pms.R` | Compare fitted `p_misseg` with live-cell effective `p_ms` across materialized seeds and `sigma_burden` caps. | Each seed's `simulation/invivo/cin/live_effective_pms` manifest and summary tables. | Task, by-seed, by-cap, plot-ready, schema, and analysis-manifest TSV files. |
| `estimate_live_effective_pms.R` | Deprecated compatibility wrapper for the migrated simulation producer. | Same CLI arguments as the producer. | Outputs from `simulation/invivo/cin/generate_live_effective_pms_outputs.R`. |

Plot the comparison output with
`vis/profile_likelihood/plot_sigma_burden_live_effective_pms.R`.
For a single command that runs the separated analysis and visualization stages,
use `runner/profile_likelihood/run_live_effective_pms_comparison.R`.
