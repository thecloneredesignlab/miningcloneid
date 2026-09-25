# Profile-likelihood runners

| File | Responsibility | Inputs | Outputs |
|---|---|---|---|
| `run_live_effective_pms_comparison.R` | Orchestrate the analysis-only cross-seed comparison followed by its visualization-only consumer. | Existing per-seed `simulation/invivo/cin/live_effective_pms` artifacts and comparison CLI settings. | Analysis tables, optional figure/visualization manifest, and a pipeline manifest. |

The runner does not fit a model or materialize missing per-seed simulations.
