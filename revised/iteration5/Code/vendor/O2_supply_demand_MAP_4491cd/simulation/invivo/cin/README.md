# In-vivo CIN simulation

| File | Responsibility | Inputs | Outputs |
|---|---|---|---|
| `invivo_cin_simulation.R` | Materialize CIN-related tables from an existing in-vivo fitted-parameter trajectory. | In-vivo simulation state supplied by the parent producer. | CIN tables used by downstream analysis and visualization. |
| `generate_live_effective_pms_outputs.R` | Compute live-cell weighted effective `p_ms` from fitted parameters and already materialized in-vivo O2, ploidy, and viability tables. | Seed `best_params.tsv`, `fit_config.rds`, and `simulation/invivo` manifest/tables. | Six scientific TSV tables plus schema and simulation manifest under `cin/live_effective_pms`. |

This folder does not create figures or render reports.
