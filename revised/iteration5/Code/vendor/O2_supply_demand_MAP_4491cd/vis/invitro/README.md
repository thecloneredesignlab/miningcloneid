# In-vitro visualization library

## `o2_supply_demand_map_invitro_plot_utils.R`

- Entrypoint: no; source as a library.
- Single responsibility: return in-vitro ggplot/patchwork objects from
  precomputed lineage, daily-count, ploidy, distribution, and flow tables.
- Required packages at call time: `ggplot2`, `scales`; selected composites may
  use `patchwork` through the caller.
- Reads files: none.
- Writes files: none.
- Runs model/simulation: no.
- Public API:
  - `ivt_ploidy_fraction_fill_scale()`
  - `ivt_plot_daily_counts()`
  - `ivt_plot_distribution_heatmap()`
  - `ivt_plot_lineage_counts()`
  - `ivt_plot_lineage_flow_density()`
  - `ivt_plot_lineage_growth()`
  - `ivt_plot_lineage_ploidy()`
- Caller: `../viz_invitro_model_O2_supply_demand_MAP_results.R`.
- Migration source: former `oxygen/code/in-vitro-utils/plotting.R`.
- Regression: API/formals contract plus seed10 PDF/PNG/manifest comparison.
