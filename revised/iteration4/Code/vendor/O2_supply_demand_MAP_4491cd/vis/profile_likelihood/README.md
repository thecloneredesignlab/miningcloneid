# Profile-likelihood visualization

| File | Responsibility | Inputs | Outputs |
|---|---|---|---|
| `plot_sigma_burden_live_effective_pms.R` | Draw the `p_misseg` versus live-cell effective `p_ms` violin/boxplot. | `sigma_burden_p_misseg_vs_live_cell_plot.tsv` and the completed analysis manifest. | One PDF and a visualization manifest. |

This visualization layer does not read fitted parameters, RDS fit objects, or
model code.
