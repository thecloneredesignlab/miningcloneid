# Fit-results simulation layer

| File | Responsibility | Inputs | Outputs |
|---|---|---|---|
| `materialize_extra_results_predictions.R` | Read completed per-seed prediction trajectories and materialize concrete ploidy and burden values without statistical comparison or plotting. It prefers `seed*/simulation/invivo/` and supports legacy `seed*/viz/` inputs. | Fit run directory with completed per-seed prediction TSVs. | Per-seed and cross-seed 1000-day ploidy/burden tables, prediction-gate values, input status, `simulation_manifest.tsv`. |
| `materialize_invivo_long_ploidy_metrics.R` | Extract the last available cohort/dose chromosome-number value at a requested horizon for every fitted seed. | Existing fit summaries and ploidy prediction TSVs. | `invivo_long_ploidy_metrics.tsv` and `simulation_manifest.tsv`. |

These scripts never fit parameters, source the protected model/optimizer, draw
figures, classify seeds statistically, or render reports.
