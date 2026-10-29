# Figure7: fitted native p_misseg condition

This adds `Figure7_fitted_p_misseg` and `supp_fig7_8_fitted_p_misseg` without replacing fixed-p products. It reuses the current locked model and original numerical functions. The original numerical entry and its authenticated SHA256 remain unchanged.

## Unified dispatch

Use `Code/hpc/run_Figure78.sh --figure=7 --condition-mode=fitted_p_misseg --stage=STAGE --campaign-root=ABSOLUTE_ROOT`. On a native local R installation set `FIGURE78_NATIVE_RSCRIPT` to the Rscript executable. The dispatch selects `run_Figure7_fitted_p_misseg.R`; supported explicit stages are `plan`, `prepare`, `smoke`, `simulate`, `acceptance`, `aggregate`, `validate`, and `plot`.

`smoke` is a separate 30-day audit-only run. It is never publication eligible. Full simulation cannot start until the same campaign's small acceptance and grid gate passes. Representative workers cover C01/2N/0.5% and C02/4N/5% in both environments, at all three candidate grid bounds and all three beta values. Full workers reuse completed representative results.

HPC submission uses `Code/hpc/submit_Figure7_fitted_p_misseg.sh --campaign-root=ABSOLUTE_ROOT --phase=small|full --submit=yes`. It uses partition red, qos xxlarge, and no array concurrency cap. Every task has its own staged locked-model source/cache. Collection is explicit and never draws or promotes results.

## Inputs and calculation

Full winning joint fits: C01 seed94 and C02 seed80. Each context retains its own complete vector (26 in-vivo or 33 in-vitro parameters). Native p_misseg is read from the vector, never overridden or replaced by a constant effective probability. Fixed O2 and pure 2N/4N starts remain the approved diagnostic conditions. In-vitro configuration comes from `fit_result.rds$ctx$invitro_cfg`.

Daily A/B: viable low daughters from current high/low mothers, divided by their shared total. First-generation and permanent-lineage diagnostics retain their current definitions. Original native culture passage and last-five segment cycling continue to cumulative day1000 without an inoculum-feasibility stop. Culture summaries use three equal protocol routes, not optimizer seeds or stochastic CI.

The same relative-only mean mask is applied per route before aggregation: mass >=1e-8 of all live cells, identifiable origin, finite means; theoretical empty groups are NA. Expected count <1 is recorded only. All three culture routes are required. Counts, output shares and CDFs retain the existing rules.

## Data and rendering

The retained condition dataset is `data/Figures/Figure7/fitted_p_misseg/`, containing its own prepared bundle, tables, postprocess layer and provenance. Raw task/routes/runtime caches remain in the separately named campaign, never a required input for drawing. Render with `draw_Figure7_fitted_p_misseg.R --campaign-root=ABSOLUTE_RETAINED_ROOT --output-root=ISOLATED_OUTPUT`.

Main PDF: one page, 14x11 inches. Supplement PDF: 24 pages, 24x13 inches (12 time profiles and 12 unblocked CDF pages). Supplement PNG previews page1. Numerical validation, independent Python verification, complete PDF inspection, and original-file checksum checks precede publication. Original fixed-p figures, supplements, Figure8, and historical ver2 are not publication targets.
