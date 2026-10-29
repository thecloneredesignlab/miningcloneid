# Figure7 conditional-mean postprocessing

This reporting-only layer affects the three conditional-mean time courses in
`supp_fig7_8`: per-cell viable low-daughter output, division rate and ploidy.
The policy is `Code/config/figure7_population_validity.json`: source-group mass
fraction of all live cells ≥1e-8, finite means and identifiable ancestry. There
is no absolute-count filter. Complete blocking with initially
empty H-to-L labels makes direct/historical groups structurally empty. Valid
zero rates are retained. All three culture routes must qualify before a median
and type-7 IQR are reported. The continuous in-vivo route is evaluated alone.
Initial empty labels and low groups that are impossible from pure4N roots are
also structurally empty. Saved log expected counts, below-one flags and count
availability are reporting fields only; they never mask otherwise valid means.

Raw model trajectories, simulation tables, scientific fingerprints and numerical
receipts remain immutable. The derived layer is retained in the single current
directory `data/Figures/Figure7/postprocess/population_validity/`. Its compact
`route_metrics.rds` contains the saved masses, log counts and unfiltered means
needed for subsequent postprocessing; neither redraw nor repeat postprocessing
requires the original HPC cache paths or an archive. `cache_sources.tsv` records
those original read-only sources as provenance, not runtime dependencies.

To repeat postprocessing from retained telemetry, choose a **new** staging
directory, then inspect `checks.tsv`, `primary_coverage.tsv`, `reason_summary.tsv`
and the three relative thresholds in `sensitivity_summary.tsv`:

```sh
Rscript Code/Figures/postprocess_Figure7_population_validity.R \
  --telemetry-root=data/Figures/Figure7/postprocess/population_validity \
  --output-root=audit/new_population_validity_candidate
```

The script refuses to overwrite an existing output directory. Only after checks
and sensitivity review should the derived directory be promoted with a fresh
backup and checksum verification. Then use the existing drawing entry:

```sh
Rscript Code/Figures/draw_Supp_Figure7_8.R \
  --products=supp_fig7_8 --output-root=audit/new_population_validity_figures
```

The renderer authenticates the original dataset and the separate postprocess
receipt, including its policy, code, telemetry and derived-table hashes. Stale or
missing validity data cause an error, never fallback to unfiltered means or a
simulation. Lines and ribbons have explicit segment groups to preserve NA gaps.

For a newly saved cache, `export_Figure7_population_cache.R` reads current native
culture route files and saved continuous telemetry. It never loads a model or
propagator. If continuous telemetry comes from an earlier retained cache, every
current mean/IQR and daily counter must match exactly before import; cache age
or matching fit IDs alone are insufficient. That one-time import is consolidated
into the compact current telemetry file above.

`test_Figure7_population_validity.R` tests threshold boundaries, structural
emptiness and valid zero rates without simulation. The independent Python
validator checks every stored route/summary mask, the 48-page output, publication
mirrors, and pixel identity of all 24 unchanged CDF pages. It emits page rasters
for visual review. None of these commands submits an HPC job.

Cell-share and low-daughter-output-share values and the CDF data are unchanged.
Figure7 A/B, Figure8 and historical `_ver2` outputs are outside the write scope.
Below-one flags use the existing saved model scale and imply no physical
extinction event. Relative sensitivity changes temporal reporting support, not
the values at retained points. See the Chinese and English calculation methods
for the scientific definitions and interpretation limits.
