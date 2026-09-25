# Joint fixed-O2 steady-state ploidy figures

This consume-only visualization module renders the materialized tables from
`analysis/joint_fixed_o2_ploidy_classification/`.

| File | Figures |
|---|---|
| `plot_joint_fixed_o2_ploidy_curves.R` | All regression-smoothed seed curves faceted by pair and fixed-O2 curve class, plus all pairs pooled by curve class. Curve class is encoded by facets; CatA/B/C is the only color encoding. Pair/class and Cat/class median curves are overlaid from analysis-owned tables. |

Both figures are written as 300-dpi PNG and vector PDF and registered in
`visualization_manifest.tsv`. This layer does not classify curves or calculate
summary statistics.
