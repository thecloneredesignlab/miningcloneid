# Joint soft-coupling stability figures

This consume-only folder reads analysis TSVs and writes PNG/PDF pairs plus
`visualization_manifest.tsv`.

| File | Figures |
|---|---|
| `joint_coupling_visualization_utils.R` | Shared light theme, non-red/green Class/Cat palettes, biological parameter order, short pair labels, ClassB reference band, 300-dpi PNG/PDF saver, and chart-map manifest writer; no executable entrypoint. |
| `plot_joint_coupling_overview.R` | Two overview matrices: median signed ratio for every pair×parameter cell, and dominant ClassA/B/C with seed agreement and pair-level ploidy category. |
| `plot_within_pair_soft_coupling.R` | Within-pair A/B/C composition, entropy, ClassB-only `<1/=1/>1` direction, violin+box ratio distributions, and union/80/90/95/strict stability matrix. |
| `plot_between_pair_soft_coupling.R` | Pair-specific 5–95% ratio forests, cross-pair consensus, all-class rankings, stable-90 coverage, threshold/objective sensitivity, feasibility/boundary diagnostics, and UpSet-style pair-set membership. |
| `plot_soft_coupling_processes.R` | Parameter×Class heatmap faceted by audited biological process plus process-level ClassA/B/C profiles. |

These scripts do not read fit objects, recompute ratios, classify seeds, or
perform statistical aggregation. Every figure has an adjacent PNG and PDF and
an entry in `visualization_manifest.tsv` describing its chart family,
analytical question, dimensions, and source tables.
