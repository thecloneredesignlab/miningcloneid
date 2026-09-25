# Joint ploidy-coupling figures

This consume-only folder renders the plot-ready tables produced by
`analysis/joint_ploidy_coupling_association/`.

| File | Figures |
|---|---|
| `plot_joint_ploidy_categories.R` | Pair assignments with explicit estimability warning, CatA/B/C archetypes, pair medians/IQRs, deterministic representative trajectories with detected drops, CatB/C diagnostics, cutoff margins, and 81-scenario sensitivity composition. |
| `plot_ploidy_coupling_associations.R` | Pair-level Cat×parameter ratio intervals, pair-balanced Cat×Class profile, descriptive six-pair Cramér's V, in-vivo/in-vitro natural-value dumbbells, and the two-pair detail underlying each Cat. |

All figures are emitted as 300-dpi PNG and PDF and registered in
`visualization_manifest.tsv`. Numerical summaries remain owned by analysis.
