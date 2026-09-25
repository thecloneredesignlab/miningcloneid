# Joint fixed-O2 steady-state ploidy classification

This analysis-only module joins already-materialized joint-fit fixed-O2
dominant-attractor curves to the exact synthetic-seed manifest and the
pair-level CatA/B/C assignments. It does not refit the model, recompute
attractors, or draw figures.

| File | Exact responsibility | Main outputs |
|---|---|---|
| `analyze_joint_fixed_o2_ploidy_classes.R` | Validates the dense fixed-O2 grid and the one-to-one seed-to-pair mapping, attaches the existing `loess_persistent_v1` curve class and temporal ploidy Cat, and computes pair- and Cat-balanced summaries. | Annotated seed/curve tables, pair class composition and dominant class, Cat-balanced class composition, pair and Cat median curves, input-quality/provenance tables, and manifest. |

The two classifications remain distinct: `curve_class` describes the shape of
steady-state dominant mean ploidy across fixed O2, whereas CatA/B/C describes
the pair's separate 0-1000 day in-vivo ploidy trajectory. Cat is used as an
annotation and visualization color, not as the fixed-O2 classifier.
