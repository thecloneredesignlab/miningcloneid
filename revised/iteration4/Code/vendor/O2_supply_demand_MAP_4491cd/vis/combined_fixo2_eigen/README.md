# Combined FixO2 eigen-attractor visualization

- `plot_fixo2_eigen_attractor_embedding_curve_class.R` is a consume-only visualizer. It reads `pooled_embedding_curve_class_analysis_manifest.tsv` plus the annotated coordinate tables produced by analysis, writes the same PNG/PDF figure set, and writes the legacy `pooled_embedding_curve_class_manifest.tsv` used by the report.

It does not join curve classes, calculate slopes, or write scientific analysis tables.
