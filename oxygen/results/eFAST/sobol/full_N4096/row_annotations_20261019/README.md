# Annotated full_N4096 summary figures

This display overlay updates 15 scope-summary/comparison figures (30 PDF/PNG files).
Per-seed/case plots and response surfaces retain their prior display versions.
No model evaluation or Sobol-index estimation is performed.

## Row annotations

1. Process: exact Figure 4 parameter-group mapping and palette copied here.
2. O2 correlation: exact Figure 4B categories, including its original magnitude gate; no reclassification.
3. Sobol O2 sensitivity: shown only in S1 and ST, using the existing ploidy-based categories without a 0.3 gate. Growth follows the same rows and annotations. This column never determines row order.

S2 has exactly two annotation columns. Their left/right halves are the first/second parameter named in the row. Earlier S2 classification statistics remain historical/source diagnostics only and are not displayed or used for ordering.

## O2 preference and row order

Low is [0,1]%, High is [3,5]%. Scores are normalized trapezoidal integrals in physical oxygen, independent of the symlog display. The entire curve is resampled as one unit: five scramble-repeat curves globally; 500 fitted-seed curves in neighborhoods, each already a five-repeat mean. At every bootstrap draw the neighborhood statistic is recomputed as the pointwise median, then integrated. This is not the median of seed-specific window averages.

There are 5,000 paired bootstrap draws, RNG seed 5826. The two-sided bootstrap sign-tail score is computed with a +1 correction and adjusted by Benjamini-Hochberg within scope/index (14 S1, 14 ST or 91 S2 rows). q<0.05 and the observed contrast determine High/Low. Both window scores nonpositive implies no resolved preference, since negative sensitivity estimates are sampling artifacts. No positive magnitude threshold is introduced. Raw estimates remain unchanged.

These are exploratory grouping diagnostics, following the reference plot bootstrap procedure, not a calibrated significance test or proof of a nonzero interaction. With five global repeats, grouping can be fragile. Neighborhood resampling describes robustness over this optimizer-solution collection, not biological confidence or numerical convergence. Original resolution/repeat diagnostics remain separately available, and no seed is excluded. Gray means no resolved O2 preference, not no effect on the output. Tiny estimated contributions can acquire a directional category because the 0.3 gate is removed; always read the magnitude and diagnostics too.

Current ordering uses only the raw physical-window scores: first all rows with High >= Low, descending by High; then all rows with Low > High, ascending by Low. There is no middle/unresolved block and no classification, q-value, convergence or peak gate in the ordering. Exact equality joins the first block; exact score ties use fixed source order. Negative estimates remain signed. Standalone scope figures use their own ploidy scores. All cross-scope comparisons use neighborhood ploidy scores. S1/ST annotation colors may interleave because the categories do not control sorting.

## Files

- plan.json: original images, numerical/source identities and protected output identities.
- original_figures/: complete pre-annotation image backups.
- classification.tsv / classification.json: group scores, bootstrap diagnostics, resampling units and original diagnostic coverage.
- row_order.tsv: all scope/index orders, including 36/91-pair orders.
- bootstrap_rows/: reusable classification checkpoints.
- figure_receipts/: displayed row identities, member order, annotation colors, data hashes, color limits and image hashes.
- validation.json / completion.json: numerical, image, annotation, sort, PDF and completeness checks.
- visual_review.json: final PDF raster pages inspected visually.

The main PDFs remain at their existing global/figures, neighborhood10pct/figures and figure_comparison_20261019 paths. Older axis/redraw receipts are historical layers; this directory records the latest summary display. The original symlog axis, oxygen boundaries, raw values and color limits are retained.

## Reproduce

```sh
python3 oxygen/code/O2_supply_demand_MAP/analysis/eFAST/sobol/redraw_full_n4096_annotated.py \
  --root oxygen/results/eFAST/sobol/full_N4096 --workers 6
```

Python, NumPy, Matplotlib and pdfplumber are required. Phases are classify, render, validate or all. Classification resumes from source-identified checkpoints. Temporary curve caching is under /private/tmp. Original scientific numerical tables are never rewritten.

## Verified outcome

All 15 summary/comparison figures (30 files) completed and passed validation.
All 2,107 pre-existing numerical-file hashes, 3,012 other oxygen-axis images and
108 response-surface images are unchanged. Annotation alignment error is zero. S1/ST have three columns and S2 has two.
Independent full-curve bootstrap integration and score-only sort checks passed, including deliberately conflicting third-column labels and signed/equal-score cases.
The shared and standalone color limits match the earlier scale source table;
see `color_limit_validation.json`. The final 12-panel group, global focal S2, neighborhood S1 and neighborhood
91-pair S2 PDF rasterizations passed visual review. Header/title spacing is
checked quantitatively for every panel.

The previous category-driven order and all 30 previous plot files are preserved under `iterations/before_score_only_order/`. Current `row_order.tsv` records the score comparison side rather than a sensitivity category.
