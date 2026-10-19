# Annotated full_N4096 summary figures

This display overlay updates 15 scope-summary/comparison figures (30 PDF/PNG files).
Per-seed/case plots and response surfaces retain their prior display versions.
No model evaluation or Sobol-index estimation is performed.

## Three row annotations

1. Process: exact Figure 4 parameter-group mapping and palette copied here.
2. O2 correlation: exact Figure 4B categories, including its original magnitude gate; no reclassification.
3. Sobol O2 sensitivity: separately recomputed for current S1, ST and S2, using ploidy curves. Growth follows the same rows and annotations. No 0.3 magnitude gate is used in column 3.

For S2 columns 1/2, the left half is the first parameter named in the row, and the right half is the second. Column 3 describes the pair itself. All 91 pairs are classified together; the focal 36-pair subset retains the same classifications and adjusted values.

## O2 preference and row order

Low is [0,1]%, High is [3,5]%. Scores are normalized trapezoidal integrals in physical oxygen, independent of the symlog display. The entire curve is resampled as one unit: five scramble-repeat curves globally; 500 fitted-seed curves in neighborhoods, each already a five-repeat mean. At every bootstrap draw the neighborhood statistic is recomputed as the pointwise median, then integrated. This is not the median of seed-specific window averages.

There are 5,000 paired bootstrap draws, RNG seed 5826. The two-sided bootstrap sign-tail score is computed with a +1 correction and adjusted by Benjamini-Hochberg within scope/index (14 S1, 14 ST or 91 S2 rows). q<0.05 and the observed contrast determine High/Low. Both window scores nonpositive implies no resolved preference, since negative sensitivity estimates are sampling artifacts. No positive magnitude threshold is introduced. Raw estimates remain unchanged.

These are exploratory grouping diagnostics, following the reference plot bootstrap procedure, not a calibrated significance test or proof of a nonzero interaction. With five global repeats, grouping can be fragile. Neighborhood resampling describes robustness over this optimizer-solution collection, not biological confidence or numerical convergence. Original resolution/repeat diagnostics remain separately available, and no seed is excluded. Gray means no resolved O2 preference, not no effect on the output. Tiny estimated contributions can acquire a directional category because the 0.3 gate is removed; always read the magnitude and diagnostics too.

Row blocks: High O2, Low O2, unresolved. High rows descend by High-window score. Low rows ascend by Low-window score. Unresolved rows descend by all-grid peak. Ties use the fixed source row order. Each standalone scope uses its own ploidy-based order. All cross-scope comparison figures use neighborhood ploidy order; the global annotation still shows its own classification, so its colors need not form contiguous blocks.

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
108 response-surface images are unchanged. Annotation alignment error is zero.
Independent full-curve bootstrap integration and sort-direction checks passed.
The shared and standalone color limits match the earlier scale source table;
see `color_limit_validation.json`. The final 12-panel group, neighborhood group
and global 91-pair S2 PDF rasterizations passed visual review.
