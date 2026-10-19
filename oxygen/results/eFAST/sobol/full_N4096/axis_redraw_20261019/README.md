# Local oxygen-axis redraw, 2026-10-19

This is the current display overlay for all 1,521 oxygen-axis figures under
`full_N4096` (PDF and PNG: 3,042 files). It includes global/neighborhood scope
summaries, all 500 seed-specific S1/ST/S2 sets, both global case copies, and the
six comparison figure groups. The 54 parameter response surfaces have parameter
axes rather than oxygen axes; their 108 image files are preserved unchanged.

## Oxygen axis and regions

- Matplotlib `symlog`, base 10, linear threshold 0.025%, linear scale factor 1.
- Real O2=0 is retained. Values from 0 to 0.025% are linear; larger values are
  logarithmic. Zero is not replaced by an epsilon.
- Major ticks: 0,0.025,0.1,0.5,1,3,5%.
- Black dashed boundaries at 1% and 3%.
- Labels above each panel: Low, Mid, High.
- Region labels use a fixed 3-point offset above the axes, including 91-row
  long figures, so they remain separated from panel titles.
- Low: [0,1]%; intermediate: (1,3)%; high: [3,5]%.
- Heatmaps use `pcolormesh` with actual physical oxygen midpoint edges, clipped
  to [0,5]%, rather than treating uniformly spaced columns as log-spaced data.
  Original values, ordering and all 201 oxygen samples are unchanged.
- Expanding low oxygen visually does not add samples between 0 and 0.025%.
  Existing oxygen-band averages/classifications retain their original physical
  oxygen-grid weighting; heatmap area on a nonlinear axis is not a band integral.
- S1/ST positive estimates use white-to-purple. Existing numerical color limits
  are retained. Original 91-pair/per-case signed plots retain their documented
  two-slope normalization; the newer focal/comparison plots retain their
  symmetric signed scales. The comparison groups still share limits across
  scopes for a given metric/output.

## Completion and validation

- `plan.json`: exact original figure inventory, all 3,042 backup identities,
  2,107 numerical-file hashes, source tables and renderer identity.
- `original_figures/`: complete pre-axis-change image backup. Files were linked
  before rendering; all new images are written by atomic replacement, preserving
  the original linked file contents.
- `source_before/`: the previous local renderer.
- `documents_before/`: previous comparison README, interpretation and receipts.
- `receipts/`: one completion record per figure; each certifies both PDF/PNG.
- `progress.json`: restartable dataset rendering progress.
- `completion.json` and `figure_manifest.tsv`: final output coverage, hashes,
  original-backup verification, unchanged numerical inputs and response surfaces.
- `main_pdf_validation.json`: geometry plus the 15 main PDF checks.
- `pdf_validation.json`: geometry plus all 1,521 PDF text/page-boundary checks.
- `visual_review.json`: selected rendered PDF pages inspected visually.

PDF and PNG use one resolved layout before export. This fixes a PDF-first
constrained-layout display issue found during visual QA while retaining
embedded TrueType text. The diagnostic iteration is preserved under
`iterations/before_pdf_layout_freeze/`.

Historical collection and earlier redraw receipts remain as provenance. Their
old figure hashes describe the corresponding backed-up versions; the current
image hashes are in this directory's completion/figure manifests.

## Reproduce or resume

Requires Python, NumPy and Matplotlib. PDF validation also requires pdfplumber.
The local renderer uses NumPy 2.2.6 and Matplotlib 3.9.4.

```sh
python3 oxygen/code/O2_supply_demand_MAP/analysis/eFAST/sobol/redraw_full_n4096_oxygen_axis.py \
  --root oxygen/results/eFAST/sobol/full_N4096 --workers 8
python3 oxygen/code/O2_supply_demand_MAP/analysis/eFAST/sobol/validate_full_n4096_oxygen_axis.py \
  --root oxygen/results/eFAST/sobol/full_N4096 --workers 4
```

Phases are `prepare`, `main`, `cases`, `finish`, or `all` (default). Completed
figures are verified and reused. Renderer identity changes require an explicit
`--refresh-renderer`; original backups are not replaced. No SSH, Slurm, SALib
analyzer, model computation or numerical source-file writes are performed.

## Verified outcome

Completed 2026-10-19: all 1,521 figures (3,042 PDF/PNG files) were regenerated.
All 2,107 numerical-file hashes and all original figure backups verified.
All 108 response-surface image files are unchanged. All 1,521 final PDFs passed
axis geometry, expected row-label, region-label collision and page-boundary
checks. The final 12-panel comparison, global 91-pair S2 and seed500 91-pair S2
PDF rasterizations were inspected visually. See `visual_review.json` for exact
reviewed identities. Model evaluations and sensitivity estimations: zero.
