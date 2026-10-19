# Local N4096 figure redraw and global/neighborhood comparison

Current summary display: [three annotations and oxygen-dependent row order](../row_annotations_20261019/README.md). This supersedes the earlier summary-figure display; individual seed/case figures retain the symlog overlay.

This directory is a local figure overlay dated 2026-10-19. It uses the two
existing scope summary tables without evaluating the model or re-estimating
any sensitivity index. The original five repeats, all 500 neighborhoods and
the 201 oxygen points remain represented. No convergence exclusions apply.

## Figures

- `global_neighborhood_s1_st_focal_s2.pdf`: the complete 12-panel group.
  Columns: global ploidy, neighborhood ploidy, global growth, neighborhood growth.
  Rows: S1, ST, and the 36 unique pairs involving p_mis_base, mu_hp or gamma_mu.
- `global_neighborhood_ploidy_s1_st_focal_s2.pdf`: the six-panel ploidy comparison.
- `global_neighborhood_growth_s1_st_focal_s2.pdf`: the six-panel growth comparison.
- `global_s1_st_focal_s2.pdf` and `neighborhood10pct_s1_st_focal_s2.pdf`:
  separate six-panel scope groups, with the same limits as the comparison groups.
- `global_neighborhood_diagnostic_fraction.pdf`: diagnostic coverage for every
  displayed cell. Global is 0/1; neighborhoods show the fraction of all 500 seeds.
- All figures also have PNG versions.

The original scope figure paths for S1 and ST are redrawn with the same positive
white-to-purple palette. Original files are preserved under `original_figures/`.
The original 91-pair S2 estimates are unchanged. New focal-pair PDFs are saved as
`../SCOPE/figures/s2_focal_pairs_by_oxygen.pdf`.

## Reading and scale policy

- Global cells are means of five repeats in one full-bound domain.
- Neighborhood cells are medians of 500 seed-specific five-repeat means, not
  results from one pooled 500-seed distribution and not averages of all runs.
- Positive sensitivity estimates are purple. Negative estimates are blue;
  undefined values are gray. Negative estimates do not represent antagonism.
- New signed color scales are symmetric around zero. The same metric/output
  has the same scale across scopes in every assembled group. Different metrics
  and outputs have separate color bars. Standalone scope panels retain their
  own ranges to show within-scope detail; compare absolute values using the
  assembled groups, not standalone color darkness.
- Cross-scope groups use neighborhood ploidy row order for each index. Standalone scope figures use their own ploidy row order; growth follows that order. S1/ST show three left annotations: Figure 4 process, original Figure 4B O2 correlation category and current ploidy-based Sobol O2 preference. S2 shows only process and O2 correlation. Rows with High >= Low descend by High-window score, followed by rows with Low > High ascending by Low-window score; there is no unresolved/middle sorting block. The third annotation is shown only for S1/ST and is not used to sort.
- Focal pairs remain 36 unique pairs involving p_mis_base, mu_hp or gamma_mu, but are sorted by oxygen preference rather than fixed focal-parameter blocks. The two halves in S2 annotations 1/2 follow the first/second parameter named in the row.
- The oxygen axis now uses base-10 symlog with a 0.025% linear threshold to
  retain O2=0. Major ticks are 0,0.025,0.1,0.5,1,3,5%. Black dashed lines at 1%
  and 3% separate Low [0,1]%, Mid (1,3)%, and High [3,5]% regions.
- Diagnostic thresholds are unchanged: five-repeat range <=0.05 and the
  absolute difference between N2048 and N4096 repeat means <=0.02.
  Passing these checks is neither a significance test nor proof of convergence.

## Source tables and provenance

- `plotted_indices.tsv.gz`: all displayed S1/ST/focal S2 values, original
  endpoint quantiles, endpoint/defined counts and convergence diagnostics.
- `focal_s2_source.tsv.gz`: the 36-pair subset, 28,944 scope/output/pair/O2 rows.
- `focal_pair_order.tsv`: canonical pair identities and the displayed orientation.
- `oxygen_band_comparison.tsv`: arithmetic oxygen averages of each scope
  summary in the documented bands. Neighborhood averages are averages of
  pointwise seed medians, not medians of seed-specific band averages.
- `focal_parameter_oxygen_points.tsv`: exact S1/ST values at 0,0.5,2.5,5% oxygen.
- `diagnostic_summary.tsv`: pass fractions and counts of negative scope summaries.
- `color_scales.tsv`: limits for every panel class.
- `redraw_manifest.json`: original/new figure hashes, numerical input hashes,
  software versions and the zero-model-evaluation receipt.
- `pdf_validation.json`: PDF boundary/text checks and the rendered pages reviewed.
- `interpretation_zh.md`: Chinese reading guide and comparison.

Collection-time manifests and the historical collection receipt are preserved.
Their old figure hashes describe the original collected figures, now backed up
here; `redraw_manifest.json` records the subsequent local overlay. Source-table
hashes are checked before and after rendering.

The later symlog redraw covers every oxygen-axis figure, including individual
seed plots. Current display hashes and validation are recorded in
`../axis_redraw_20261019/`; pre-axis-change comparison images and documents are
backed up there. The earlier redraw receipt remains historical provenance.

## Reproduce

Requires Python with NumPy and Matplotlib. This redraw used NumPy 2.2.6 and
Matplotlib 3.9.4; exact runtime metadata is in the receipt.

```sh
python3 oxygen/code/O2_supply_demand_MAP/analysis/eFAST/sobol/redraw_full_n4096_annotated.py \
  --root oxygen/results/eFAST/sobol/full_N4096
```

The script has no SSH, Slurm, SALib analyzer, model evaluator, or network calls.

## Chart contract

The selected surface is standalone scientific PDF/PNG output. Heatmaps compare
the same parameter identities over 201 ordered oxygen values. Paired scope
panels share scales for a fixed metric/output, while output-specific scales
remain separate. Purple/blue/neutral are the only data palette roots. All raw
estimates are retained; the companion diagnostic map supplies the separate
uncertainty information. Final QA uses rasterized PDF pages, PDF text checks,
pair/row completeness checks and unchanged numerical-source hashes.

Current classification scores and exact orders are in `../row_annotations_20261019/classification.tsv` and `row_order.tsv`. Classification uses normalized trapezoid windows; the earlier `oxygen_band_comparison.tsv` retains its original arithmetic summaries. Third-column labels have no 0.3 magnitude gate and do not indicate effect direction or proof of a nonzero contribution. Current image hashes and QA are in the annotation audit; earlier manifests remain historical.
