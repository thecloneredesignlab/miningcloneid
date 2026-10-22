# Figure7 / Figure8: canonical code and retained-data layout

All entry points now require the private 12-file locked model at
`Code/vendor/figure78_model` for source authentication. The initializer definition
is authenticated from that private copy, not loaded from a frozen data snapshot;
no numerical model is executed by drawing. See [entry guide](Figure78_entrypoints_EN.md).

Current scope note (2026-10-17): this document covers drawing only. Complete simulation/data generation is now provided by `run_Figure7.R` / `run_Figure8.R`; `data_Figure7.R` / `data_Figure8.R` call that real pipeline, not a validate-only stub. See [complete guide](Figure78_pipeline.md).

Each figure now has one plotting implementation and one data root. Dates and Vxx run directories
are not used to select results. Rendering never fits, simulates, integrates, re-aggregates seeds,
reads archives, or restores old directories.

## Commands

Run from revised/iteration5; this Mac's explicit Rscript path is /usr/local/bin/Rscript.

```sh
Rscript --vanilla Code/Figures/draw_Figure7.R
Rscript --vanilla Code/Figures/draw_Figure7.R --analysis=history
Rscript --vanilla Code/Figures/draw_Figure7.R --analysis=all
Rscript --vanilla Code/Figures/draw_Figure8.R
Rscript --vanilla Code/Figures/draw_Figure78_all.R --output-root=/absolute/path/to/trial
Rscript --vanilla Code/Figures/draw_Figure7.R --preflight-only=yes
Rscript --vanilla Code/Figures/draw_Figure8.R --preflight-only=yes
```

Figure7 defaults to eight instantaneous products; history selects eight lineage products and all
selects both analyses. Figure8 has ten products. The complete set contains 26 canonical products
and 54 PDF pages (four matrix supplements have eight pages each). Matrix component files remain.
Published history filenames retain figure7_ver2 / supp_fig7_*_ver2 for manuscript compatibility;
these identify a distinct scientific analysis, not an alternative software/data release.
Use --analysis; the old --version option is only a compatibility alias.
Individual supplement and old group scripts are thin wrappers over the same implementation.
Use --products=assembled_fig7 or --products=assembled_fig8 for a single main figure.

Use the drawing entries with --preflight-only=yes for validation without simulation.
The data_Figure7.R / data_Figure8.R entries now run the complete current forward pipeline;
do not invoke them as retained-data validators. See Figure78_pipeline.md for explicit stages.
Old numerical sources/configurations authenticate preserved published results only.
Newly generated full datasets use their own current provenance, without migration/frozen-worker manifests.

## Layout

- data/Figures/Figure7/tables: instantaneous, lineage, and diagnostic tables.
- data/Figures/Figure7/provenance: original bundles, receipts, frozen sources/configurations.
- data/Figures/Figure8/tables: current main and supplementary summaries.
- data/Figures/Figure8/seed_data: in vivo / in vitro seed snapshots.
- data/Figures/Figure8/provenance: parameters, passage rules, bundles, receipts, frozen sources.
- Code/Figures/util/analysis/figure7_plotting.R and figure8_plotting.R: one implementation each.
- figure78_plot_common.R, render_only.R, render_plots.R, render_io.R, provenance.R: shared plumbing.
- Code/config/figure7.json and figure8.json: version-free current configurations.

All scientific definitions, parameter values, finite-time/passage protocols, quantiles, filtering,
and visual style are unchanged. Instantaneous flux/growth share the current summary table;
the original flux reference is retained only for provenance authentication. Fixed-p mixture
tables use pooled_fixed_low/high names rather than outdated supplement numbers.

## Provenance and migration

Under Code/config/manifests, figure78_current_plot_inputs.tsv pins semantic groups, logical names,
relative paths, sizes and SHA-256. figure78_data_migration.tsv maps the original 56 unchanged
files. figure78_frozen_sources.tsv maps original sources to byte-identical frozen copies.
figure78_bundle_bindings.tsv preserves ordered original RDS code/model fingerprint bindings.
figure78_layout.json interprets original receipt paths; it is never a filesystem fallback.
figure78_plot_contract.tsv consolidates the scientific roles, presentation, uncertainty and page
counts of all 26 current products.

Original source/code receipts, RDS parameters and fingerprints are never rewritten or re-signed.
Frozen sources are authenticated, not imported as plotting code. Original model and Gaussian
initializer identities are checked against the RDS. Rendering needs no old version directories,
archive restoration, fitted input files, or external model source installation. Preflight verifies
every retained file for each selected figure and fails on missing/corrupt files without refitting.

## Outputs and validation

Default outputs are iteration5/Figures and manuscript/Figures. Use --output-root for isolation.
Rendering stages outputs, protects non-target images, backs up selected targets, and publishes
only after checks pass. Audits/caches/backups go under output-root/audit/figure78_render_only,
never the scientific data directories.

Options: --preflight-only=yes, --list-products=yes, --formats=pdf,png, --data-root, --audit-root.
Unknown options/products, conflicting analysis selections, and archive paths/symlinks are rejected.
Matrix PDF packaging requires pypdf; FIGURE78_PYTHON can select an installed Python. No installation
or HPC submission is automatic. validate_Figure78_render_only.py exercises real grouped/individual
commands, output coverage, PDF pages, PNG identity, mirrors, input failures, and non-target guards.
Migration evidence is under audit/figure78_canonical_consolidation_20261016.
