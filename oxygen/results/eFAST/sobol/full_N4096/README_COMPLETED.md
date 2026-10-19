# Completed recovery (2026-10-18)

All 500 neighborhoods and the global domain are complete at N4096, five repeats
and 201 oxygen points, with all 91 parameter pairs retained. The 54 conditional
response surfaces and final report completed successfully. Convergence remains
a diagnostic; no endpoint was excluded.

- [Recovery jobs, scope and validation](recovery_20261018/retry5734/README.md)
- [Global figures](global/figures) and [global source summary](global/summaries/indices_scope_summary.tsv.gz)
- [Neighborhood figures](neighborhood10pct/figures) and [neighborhood source summary](neighborhood10pct/summaries/indices_scope_summary.tsv.gz)
- [Response surface inventory](response_surface_inventory.tsv)
- [Method and interpretation limits](analysis_notes.md)
- [Verified full collection receipt](artifact_collection_receipt.json)

The verified local/HPC collection contains 19,828 files (1,619,105,148 bytes).
Git contains primary figures, scope summaries, all 54 response figures/source
tables, reference parameter sets and audit receipts. Detailed per-endpoint
files can be retrieved with collect_full.py using collection_manifest.tsv.

## Local figure display update (2026-10-19)

All current oxygen-axis figures now use a base-10 symlog axis, retaining real
O2=0 with a linear interval from 0 to 0.025%. Dashed boundaries at 1% and 3%
mark Low/Mid/High oxygen. Numerical tables, repeats and convergence diagnostics
are unchanged. The 54 parameter response surfaces are unchanged.

- [Current figures and complete redraw audit](axis_redraw_20261019/README.md)
- [Global/neighborhood comparison figure groups](figure_comparison_20261019/README.md)

The earlier collection receipt and file hashes above describe the original
collected artifacts. Current image identities and verification are recorded in
`axis_redraw_20261019/completion.json` and `figure_manifest.tsv`.
