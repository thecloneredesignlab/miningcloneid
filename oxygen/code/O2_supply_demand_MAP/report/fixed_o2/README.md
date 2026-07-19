# Fixed-O2 report assembly

This directory assembles narrative artifacts from materialized tables and
figures. It does not load fitted objects, source the model, run simulations, or
draw scientific figures.

## Files

- `render_fixed_o2_report.R`: assembly-only CLI. It accepts a unified
  `--results_dir`, calls the existing HTML renderer, and writes
  `OUT_DIR/REPORT_BASENAME.html`.
- `fixed_o2_report_utils.R`: compatibility report helpers extracted from the
  former fixed-O2 monolith. They write section summaries from already computed
  tables and invoke the standalone HTML assembler used by the legacy runner.

The canonical HTML implementation remains
`report/render_fixo2_invivo_report.R`; this staged wrapper gives it an explicit
report-only entrypoint.
