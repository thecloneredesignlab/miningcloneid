# Report layer

`report/` is a consume-only presentation layer. It reads fitted summary tables,
materialized analysis tables, visualization manifests, and existing figures,
then assembles HTML/PDF reports. It does not source model/simulation code, read
`fit_result.rds`, recompute fitted quantities, or create analytical plots.

The upstream order is:

```text
fit outputs -> simulation -> analysis -> vis -> report
```

## Report families

| Path | Responsibility |
|---|---|
| `render_fit_report.R` and `fit_result_report.Rmd` | General in-vivo/joint seed report assembly from existing tables and figures. |
| `render_invitro_fit_report.R` | In-vitro seed report assembly. |
| `render_fixo2_invivo_report.R` | Stable compatibility entry for fixed-O2 report assembly. |
| `run_provenance_report.R` | HTML-formatting helper for recorded commands, arguments, configuration snapshots, and provenance tables. |
| `fit_results/` | Extra-results and paired-sigma report assemblers. |
| `fixed_o2/` | Canonical fixed-O2 report assembler and report-only helpers. |
| `fixed_o2_eigen/` | Fixed-O2 eigen-attractor report. |
| `multi_warmup/` | Multi-warmup collected-results report. |
| `parameter_landscape/` | Parameter-landscape clustering and contribution reports. |
| `process_fingerprints/` | Process, ploidy-regime, medium-O2, and O2-ploidy coupling reports. |

PDF-to-image conversion used only to embed an existing PDF preview is allowed
presentation staging; it is not analytical figure generation.

## Consume-only regression

`oxygen/tests/testthat/test-report-consume-only-boundary.R` renders a report
from materialized tables and figures without any fit RDS object and checks that
upstream input hashes do not change. The global architecture check also rejects
report code that sources upstream executable layers or draws plots:

```bash
Rscript oxygen/tests/check_o2_layer_boundaries.R
```

Concrete purpose, required inputs, outputs, functions, and direct tests for
every report file are listed in `../docs/CODE_FILE_REGISTRY.md` and
`../docs/code_file_registry.tsv`.
