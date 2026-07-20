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
| `migrate_existing_html_report_lightboxes.R` | Idempotent maintenance entrypoint that adds the shared image viewer to existing HTML reports and records a TSV migration manifest under `oxygen/results/analysis/`. |
| `run_provenance_report.R` | HTML-formatting helper for recorded commands, arguments, configuration snapshots, and provenance tables. |
| `fit_results/` | Extra-results and paired-sigma report assemblers. |
| `fixed_o2/` | Canonical fixed-O2 report assembler and report-only helpers. |
| `fixed_o2_eigen/` | Fixed-O2 eigen-attractor report. |
| `multi_warmup/` | Multi-warmup collected-results report. |
| `parameter_landscape/` | Parameter-landscape clustering and contribution reports. |
| `process_fingerprints/` | Process, ploidy-regime, medium-O2, and O2-ploidy coupling reports. |
| `joint_coupling/` | Static HTML assembly of the joint ClassA/B/C stability and CatA/B/C ploidy-association tables and figures. |

PDF-to-image conversion used only to embed an existing PDF preview is allowed
presentation staging; it is not analytical figure generation.

## Shared image viewing

Every canonical image-bearing HTML writer calls
`o2sd_inject_report_image_lightbox()` after it has completed its existing
report assembly. The injected display layer leaves the report body, navigation,
captions, links, and image payloads unchanged. Clicking an existing report image
opens it in the current page with button, wheel, keyboard, and drag controls for
zooming and panning. `Fit` never enlarges an image above its natural size; the
corrected zoom-out floor is `max(Fit scale × 0.35, 0.025)`.

To migrate already generated reports without rerunning analysis or fitting:

```bash
Rscript oxygen/code/O2_supply_demand_MAP/report/migrate_existing_html_report_lightboxes.R
```

The command scans `oxygen/results/**/*.html`, skips reports without images,
skips files already carrying the current viewer, and writes
`oxygen/results/analysis/report_image_lightbox_migration.tsv`.

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
