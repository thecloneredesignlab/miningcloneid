# Stage 2 regression: simulation, analysis, and visualization split

Stage 2 replaced visualization-time numerical reconstruction with explicit
post-fit simulation and analysis products. Tests were run against existing
completed fits; no fitting was repeated.

## In-vivo seed366

- Legacy output: `/tmp/o2_invivo_seed366_legacy.z733qC`
- Simulation output: `/tmp/o2_invivo_producer_full.wb5MIp`
- Pure visualization output: `/tmp/o2_invivo_consumer_full.9tud1q`
- Scientific tables: 31/31 legacy TSV files have identical names, rows,
  columns, text, and numeric values; maximum absolute difference is zero.
- Largest table checked: 1,065,064 rows.
- Figures: 28/28 expected names are present. All 4 PNG files are byte-identical.
  At 96 dpi, 20/24 PDF first pages are pixel-identical. The remaining four
  differ only because the corrected manifest supplies `fit_dir=seed366`
  instead of the legacy temporary-directory label `fit_dir=tmp`.

## In-vitro seed10

- Simulation output:
  `/tmp/oxygen_stage2b_invitro/simulation_seed10_domain_modules`
- Analysis output:
  `/tmp/oxygen_stage2b_invitro/analysis_seed10_domain_modules`
- Pure visualization output:
  `/tmp/oxygen_stage2b_invitro/stage2_viz_domain_modules_seed5826`
- Scientific tables: 10/10 tables agree at tolerance `1e-12`.
- Fixed-seed figures: 11/11 PNG files are byte-identical. PDF differences are
  limited to generated document metadata.
- The visualization smoke fixture contains no `fit_result.rds` or
  `best_params.tsv` and still produces all 22 expected plot files plus its
  manifest.
- The diagnostics condition proxy differs by approximately `4e-14` relative
  after a TSV write/read round trip; all other diagnostics tables agree at
  tolerance `1e-12`.

## Perturbation workflows

- Mixed-ploidy simulation: all five legacy scientific files, including the
  compressed trajectory table, are byte-identical.
- Factorial interaction simulation: all four untreated scientific files are
  byte-identical.
- The mixed-ploidy pure visualization produces four figures plus a manifest.
- The factorial pure visualization produces 59 artifacts plus a manifest.
- Compatibility wrappers were exercised with figures enabled and disabled,
  from the repository root and from an unrelated working directory.

## Joint post-fit products

- The nine in-vivo versus in-vitro comparison PNG files are byte-identical to
  the former report-embedded implementation.
- Joint parameter analysis table SHA-256:
  `71aaedc6eb948d8abb19d55fc4fad1225548626fc354961a2f49e4b8aa777329`.
- Joint parameter PNG SHA-256:
  `4499746e19111a97afd7889f8e1b0f3afba1e624b5fa30fbeb94598deff00291`.
- The same hashes were reproduced after the Stage 5 util cleanup.

## Report consumer boundary

A real joint seed3 report was rendered from a temporary fixture containing 142
materialized inputs and no `.rds` files:

`/tmp/o2_report_consume_smoke.seed3.IdWsye/seed3/report_smoke/consume_only_seed3.html`

The HTML is 7,436,982 bytes and contains 34 figure cards. All input SHA-256
values were unchanged after rendering.

## Test gates

- `test-invivo-stage2-boundaries.R`
- `test-invitro-stage2-boundaries.R`
- `test-perturbation-stage-split.R`
- `test-stage2-postfit-orchestration.R`
- `test-report-consume-only-boundary.R`
- `test-stage5-util-consolidation.R`
- full `oxygen/tests/run_unit_tests.R`
- `Rscript oxygen/tests/check_immutable_o2_core.R --full`

The protected `model/` and `optimizer/` tree passed the complete 18-file
size/mtime/SHA-256 check after every real-model smoke.
