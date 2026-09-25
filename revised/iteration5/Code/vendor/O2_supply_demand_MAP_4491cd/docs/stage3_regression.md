# Stage 3 regression: analysis-only responsibilities

Stage 3 moved numerical reconstruction and plotting out of `analysis/`.
Canonical analysis entrypoints now consume completed fit summaries or
materialized simulation tables and write analytical/plot-ready tables only.
Historical public paths are retained as explicitly labelled compatibility
orchestrators where external callers require them.

## Fit-result workflows

The extra-results workflow was split into fitted-parameter prediction
materialization, cross-seed analysis, visualization, report assembly, and one
runner. The same pattern was applied to paired joint-sigma, sigma-burden, and
long-ploidy seed-selection workflows.

Validated evidence:

- Existing ten-seed extra-results fixtures reproduced all nine core TSV
  contracts, including the 10-row by 329-column seed summary, the 210-row by
  24-column parameter table, and four prediction tables.
- Paired joint-sigma analysis reproduced four reference tables.
- Sigma-burden comparison reproduced three reference tables.
- Long-ploidy materialization and selected-seed outputs matched their legacy
  fixture.
- `test-fit-results-stage-split.R` passed 98 expectations at the time of the
  component gate.

## Process, ploidy, O2, and coupling analyses

Process-fingerprint, ploidy-regime, medium-O2-window, and O2-ploidy
event-coupling workflows were separated into numerical producers under
`simulation/process_fingerprints/`, table-only analyses under
`analysis/process_fingerprints/`, pure figure consumers under
`vis/process_fingerprints/`, and consume-only report assemblers under
`report/process_fingerprints/`.

Validated evidence:

- All four historical entrypoints completed against mock/full fixtures through
  their compatibility orchestrators.
- `test-invivo-process-analysis.R` passed after the split.
- `test-process-fingerprint-stage-split.R` passed its layer, path, and
  source-without-main checks, including a dynamic assertion that simulation and
  analysis helpers do not create figure or report directories.

## Live-cell effective missegregation analysis

The numerical live-effective-`p_ms` producer now lives at
`simulation/invivo/cin/generate_live_effective_pms_outputs.R`. The comparison
under `analysis/profile_likelihood/` consumes its manifests and tables and the
figure is produced only under `vis/profile_likelihood/`.

For real in-vivo seed366:

- five scientific output tables matched the legacy path at tolerance `1e-12`;
- all non-path context fields matched exactly;
- the deprecated estimator path reproduced the same outputs;
- by-seed scientific columns and the comparison summary matched the frozen
  fixture;
- the first-page figure raster SHA-256 was
  `af9e15dac6d89eb30cb7475bf75823caf276980eb533a5ba0e5816d595daa4bb`;
- `test-profile-likelihood-stage-split.R` passed 34 expectations.

## Stage status

Component-level real-data and boundary regressions above passed. The final
aggregate architecture check, complete unit suite, immutable-core check, and
repository-wide syntax/whitespace gate are intentionally recorded as pending
until all five stages stop changing:

```bash
Rscript oxygen/tests/run_o2_reorganization_regression.R
```

No aggregate pass should be inferred from this document until that final
command is recorded as successful in the handoff.
