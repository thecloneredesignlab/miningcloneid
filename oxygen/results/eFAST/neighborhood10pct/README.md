# Per-endpoint in-vivo eFAST neighborhoods

This folder contains a separate sensitivity analysis around each of the 500
independent **in-vivo** optimizer endpoints. It does not use in-vitro or joint-fit bounds.

For parameter j and endpoint s, with the original natural bounds L and U:

```
lower(s,j) = max(L(j), best(s,j) - 0.10 * (U(j) - L(j)))
upper(s,j) = min(U(j), best(s,j) + 0.10 * (U(j) - L(j)))
```

Thus the untruncated total width is 20% of the original natural range.
Log10 parameters are log-uniform within these natural bounds; other parameters
are uniform. Parameters are varied independently and no local refitting is done.
Fourteen parameters affect the fixed-oxygen operator. The oxygen grid is the
same 201 points, 0–5% by 0.025%, as Figure 4. Each endpoint uses N=513, M=4
and **five independent phase repetitions**, regardless of pilot convergence.

## Completed checks on 2026-10-08

- All 500 endpoints and 7,000 neighborhood bounds passed the input audit.
  All 1,500 reference comparisons passed (maximum absolute error
  `6.6489036498751375e-12`). Six parameters per endpoint require clipping at the
  median; every clipping event is documented in `neighborhood_bounds.tsv`.
- All five regression checks passed locally and inside the specified SIF on
  hpctpa3pc0028. A separate synthetic end-to-end check verified median versus
  mean aggregation, labeled source arrays, and rendered PDF layout.
- The timing pilot completed for seed25 (best fit), seed464 (worst fit), and
  seed165 (most clipped), using N=129 and O2=0,0.5,2.5,5%, five phases each.
  All 15 designs and 108,360 operator evaluations passed validation.
- `pilot/runtime_projection.json` estimates 33.13 days (median throughput) to
  35.17 days (p90 operator time) for the full 500-endpoint, 16-worker analysis.
  This is an extrapolation from four oxygen points and excludes some analysis
  overhead; the full-grid pilot will refine it. The 90.7 GiB raw-output estimate
  comes from the existing global N=513 compressed outputs.
- `pilot/smoke_seed*_convergence.tsv.gz` retains per-endpoint five-phase means
  and ranges. `pilot/smoke_seed*_phase_indices_N129.npz` retains the exact
  indices with axes `(phase, index [S1,ST], parameter [ACTIVE order in efast.py],
  O2 [oxygen array], output [ploidy,growth])`.

At this dated checkpoint the background pipeline has entered the global
R3–R5 extension, followed by the ten-endpoint resolution pilot. The full
500-endpoint neighborhood calculation has **not** started and its convergence
has **not** been established. Controller PID: 1148713; log on HPC:
`oxygen/results/eFAST/neighborhood10pct/logs/pipeline_20261008.log`.

## Execution and diagnostic policy

On 2026-10-08 the user changed convergence from an execution gate to a
diagnostic. After the ongoing ten-endpoint pilot, the pipeline continues to
all 500 endpoints at N=513 with five phases regardless of the pilot result.
Unstable endpoints remain in the full analysis. The existing controller can
finish its loaded pilot code; its next, separate full-stage Python process
loads the updated diagnostic-only policy. `execution_policy.json` records
this policy, the input hash and the pilot result used at full-stage startup.

`run_neighborhood.sh pipeline` audits 500 endpoints and validates their original
operators at 0, 2.5 and 5% oxygen. It runs three representative endpoints at
N=129 and four oxygen points, five phases each, and records runtime projections.
It then extends the global designs to five phases and checks ten representative
endpoints at N=129,257,513 on the full oxygen grid, five phases each.
The pilot covers the best/worst fitting endpoints, clipped neighborhoods,
seed25, and farthest points in normalized transformed parameter space.

The pilot records whether every output/index/endpoint has all cells defined,
a p90 five-phase index range <=0.10, and a p90 absolute change from the N=257
mean to N=513 mean <=0.05. These are declared practical tolerances in
variance-fraction units, not mathematical guarantees. The historical filename
`pilot/convergence_gate.json` is retained, but `needs_review` no longer blocks
full execution. Full execution still requires valid audited inputs and at
least 200 GiB free disk space.

After each full endpoint finishes, the runner updates
`summaries/seed_convergence_diagnostics.tsv` (endpoint/output/index checks)
and `summaries/seed_convergence_status.tsv` (one status per completed endpoint).
`repeat_status` describes the five-phase range. `resolution_status` describes
the N=257-to-513 change only for pilot endpoints; the other endpoints are
explicitly `not_tested`. Overall `convergence_status` is `passed` only when
both checks pass, `not_passed` when a check fails or is undefined, and
`resolution_not_tested` when repeats pass but resolution was not tested.
These statuses never remove an endpoint from the analysis.

No Slurm submission is used. The runner enforces hpctpa3pc0028 and the validated
SALib 1.5.2 SIF checksum. Only one neighborhood controller may run per root;
`EFAST_WORKERS` controls fork workers and defaults to 16.

## Outputs and interpretation

- `seed_manifest.tsv`, `neighborhood_bounds.tsv`, `operator_input_audit.json`,
  and `input_manifest.json`: all endpoints, clipping, references and input hashes.
- `bounds/seed*.tsv`: per-endpoint bounds on HPC.
- `pilot/runtime_projection.json`, `pilot/convergence_diagnostics.tsv`,
  `pilot/convergence_gate.json`: measured costs and convergence checks.
- `execution_policy.json`: diagnostic-only full-run policy and input provenance.
- `summaries/seed_convergence_status.tsv`, `seed_convergence_diagnostics.tsv`:
  incrementally updated per-endpoint status and per-output/index diagnostics.
- `runs/seed*/runs/N513_R*/`: ordered samples, raw outputs, indices and completion
  receipts on HPC. A complete design is reused only after its hashes match.
  Partial designs are recomputed; completed designs and completed endpoint
  summaries are retained. Evaluation is streamed in chunks to bound memory.
- `summaries/phase_indices_*.npz`: exact phase indices in 100-endpoint shards with labeled axes
  `(fit_seed, phase, index, parameter, O2, output)`; five phases stay together.
- `summaries/seed_mean_indices_*.tsv.gz`: per-endpoint five-phase means, SD,
  ranges, and valid-phase counts. No pooling of raw outputs across endpoints.
- `convergence.tsv`: across-endpoint mean, median, Q25, Q75 and valid-seed counts.
- `summaries/seed_parameter_bands.tsv.gz`, `rank_consistency.tsv`,
  `mechanism_band_summary.tsv`, `mechanism_support.tsv`: window scores, ranks,
  group means and fraction of evaluable endpoints supporting each mechanism claim.
- `figures/efast_four_panel.pdf/png`: medians across endpoints, separate S1/ST
  panels, separate ploidy/growth scales, three row annotations. The original
  Figure 4B association group remains; the new sensitivity group follows ploidy
  S1 or ST. Complete endpoint curves are bootstrap units after phase averaging.
- `raw_output_inventory.tsv`, `collection_manifest.tsv`, `interpretation.md`:
  raw-output hashes, collected artifact hashes and interpretation.

Raw neighborhood samples/outputs and per-run caches stay on HPC; compact source
tables, labeled index arrays, figures and provenance are tracked in Git.
Projected raw compressed outputs are about 90.7 GiB before pilot data. The full
design requires 3,608,955,000 operator evaluations. Runtime projections are
empirical estimates and can change with node load and parameter values.

Indices from numerically constant trajectories are undefined and recorded as
NaN with variance and valid-count diagnostics. No failed sample row is removed.
Optimizer endpoints are not independent posterior draws. Bootstrap groups and
support fractions describe these fitted endpoints, not posterior probabilities.
Total-effect group means overlap through interactions and are not additive
group variance fractions. The existing Spearman panel retains direction.
