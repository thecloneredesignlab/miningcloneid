# Incremental Sobol refinement approved 2026-10-18

## Scientific scope

- Global: all 14 factors, all 91 pairs, all 201 oxygen points; N8192,
  N16384 and N32768, retaining the original five scramble seeds.
- Neighborhoods: 497 seeds, exactly 55,089 selected seed/oxygen combinations,
  N8192 and five repeats. A target is the union of nonpassing S2 diagnostics
  for six pairs and either output. Individual estimates are never selected by
  their sign or magnitude. The pair list was fixed using a scope-median peak
  threshold of 0.01, which is a budget rule, not a significance test.
- All original N4096 results remain read-only and all original 500 endpoints
  remain in the ensemble. Convergence never excludes an endpoint.

## Exact reuse

The immutable `selection.json` contains source-table checksums and case masks;
`targets.tsv.gz` records each target and the diagnostics that triggered it.
Preparation in the original SALib 1.5.2 SIF proves all five N4096 unit-design
prefixes bit-for-bit. Case audits reconstruct original natural-space sample
identities and validate their 150 completed model blocks. Production model
tasks evaluate only base samples 4097 onward. Statistical processing joins
the old and new ordered outputs before calling the unchanged SALib analyzer.

An actual-model smoke task writes the first eight new oxygen slices into the
production checkpoint. Its results are reused by the first production task;
the smoke test never evaluates an old sample. Synthetic unit tests also check
prefix nesting, block interleaving, no-tail overlap and estimator parity.

## Storage and recovery

Only new model tails are stored. Full old/new concatenations are assembled in
memory during analysis and never duplicated on shared storage. Preparation
requires 650 GiB free; workers preserve a 64 GiB reserve. Model outputs and
solver flags retain lossless precision. Model checkpoints certify small oxygen
chunks; analysis checkpoints certify every completed repeat/oxygen/output cell.

## Outputs

- `../global_refinement/N8192`, `N16384`, `N32768`: global summaries and plots.
- `../neighborhood_refinement/N8192`: matched target statistics, source tables,
  targeted-only plots, and a 500-endpoint mixed-resolution summary.
- `N_used.tsv.gz` explicitly identifies N4096 versus N8192 cells. Mixed plots
  must never be described as all 500 neighborhoods uniformly sampled at N8192.
- Existing conditional response surfaces are reused; their model values do
  not change with Sobol N.

## Submission

Code: `oxygen/code/O2_supply_demand_MAP/analysis/eFAST/sobol/refinement_run.py`,
`submit_refinement.py`, `run_refinement.sbatch`, `test_refinement.py`.

Host submission (RED Python 3.6; use the committed code revision):

```sh
python3 oxygen/code/O2_supply_demand_MAP/analysis/eFAST/sobol/submit_refinement.py --revision CODE_COMMIT
```

Every analysis job requests xxlarge, one CPU and an initial 12 hours, without
a specified node or array concurrency cap. Each parent job has a separately
registered monitor using the existing RED home-directory helpers: add six
hours when less than six hours remain, polling every 300 seconds until the
target ends. The monitor itself uses the previously established small QoS.

`slurm/jobs.json`, `launch_receipt.json`, and `time_extension/submission_jobs.tsv`
are the authoritative submission records. A submitted array is not completed
computation. Audits, smoke, model, statistics and reports have separate states.

## Submitted jobs (2026-10-18)

| Stage | Job | Extension monitor |
|---|---:|---:|
| prepare | 21140536 | 21140550 |
| audit_global | 21140537 | 21140551 |
| audit_neighborhood | 21140538 | 21140552 |
| smoke | 21140539 | 21140553 |
| model_g8192 | 21140540 | 21140554 |
| analyze_g8192 | 21140541 | 21140555 |
| model_g16384 | 21140542 | 21140556 |
| analyze_g16384 | 21140543 | 21140557 |
| model_g32768 | 21140544 | 21140558 |
| analyze_g32768 | 21140545 | 21140559 |
| report_global | 21140546 | 21140560 |
| model_n8192 | 21140547 | 21140561 |
| analyze_n8192 | 21140548 | 21140562 |
| report_neighborhood | 21140549 | 21140563 |

Code revision: `0f7f3a488b99722a98afec4717b7f50ab95409be`. All 14 parent jobs have a registered
six-hour-threshold / six-hour-increment monitor. Initial live verification
found preparation running and downstream stages pending on dependencies.
All 14 monitors were RUNNING. This is a submission receipt, not completion.

## Initial validation

As of 2026-10-18T21:50:40.531226Z: preparation and all 501 case audits completed
with exit code 0:0; all 11 tests passed in the original SIF; all five N4096
prefixes matched exactly. The actual-model smoke task was RUNNING and its
new-sample outputs will be reused. All 14 extension monitors were RUNNING;
no extension had yet been needed and no failure receipt was present.

Static publication receipts: `submission_receipt.json` and
`submission_verification.json`. Live runtime manifests and events remain on HPC.
