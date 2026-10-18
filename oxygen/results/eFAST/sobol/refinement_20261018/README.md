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
