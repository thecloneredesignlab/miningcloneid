# Four-endpoint eFAST follow-up: separate N and M experiments

Approved scope (2026-10-18): current branch; SALib 1.5.2 in the existing SIF;
xxlarge, 12 hours, no requested compute node, no array concurrency throttle.

| Configuration | Endpoints | Repeats | Main trajectories | Operator evaluations |
|---|---|---:|---:|---:|
| N8193_M4 | seed25, seed138 | 5 | 140 | 230,551,020 |
| N4097_M8 | seed101, seed174 | 5 | 140 | 115,289,580 |

Each trajectory retains all 14 jointly varied parameters, one focal FAST
frequency assignment, and all 201 oxygen points (0 to 5%, step 0.025%).
The fitted optimum +/-10% of the full natural parameter span, clipped bounds,
original uniform/log-uniform transforms, five phase seeds, canonical fixed-O2
operator and certified targeted Perron recovery are unchanged. Both outputs
share the operator solve. Counts exclude smoke/reference/recovery checks.

## Comparison contract

- N8193_M4: primary N4097_M4 -> N8193_M4 resolution comparison for seeds 25/138.
- N4097_M8: primary N4097_M4 -> N4097_M8 frequency/harmonic-setting robustness
  comparison for seeds 101/174. This is not an N convergence step and does not
  imply that M8 is intrinsically more accurate.
- N2049_M4 provides historical context in both groups. M4 and M8 phases are
  never pooled into ten repeats. Equal phase seeds do not make designs nested.
- Both N and M enter preparation signatures, task metadata, paths, summaries,
  FFT frequency selection and pre/post-recovery index-impact diagnostics.
  Main omega = floor((N-1)/(2*M)); complementary frequency budget is 128 for
  N8193_M4 and 16 for N4097_M8 (enough for 13 complementary parameters).
- Practical diagnostics: five-phase range R90 <= 0.10; absolute difference of
  phase means D90 <= 0.05, each across an O2 window. Report maxima and their
  locations as well as p90 for the full [0,5], low [0,1], and high [3,5] grids.
  M-comparison D90 quantifies setting dependence, not N convergence.
- Keep all 14 parameters and both S1/ST outputs. Highlight ploidy ST for
  buffer_beta, buffer_smax, p_mis_base, mu_hp and O2_crit. Do not omit a seed
  because it fails convergence. These selected diagnostic seeds do not
  represent the distribution of all 500 fitted endpoints.

## Reproduction and scheduling

Code: `oxygen/code/O2_supply_demand_MAP/analysis/eFAST/resolution_nm.py`,
`test_resolution_nm.py`, `submit_resolution_nm.sh`, `run_resolution_nm.sbatch`.

After committing/pushing code and synchronizing the same commit on HPC:

```bash
bash oxygen/code/O2_supply_demand_MAP/analysis/eFAST/submit_resolution_nm.sh prepare
bash oxygen/code/O2_supply_demand_MAP/analysis/eFAST/submit_resolution_nm.sh submit
```

For each configuration the submitter queues one smoke job, a 1-140 main
trajectory array depending on `afterok:smoke`, a 1-2 seed-summary array, and one
final comparison job. All use `--qos=xxlarge --time=12:00:00`, without a node
request or `%` concurrency suffix. Summary/final jobs use `afterany` so that
missing or failed inputs are reported; their validators require all expected
inputs before declaring completion. A main worker independently requires a
passed smoke receipt bound to its exact plan. Submission is idempotent using
saved job IDs and a submission lock. An absent/failed dependency requires
inspection; the runner does not silently create another attempt.

The compute-node smoke runs every new test, evaluates a full-N trajectory at
0/2.5/5% O2 for both endpoints in its configuration, and checks the unchanged
50/100-digit recovery solver on an actual sampled canonical matrix. A failed
smoke prevents its main array from starting. Main tasks use 1 CPU/4 GiB,
seed summaries 1 CPU/8 GiB, final comparisons 4 CPUs/32 GiB, smoke 1 CPU/4 GiB.

Image:
`/share/lab_crd/taoli/Docker/o2_supply_demand_map_r442_hpc_exact_salib152_20260923.sif`

SHA256: `0b60f6cdab8a91f6bbeea1ab6cf02dd79660a295f24bc6e8ad6c5fed04d982ad`.

Each configuration has its own input/baseline hashes, task and summary
manifests, submission receipts, smoke outputs, trajectory outputs, five-phase
NPZ summaries, repeat/comparison tables, O2 profiles, and PDF/PNG ploidy-ST and
growth-ST comparison figures. The archived preceding pilots remain intact.

## Status

Local configuration, cache identity and submission tests passed (5); four
SALib-dependent tests are reserved for the approved SIF. All 18 frozen
scientific source hashes and Python 3.9 syntax match. A submission receipt is not
scientific completion. Check smoke, scheduler terminal states/exit codes,
trajectory receipts, summaries and completion.json separately.
