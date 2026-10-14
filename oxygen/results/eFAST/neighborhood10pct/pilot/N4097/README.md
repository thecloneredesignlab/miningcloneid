# Six-endpoint N4097 eFAST pilot

Approved scope: `seed25, seed101, seed138, seed174, seed66, seed311`, N=4097,
M=4, five independent phase repeats, all 14 sampled parameters, and O2
0:0.025:5% (201 points). This is a diagnostic subset of the original 500
in-vivo fitted endpoints. The neighborhood is each fitted optimum plus/minus
10% of the full natural parameter span, clipped to global bounds, retaining
the original uniform/log-uniform transforms.

## Status and execution

**Submitted on 2026-10-14.** Scheduler audit at 12:25:36 EDT verified:

| Role | Job ID | Tasks | Dependency |
|---|---|---|---|
| SIF smoke | 20863802 | 1 | None |
| N4097 trajectories | 20864027 | 420 | afterok:20863802 |
| Seed summaries | 20864028 | 6 | afterany:20864027 |
| Comparison and figures | 20864029 | 1 | afterany:20864028 |

At that snapshot, all jobs were pending; smoke had not passed and numerical
analysis was incomplete. Smoke was waiting on `QOSMaxCpuPerUserLimit`.
`slurm/submission_validation.json` validates resources and dependencies,
not numerical completion. `collection_receipt.json` verifies ten collected
files, 18 frozen code files, 40 preceding source files, and 18 baseline arrays.
Refresh scheduler state before treating this historical snapshot as current.

The initial code commit prepares this directory; submission is established
only by `slurm/submission_jobs.tsv` and live Slurm accounting. A passed
`smoke/smoke_receipt.json` establishes the compute-node container preflight.
`completion.json` is written only after all 420 trajectories, six summaries,
and final tables/figures complete. Numerical recovery records are retained
in the summary receipts. Convergence failures do not exclude an endpoint.

Smoke job **20863802** was queued after source commit
`bb4f60998134ace2a3d150257030521a350a87bb` reached HPC. It initially waited
with `QOSMaxCpuPerUserLimit`. `queue_resolution_4097_after_smoke.sh` allows
the 420 main trajectories to queue with an `afterok` dependency on smoke;
they run only after smoke succeeds. It validates the frozen plan, comparison
sources, image identity and smoke job state. Its source hash and submission
commit are recorded in `slurm/deferred_submission.json`. Both submitters
share a lock and duplicate-submission marker. This changes queue timing;
the frozen numerical sources and mandatory compute-node preflight stay intact.

- SIF: `/share/lab_crd/taoli/Docker/o2_supply_demand_map_r442_hpc_exact_salib152_20260923.sif`.
- Runtime: SALib 1.5.2 in the existing Python environment.
- Main array: 420 tasks = 6 endpoints x 5 phases x 14 complete focal trajectories.
- Main task: 1 CPU, 4G; summary: 1 CPU, 8G; final comparison: 4 CPU, 32G.
- All submitted jobs: `qos=xxlarge`, `time=12:00:00`, no node restriction,
  no array concurrency suffix.
- Design evaluations: 345,868,740 parameter-vector/oxygen rows, excluding
  reference validation and smoke checks. Both outputs share the same solve.
- Raw sampled designs and per-trajectory outputs stay on HPC in `runs/`.
  This repository collects compact manifests, receipts, indices and figures.

## Inputs and comparisons

`slurm/submission_plan.json` records source/input/image hashes, Git commit,
resources and selection. `seed_selection.tsv` documents each endpoint.
`baseline_indices.npz` freezes the selected N513/N1025/N2049 five-phase
indices from the completed preceding pilot. Source files and their original
receipts are hashed before freezing; finalization rechecks those sources.
The seed-specific five phase seeds are inherited across resolutions. The
FAST designs at different N are not nested samples.

## Primary question: ploidy ST convergence

All 14 parameters are checked. Focused plots and tables emphasize
`buffer_beta`, `buffer_smax`, `p_mis_base`, `mu_hp`, `O2_crit`.

1. N4097 repeat range: `max(ST across five phases) - min(ST across five phases)`.
2. N2049 to N4097 difference: absolute difference of the two five-phase means.
3. Calculate p90 and maximum separately within overall [0,5], low [0,1],
   and high [3,5]% O2; save the O2 position of each maximum.
4. Repeat-range p90 threshold: 0.10. Mean-change p90 threshold: 0.05.
   These are absolute index units, not relative percentages.
5. `passed_by_threshold`: p90 and maximum pass. `localized_review`: p90
   passes, maximum fails. `not_passed`: p90 fails. Missing/undefined data
   remain `undefined_or_incomplete`, never zero.

## Outputs after completion

- `summaries/seed*_N4097.npz`: five-phase S1/ST and output variances.
- `repeat_diagnostics.tsv`, `resolution_diagnostics.tsv`: all inputs/outputs.
- `key_ploidy_ST_*_diagnostics.tsv`: five focused parameters.
- `oxygen_curve_summary.tsv.gz`: means, SD, full phase ranges and counts.
- `band_rank_mechanism_*.tsv`: low/high contrasts, parameter ranks and
  descriptive mechanism checks. Process means are not formal group ST.
- `figures/ploidy_ST_resolution_comparison.pdf` and growth companion:
  N513/1025/2049/4097 means with full five-phase min/max ribbons (not CIs).
- `figures/ploidy_ST_convergence_overview.pdf`: error metrics for N2049,
  N4097 and their change; cells show diagnostic errors, not sensitivity.

Selected endpoint pass fractions cannot be extrapolated to all 500 seeds.
Any recommendation for N8193 or a separate frequency-setting comparison
requires reviewing these results first; neither is part of this submission.
