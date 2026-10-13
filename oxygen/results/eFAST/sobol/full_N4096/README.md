# Complete N4096 second-order Sobol run

Submitted 2026-10-13 on branch `soft_couping_eFAST`.
Code commit: `471f6a84afbf2363078a9cefce12a75cd0ee691f`.

## Scope

- One global independent in-vivo full-bound domain and every neighborhood seed1–seed500.
- Neighborhood radius ±10% of the full natural parameter span, clipped to global bounds;
  original uniform/log-uniform distributions retained.
- Fourteen active factors, all 91 pairs, N4096 and five independently scrambled repeats.
- Fixed oxygen 0–5% by 0.025%, all 201 Figure 4 reference points.
- Two outputs from each unchanged operator: dominant mean ploidy and leading-eigenvalue net live-cell growth.
- 75150 evaluation tasks, 61870694400 operator calls.
- No convergence gate or seed exclusion. N2048 prefixes provide diagnostics without additional operator calls.
- xxlarge,12h; no specified node or `%` array concurrency limit.

## Submitted pipeline

| Stage | Job | Tasks |
|---|---:|---:|
| Numerical tests and snapshot cache warmup | 20839420 | 1 |
| All-domain reference/configuration validation | 20839421 | 501 |
| Configuration and storage validation report | 20839422 | 1 |
| Complete operator evaluation | 20839423 | 75150 |
| Per-domain indices and diagnostics | 20839424 | 501 |
| Separate global and 500-neighborhood reports | 20839425 | 1 |
| Representative conditional response surfaces | 20839426 | 54 |
| Final figures, tables and inventory | 20839427 | 1 |

Scheduling dependencies and submission timestamps are in `slurm/submission_jobs.tsv`.
Bootstrap and the 501-domain validation concern correctness/storage, not convergence.
The later `slurm/scheduler_requeue_status.json` records 23 node-initialization
(`Prolog`) failures and their current `COMPLETING` state during Slurm automatic
requeue. Representative jobs have `Requeue=1` and `Restarts=1`; no duplicate
array was submitted. Accounting attempt failures are kept distinct from the
current queue state and from numerical convergence diagnostics.
Results are complete only after `status.json` says complete and the final report
job is COMPLETED with exit 0:0. This README records submission, not analysis completion.

## Storage and recovery

HPC root:
`/share/lab_crd/taoli/Project/soft_couping_eFAST/oxygen/results/eFAST/sobol/full_N4096`.

The job uses a frozen read-only code snapshot. Float64 outputs are losslessly
compressed; solver flags, per-oxygen eigenpair QC and exception matrices remain
auditable. Each task checkpoints every eight oxygen conditions and resumes
only after identity/axes/completed-point checks. Completed blocks are never
discarded to fix a missing block. The cluster's own CPU allocation limits still
apply; this submission imposes no array concurrency cap.

Individual sources/figures will appear under `global/full/global` and
`neighborhood10pct/full/seed*`; ensemble summaries and figures under each scope's
`summaries` and `figures`. Undefined indices remain undefined, with endpoint
counts explicit. Every fitted endpoint remains included regardless of diagnostics.

After final completion, collect deliverables and verify source checksums:

```bash
python3 oxygen/code/O2_supply_demand_MAP/analysis/eFAST/sobol/collect_pilot.py \
  --remote-dir /share/lab_crd/taoli/Project/soft_couping_eFAST/oxygen/results/eFAST/sobol/full_N4096 \
  --local-dir oxygen/results/eFAST/sobol/full_N4096
```

Raw repeated float64 arrays and exception matrices stay on HPC; summary receipts
pin the block/repeated-array checksums. See the code README for full methodology.
