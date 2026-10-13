# N=1025/2049 five-phase neighborhood eFAST pilot

## Submission on 2026-10-13

Scientific code commit: `47055427002b86bc8796d4fb3ee0df218d3ac530`.
Smoke job **20833826** completed with exit `0:0`; all nine tests ran and passed
in the requested SIF. Real full-N trajectories at three oxygen values passed
for seed174/N1025 and seed494/N2049, together with the actual-matrix 50/100-digit
recovery validation (double residual `2.7755575615628914e-17`).

Pilot trajectory array **20833838** (`1-1680`), seed/N summary array **20833839**
(`1-24`, afterany:20833838), and comparison job **20833840**
(`afterany:20833839`) were submitted. The recorded controller snapshot confirmed
all 1,680 trajectories running; summaries and comparison were waiting on their
dependencies. This is a submission snapshot, not a completed-analysis claim.
Actual controller settings confirmed xxlarge, 12 hours, no requested node and
no array concurrency throttle.

`slurm/submission_jobs.tsv`, `slurm/resource_verification.json`,
`smoke/smoke_receipt.json` and `collection_receipt.json` record the evidence.
All ten collected baseline/submission/preflight artifacts matched HPC SHA256
hashes in `artifact_collection_manifest.tsv`. Subsequent provenance commits do
not alter the frozen scientific code. New high-N outputs and figures will be
generated on HPC by the background jobs and require collection after completion.

This pilot evaluates the separate **in-vivo** fixed-oxygen operator at 12 fitted
endpoints. It retains all 14 active parameters and the original 201 oxygen values
(0–5%, step 0.025%). Outputs are dominant mean ploidy and asymptotic net live-cell
growth, the leading eigenvalue. N is the number of samples **per focal FAST
trajectory**, not a chromosome count.

## Design and endpoint selection

Each endpoint retains its audited natural-scale neighborhood:
`[max(L, best - 0.10*(U-L)), min(U, best + 0.10*(U-L))]`.
The original identity/log10 transforms and independent transformed-uniform
distributions are inherited. No new optimum, parameter range, model matrix,
oxygen grid, or package environment is introduced.

The original ten representatives are seed25, seed464, seed165, seed2, seed214,
seed138, seed66, seed395, seed101 and seed311. They cover fit quality, clipped
bounds and normalized parameter space. Add seed174, the largest pooled ploidy
ST five-phase range p90 outside those ten, and seed494, the largest key-parameter
ploidy ST phase-range p90 among the remaining endpoints. The selection script
recomputes those criteria from all 500 N=513 source arrays and refuses a different
selection. `seed_selection.tsv` and `selection_diagnostics_all500.tsv` document
the endpoints and their criteria.

Use SALib 1.5.2, M=4, N=1025 and 2049, five repeats each. Reuse the phase RNG
seed `202610080 + endpoint_number*10 + repeat - 1`. Each N is sampled afresh;
larger designs are not interpolations or subsets of smaller designs.

The frozen original N=513 five-phase indices for these 12 endpoints are copied
to `baseline_N513.npz` after verifying the original collection hashes. Their
source paths/hashes are in `slurm/submission_plan.json`. Original ensemble
outputs and figures remain the reference 500-endpoint analysis.

## Reproduce and submit on RED

From the current branch in `/share/lab_crd/taoli/Project/soft_couping_eFAST`:

```bash
code=oxygen/code/O2_supply_demand_MAP/analysis/eFAST
bash "$code/submit_resolution_pilot.sh" prepare
bash "$code/submit_resolution_pilot.sh" smoke
# After the smoke receipt says passed:
bash "$code/submit_resolution_pilot.sh" submit
```

The submitter checks the passed compute-node smoke against the immutable plan
hash before submitting the pilot. All jobs request `--qos=xxlarge` and
`--time=12:00:00`, with no specified node and no array concurrency limit:

| Role | Number | CPU / memory | Unit |
|---|---:|---|---|
| Smoke | 1 | 1 / 4G | Actual model at full N, three oxygen values, plus 50/100-digit validation |
| Trajectory array | 1,680 | 1 / 4G | seed × N × phase × focal parameter, full oxygen grid |
| Summary array | 24 | 1 / 8G | seed × N, 70 trajectory caches |
| Comparison | 1 | 4 / 32G | All pilot endpoints and N=513 baseline |

The pilot contains 519,014,160 new fixed-oxygen operator evaluations. Summaries
run after all trajectory jobs terminate (`afterany`); they require every relevant
completion receipt. Comparison requires all 24 valid summaries. Thus a failed
trajectory is diagnosed and cannot become a silently incomplete result. Do not
re-submit an active array. For a failed unit, invoke the worker for its exact task
ID using the same frozen plan and wrapper; verified completed outputs are reused.

Every compute job receives private node-local copies of the audited Rcpp and
Perron caches. The required image is
`/share/lab_crd/taoli/Docker/o2_supply_demand_map_r442_hpc_exact_salib152_20260923.sif`.
SIF, source, fit/configuration, bounds, samples, baseline and completion receipts
are hash checked. Only invalid eigenvectors on otherwise complete canonical
trajectories receive a targeted 50/100-digit positivity-preserving recovery on
the unchanged matrix. Residual, precision agreement and the leading-root bound
are checked. Recovery proofs retain original rows, matrices and index impacts.
No output rows are dropped; constant-output indices stay NaN with variance and
defined-phase counts preserved.

## Diagnostics and figures

The five key parameters are **buffer_beta, buffer_smax, p_mis_base, mu_hp and
O2_crit**. They receive focused tables, while the same diagnostics cover all 14
sampled parameters and both S1/ST and outputs.

- `repeat_diagnostics.tsv`: per-parameter five-phase range p90/max, maximum's O2,
  SD p90, band-mean range/SD and peak variability, separately for 0–5%, 0–1% and
  3–5%. The practical phase-range p90 threshold is 0.10 absolute index units.
- `resolution_diagnostics.tsv`: absolute differences between five-phase mean
  curves at 513→1025 and **1025→2049**, with p90/max, peak/O2 changes and band
  deltas. The practical delta p90 threshold is 0.05 absolute index units.
- `key_ploidy_ST_*_diagnostics.tsv`: focused subsets for the five key parameters.
- `oxygen_curve_summary.tsv.gz`: plotted curve sources, SD, range, phase counts.
- `band_rank_mechanism_comparison.tsv` and `band_rank_mechanism_stability.tsv`:
  parameter ranks, high-minus-low contrasts and Noemi's three-condition mechanism
  signature at each N and its stability across N. Process scores are means of
  their parameter sensitivities; no named pairwise interactions are inferred.
- `figures/ploidy_ST_resolution_comparison.pdf` and the growth companion: each
  row is a selected endpoint, each column a key parameter; curves are five-phase
  means, shaded envelopes are full phase ranges. Gray N513, purple N1025, orange
  N2049. Oxygen windows 0–1% and 3–5% are lightly shaded.

`passed_by_threshold` requires p90 and maximum below threshold across every
expected oxygen point. `localized_review` means p90 passes but the maximum does
not. `not_passed` means p90 fails; `undefined_or_incomplete` means the five-phase
comparison is not defined throughout the window. A reported `p90_pass` alone is
not a claim that all oxygen values converge. Thresholds are diagnostics, never
endpoint-exclusion gates or formal error bounds.

These 12 endpoints are deliberately selected, including two difficult cases.
Their pooled median or support fraction must not replace the original 500-seed
figure, oxygen classification or prevalence estimate. Successful pilot checks
would justify a separately authorized larger-N extension; this submission does
not launch one automatically.

Raw trajectory outputs, sampling designs, recovery matrices, logs and scratch
caches stay on HPC. Compact sources, diagnostics, figure exports and receipts
are allowed in Git; collect them only after verifying completion and hashes.
