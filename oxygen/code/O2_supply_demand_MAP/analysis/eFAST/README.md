# Independent in-vivo eFAST companion to Figure 4B

For the separate full-bound and per-endpoint **second-order Sobol pilot** and
the parameter/equation audit, see [sobol/README.md](sobol/README.md).

`efast.py` generates SALib 1.5.2 FAST trajectories and analyzes the two outputs
of the **fixed-oxygen** operator. `evaluate_fixed_o2.R` sources the canonical
`simulation/o2/fixed_o2/run_fixed_o2_simulation.R` and calls its
`fixo2_dominant_attractor_one` for every parameter vector and oxygen value.
It verifies the selected fitted endpoint (seed25 for the global design) against
three Figure 4 reference points before evaluating a design. The leading real eigenvalue is the
asymptotic net live-cell growth rate in day^-1; the normalized leading
eigenvector determines dominant mean ploidy.

## Input and sampling contract

- Independent **in-vivo** fitting run: `parameter_table.csv` supplies exact
  transformed optimizer bounds; `seed25/best_params.tsv` and
  `seed25/fit_config.rds` supply fixed model settings and reference validation.
- Figure 4 data: `fixed_o2_dominant_ploidy_201grid.tsv` verifies the 201-point
  oxygen grid from 0 to 5% in 0.025% increments and the seed25 reference;
  `continuous_ploidy_spearman_by_o2.tsv` supplies the existing directional
  correlation panel; `parameter_function_groups.tsv` supplies Figure 4 row
  order. Both are copied to the result folder with recorded SHA-256 values.
- Fourteen parameters present in the fixed-oxygen operator are sampled
  independently. Transformed `log10_*` bounds are sampled uniformly in log10
  space (log-uniform on natural scale); identity bounds are sampled uniformly.
  These priors cover documented fit ranges. The 500 optimizer solutions do not
  constitute a FAST design and are not resampled as one.
- `o2_S0`, `kappa_O`, `eta_o2`, and `k_clear` are Figure 4B parameters whose
  direct effects are structurally absent when oxygen is externally fixed;
  their sensitivity indices are undefined for this 14-factor design. The
  eFAST heatmaps omit these parameters instead of displaying empty rows.
- `O2_crit` and `n_O` control the shared oxygen-stress response used by growth,
  death, and missegregation functions. They are reported separately from the
  death-specific `mu_hp` and `gamma_mu` when comparing mechanisms.
- Whole-genome doubling `p_wgd` is summarized separately from missegregation;
  it is a genome-change process rather than a missegregation probability.
- The 201 oxygen points use the same FAST parameter vectors within each
  resolution/replicate. FAST trajectory order is preserved, including if
  numerical evaluation fails: failures stop analysis rather than silently
  removing rows.

## Execution

### Six-endpoint N4097 resolution pilot

Completed-result collection and independent interpretation use
`collect_resolution_4097_results.py` and `report_resolution_4097.py`.
The collector runs with Python 3.6+ standard libraries; the local report
requires NumPy>=2.0 and Matplotlib. They read finished results and never
execute the model or submit jobs. Complete commands and the detailed
interpretation are in
[`N4097/interpretation_zh.md`](../../../../results/eFAST/neighborhood10pct/pilot/N4097/interpretation_zh.md).
The collector's `audit --repo-root ...` command runs on HPC; `manifest --root ...`
writes an explicit compact-artifact allowlist, and `verify-local --root ...`
checks the collected files against their original completion hashes. Original
submission receipts are retained separately from the completed-results manifest.

`resolution_4097.py` adds the approved N=4097 pilot for **seed25, seed101,
seed138, seed174, seed66, seed311**, with all 14 inputs, M=4, five phases,
and the full 201-point O2 grid. It reuses the unchanged `resolution_pilot.py`
evaluation, trajectory splitting, certified numerical recovery and seed summary.
The original N1025/N2049 code and receipts remain valid.

After committing/pushing and synchronizing the branch on RED:

```bash
bash oxygen/code/O2_supply_demand_MAP/analysis/eFAST/submit_resolution_4097.sh prepare
bash oxygen/code/O2_supply_demand_MAP/analysis/eFAST/submit_resolution_4097.sh smoke
# Once smoke/smoke_receipt.json reports passed:
bash oxygen/code/O2_supply_demand_MAP/analysis/eFAST/submit_resolution_4097.sh submit
```

If smoke is waiting for xxlarge capacity, the alternative
`queue_resolution_4097_after_smoke.sh` queues the same arrays with an
`afterok` dependency on smoke. Evaluation still requires a successful smoke
job. Its source hash and submission commit are recorded separately so the
already prepared scientific source hashes remain intact.

The submitter requests `--qos=xxlarge --time=12:00:00`, uses `--array=1-420`
without a concurrency suffix, and requests no named node. Each 1-CPU/4G task
evaluates one complete 4097-vector focal FAST trajectory at all 201 O2 points.
Six 1-CPU/8G seed summaries and one 4-CPU/32G final comparison use `afterany`
dependencies and require complete, validated inputs before calculating results.
Missing trajectories stop aggregation; convergence failures are retained as diagnostics.

Results are separate under
`oxygen/results/eFAST/neighborhood10pct/pilot/N4097/`. Frozen comparison
arrays cover N513, N1025 and N2049 for these six endpoints. Finalization
exports all-output S1/ST tables, full phase ranges and means, key-parameter
ST curves, and a ploidy ST stability overview. The primary resolution check
is **2049 to 4097**. Across overall [0,5], low [0,1], and high [3,5] O2 windows,
the phase-range p90 threshold is 0.10 and mean-curve change p90 threshold is
0.05, in absolute index units. Local maxima are checked separately; p90
passes with a failing maximum are labeled `localized_review`. These are
practical diagnostics, not confidence intervals or a population estimate.

### Original global analysis

On `hpctpa3pc0028`, with the SALib 1.5.2 SIF and direct shell access:

```bash
bash oxygen/code/O2_supply_demand_MAP/analysis/eFAST/run_efast.sh pilot
bash oxygen/code/O2_supply_demand_MAP/analysis/eFAST/run_efast.sh full
EFAST_WORKERS=16 EFAST_N=513 bash oxygen/code/O2_supply_demand_MAP/analysis/eFAST/run_efast.sh extend
```

The combined publication figure can be redrawn without recomputing eFAST:

```bash
python3 oxygen/code/O2_supply_demand_MAP/analysis/eFAST/efast.py plot \
  --out-dir oxygen/results/eFAST \
  --figure4-dir oxygen/results/eFAST \
  --figure4-layout-dir /path/to/HypoxiaLTEEFigures/revised/iteration5/data/Figures/Figure4 \
  --combined-only
```

`--figure4-layout-dir` copies the iteration5 Figure 4B oxygen classification,
row ranking, parameter-process annotations and palette into the result folder.
It can be omitted on later redraws because the copied source tables are tracked
with the eFAST results. The combined figure uses the Figure 4B oxygen-group
classification after filtering to the 14 sampled parameters as one row
annotation. A second oxygen annotation applies the same Figure 4B decision
rule to eFAST: normalized sensitivity AUC in Low [0,1] and High [3,5], a global
peak strictly greater than 0.3, and a phase-repeat bootstrap contrast at BH
q < 0.05. The eFAST classification uses dominant mean ploidy only. Panel A is
ordered by the ploidy S1-derived group and peak; panel B is ordered by the
ploidy ST-derived group and peak. The growth heatmaps follow these ploidy-based
orders and do not affect classification. The bootstrap unit is one complete
phase-repeat curve. The original figure used two repeats; the extended runner
uses five. This classification remains a repeat-consistency diagnostic rather
than a high-precision uncertainty estimate. S1 uses white to deep purple;
ST uses white to orange. Ploidy and growth each use their own data-driven
color maximum and horizontal colorbar. The oxygen axis uses a base-10 symmetric-log scale with a 0.025%
linear threshold so the observed 0% point remains visible. Black dashed
vertical lines at 1% and 3% oxygen mark the low-window upper boundary and the
high-window lower boundary, respectively.

The global pilot uses N=129, M=4, five phase seeds and four Figure 4 oxygen values.
The full analysis uses N=129,257,513, M=4, five independent phase seeds at each
resolution and all 201 oxygen values. Each design contains `14*N` model vectors.
The `extend` mode brings `EFAST_N` (513 by default) to five phases. The original
R1/R2 designs are retained after verifying metadata and sample hashes; R3–R5
are added. Existing metadata mismatches stop before overwriting a design.
`EFAST_WORKERS` controls local fork workers (default 4). The runner enforces
the named node and does not submit a Slurm job. The full run is restartable:
an existing nonempty `outputs.tsv.gz` is retained. For a fresh rerun of one
design, remove only that design's output first.
The runner clears inherited R settings so that `Matrix` and other packages
come from the validated SIF rather than an incompatible HPC home library.

## Per-fit-seed neighborhoods

`neighborhood.py` and `audit_neighborhood_inputs.R` add the separate analysis
of all 500 in-vivo endpoints. Natural bounds are `best +/- 0.10*(U-L)`, clipped
to the original bounds; the original log10/identity distribution is retained.
Every endpoint has its own FAST trajectories and five phase repeats. The
runner performs a runtime pilot and a resolution pilot before the
N=513 full run. Convergence checks are diagnostic only and do not block any of
the 500 audited endpoints. See [the neighborhood protocol](../../../../results/eFAST/neighborhood10pct/README.md)
for selection, convergence tolerances, output axes and interpretation limits.

```bash
# Execute directly on hpctpa3pc0028 using the existing SALib 1.5.2 SIF.
EFAST_WORKERS=16 bash oxygen/code/O2_supply_demand_MAP/analysis/eFAST/run_neighborhood.sh pipeline
# Or inspect/run individual stages: audit, smoke, convergence, full.
python3 oxygen/code/O2_supply_demand_MAP/analysis/eFAST/test_neighborhood.py
```

The pipeline audits all endpoint parameter tables/configs and checks 1,500
original Figure 4 reference outputs. The pilot saves convergence diagnostics
and continues to the 500-endpoint full run regardless of convergence status.
Full execution retains input integrity and free-disk checks. Per-endpoint
diagnostic tables distinguish five-phase stability from resolution stability;
resolution is untested outside the multi-N pilot. Completed designs are reused only
when their completion receipts, metadata and output hashes match. Each design
is evaluated in bounded chunks; an interrupted design is recomputed in full.
The kernel lock prevents two neighborhood controllers from sharing an output root.

### Slurm trajectory arrays

```bash
# Compile a private template using the same SIF on a compute node.
bash oxygen/code/O2_supply_demand_MAP/analysis/eFAST/build_neighborhood_rcpp_template.sh
bash oxygen/code/O2_supply_demand_MAP/analysis/eFAST/submit_neighborhood_slurm.sh
# Retry only missing/failed tasks after the preceding Slurm attempt has terminated.
bash oxygen/code/O2_supply_demand_MAP/analysis/eFAST/submit_neighborhood_slurm.sh --retry
python3 oxygen/code/O2_supply_demand_MAP/analysis/eFAST/test_neighborhood_slurm.py
```

The Slurm workflow starts independently of the direct convergence pilot.
It submits 35,000 trajectory tasks (500 fit seeds x five phases x 14 focal
parameters), a dependent 500-task seed-summary array, and a final aggregation
job. Every trajectory contains all 14 simultaneously varied parameter
columns, N=513 ordered vectors and all 201 oxygen points; its focal parameter
identifies the FAST frequency block. Sampling and indices are unchanged.
Trajectory tasks request one CPU, 4G and 12 hours; seed summaries request one
CPU and 8G; final aggregation requests four CPUs and 32G. All use xxlarge and
12 hours, with no specified compute node and no `%N` array concurrency limit.
The scheduler's account/QoS limits still apply.

Submission validates the immutable SIF once, hashes inputs/code and writes
`slurm/submission_plan.json`, `task_manifest.tsv`, `submission_jobs.tsv` and
`execution_backend.json`. Workers recheck the recorded code/input hashes and
SIF size/mtime. Per-seed preparation and per-trajectory locks prevent races.
The SIF defaults to forced C++ recompilation, so the array worker explicitly
sets `MININGCLONEID_RCPP_REBUILD=FALSE` and binds a unique node-local Rcpp cache
over the model's cache path. Each cache is populated from a source/image-hashed
template built with the same SIF. Cached wrappers validate the required backend;
any compilation fallback is confined to that task. This avoids a shared
sourceCpp lock across thousands of jobs. Local caches are removed on task exit.
Completed pilot phases are copied and reused after receipt validation.
Raw trajectory outputs are assembled in original FAST row order, with hashes
and complete row counts, before averaging five phases within each seed.
No shared summary tables are written by concurrent array workers.

`afterany` dependencies allow downstream jobs to diagnose missing tasks.
Incomplete seeds/final collections fail explicitly; partial data never become
a completed full analysis. Task receipts and `slurm/status.json` identify
failed/missing task IDs for targeted resubmission. The direct pilot can finish;
its next serial full stage sees the backend marker and hands over to Slurm.
`--retry` archives the previous submission plan, appends new job IDs and keeps
completed full-phase reuse receipts. Existing valid raw trajectory outputs
are reanalyzed rather than reevaluated if their analysis cache version changes.

Checks cover natural/log bounds, clipping, five reproducible phase seeds,
agreement with SALib point estimators and a known additive model, rejection of
incomplete trajectories, and retention of undefined indices in phase summaries.

## Outputs

`oxygen/results/eFAST/` contains `parameter_ranges.tsv`, the copied Figure 4B
correlation source, full FAST samples and operator output tables under `runs/`,
`indices.tsv` with S1/ST per replicate, `convergence.tsv` with replicate ranges
and change from the highest resolution, parameter and mechanism summaries for
0–1% and 3–5% oxygen (including separate replicate band summaries), four index heatmaps, a directional
correlation panel, and `interpretation.md`. The sample metadata records seeds,
distributions, input hashes, oxygen points, and SALib version.
`spectral_gap_qc.tsv` records how often the two leading eigenvalues nearly
coincide, which is relevant to stability of asymptotic ploidy.
The top-level `environment.tsv` records the evaluator image and Git commit;
`analysis_manifest.tsv` records the postprocessing commit and summary hashes.
`figure_redraw_manifest.tsv` records the combined figure's filtered row order,
exact color limits, source-table hashes and output hashes.
`efast_o2_sensitivity_classification.tsv` records the eFAST window scores,
bootstrap tests, BH-adjusted q values, classifications and panel-specific row
orders.
S1 and ST are variance contributions without a sign; the correlation panel retains
direction. S1/ST are conditional on these chosen input ranges/distributions,
not a posterior uncertainty decomposition. Summing ST across parameters
double counts shared interactions.
## Numerical recovery of failed neighborhood points

`recover_eigen.py`, `repair_eigen_points.R`, and `perron_high_precision.cpp`
recover invalid fixed-O2 evaluations without changing the FAST design or the
model matrix. The C++ solver uses the SIF's existing Rcpp and BH headers;
no additional Python environment or package installation is required.

For a failed trajectory, preserve the original partial gzip output and every
valid row. Rebuild only its invalid points' original canonical matrices, then
use positivity-preserving Noda inverse iteration at 50 and 100 decimal digits.
Accept a point only when both precisions agree, the vector is nonnegative and
normalized, the original-matrix residual is small, and the leading-root upper
bound meets the recorded tolerance. No matrix rates or graph edges are
regularized. The spectral-gap diagnostic retains the canonical double-precision
spectrum; ploidy and the leading eigenvalue use the certified recovered pair.

The method follows the matrix form of Noda iteration described in
[Jia, Lin and Liu](https://arxiv.org/abs/1309.3926).
Each `recovery/numerical_recovery.json` records source hashes, matrix artifacts,
unchanged-row counts, original/corrected values and the effect on SALib S1/ST.
Unchecked pre-recovery indices are used only for this numerical impact check.
Successful scientific results always pass the standard validity checks.

After validation in the required SIF and a complete-trajectory pilot, run:

```bash
bash oxygen/code/O2_supply_demand_MAP/analysis/eFAST/submit_eigen_recovery.sh
```

This submits only failed/missing absolute task IDs, then only missing seed
summaries and the final aggregate. All use `qos=xxlarge`, 12 hours, no requested
compute node, and no array concurrency cap. All 500 fit endpoints remain
included; convergence remains diagnostic only. A compiled Perron template is
copied into a private node-local cache for each task. Large raw outputs and
matrix proof artifacts remain on HPC; compact source tables, manifests, figures
and reproducible scripts are collected into this branch.

The versioned `*_safe` recovery files handle points where row-pivoted linear
solves lose positivity at finite precision. For the shifted nonsingular
M-matrix, they use Gaussian elimination without row exchanges and a larger
shift safeguard (`1e-25` times the matrix scale). The original matrix stays
unchanged. The same 50/100-digit agreement, nonnegativity, leading-root bound
and residual checks apply. Original recovery source files remain immutable so
their successful receipts can still be verified. The M-matrix factorization
property is described by [Plemmons (1981)](https://doi.org/10.1016/0024-3795(81)90091-4).
`audit_numerical_recovery.py` collects point corrections and S1/ST changes from
both versions into compact, hash-checked source tables.

After final aggregation and the numerical audit have completed, redraw the
neighborhood figure and collect its compact source data:

```bash
python3 oxygen/code/O2_supply_demand_MAP/analysis/eFAST/redraw_neighborhood.py \
  --out-dir oxygen/results/eFAST/neighborhood10pct
python3 oxygen/code/O2_supply_demand_MAP/analysis/eFAST/collect_neighborhood_results.py \
  --out-dir oxygen/results/eFAST/neighborhood10pct --archive /tmp/efast_deliverables.tar
```

Run scientific postprocessing inside the specified SALib 1.5.2 SIF. The redraw
requires a completed collection, preserves the sensitivity source-table hashes,
updates figure hashes, and records `figure_rendering_receipt.json`. The original
calculation's checksummed submission plan remains unchanged. Its source commit
and the later rendering commit are recorded separately, so a palette change
does not imply that the model evaluations were recomputed. The collector
requires all 35,000 tasks, all 500 seed summaries and a passed numerical audit;
it checks every collection-manifest hash before archiving. Raw trajectories,
matrix proof artifacts, compiled caches and execution logs remain on HPC.

## N8193/M4 and N4097/M8 follow-up

`resolution_nm.py` and `submit_resolution_nm.sh` run five phases for seed25/138
at N8193, M4 and seed101/174 at N4097, M8. The first pair tests N resolution;
the second tests M-setting robustness at fixed N. All stages include M in
identity and FFT calculations, including targeted-recovery diagnostics.
See [the follow-up protocol](../../../../results/eFAST/neighborhood10pct/pilot/N8193_M4_and_N4097_M8/README.md)
for settings, outputs and scheduling. The previous scientific modules remain
unchanged and are validated by their existing SHA-256 manifests.
