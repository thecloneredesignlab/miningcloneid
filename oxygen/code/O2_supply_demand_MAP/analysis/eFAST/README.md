# Independent in-vivo eFAST companion to Figure 4B

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
than a high-precision uncertainty estimate. All four heatmaps use the same white-to-deep-purple palette, while
ploidy and growth each use their own data-driven color maximum and horizontal
colorbar. The oxygen axis uses a base-10 symmetric-log scale with a 0.025%
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
