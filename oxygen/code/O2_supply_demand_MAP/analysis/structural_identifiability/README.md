# Culture identifiability preparation

This directory contains the WP0 no-optimization replay, WP1 parameter mapping
and bounded cost benchmark, the first culture profile pilot, and the approved
follow-up optimization review. Each stage has its own frozen implementation,
state, validation and accounting. The legacy in-vivo profile runner is not
used for these culture studies. The original plans and the approved
optimization-review plan are under `docs/`.

`replay_culture_baseline.R` consumes an explicitly prepared, isolated study
directory. It sources the frozen production model and culture helpers, uses
the archived full-precision parameters/configuration, and evaluates the
forward model exactly once. It does not source a fitting backend or call an
optimizer. The archived `best_components` are comparison targets, not replay
predictions. An independent implementation recomputes all three score streams
from observations and selected predictions.

Required study inputs:

```text
source/seed144/fit_result.rds
source/seed144/parameter_table_input.csv
source/code/oxygen/code/O2_supply_demand_MAP/{model,util}/
source/inputs/oxygen/ploidyOxygen/data/fit_objects/{fit_data,jobs_2N,jobs_4N}.Rds
source/inputs/oxygen/data/g0g1_ploidy_density_grid.csv
manifests/replay_settings.json
```

Freeze hashes, source commit, tolerances, weights, and execution settings
before running. The prepared 2026-11-09 study is:

```text
oxygen/results/structural_identifiability/20261109_seed144_wp0/
```

It archives the complete seed144 production bundle, source commit
`1a5d861c9c4d2af72ea0d07def0023831cc531e4`, the iteration5 comparison
manifests, and the single-job HPC launch script. Its report records the
historical container-identity caveat and the initial runtime-entrypoint
failure. Production source and input files are bound read-only during replay;
the C++ cache is mounted separately under the study's output directory.

Use a new study root for a new replay. The R script refuses an already
completed root. In the HPC-exact SIF, set
`RSCRIPT_BIN=/usr/local/bin/o2-hpc-exact-rscript` before invoking the repository
`scripts/agentRrunner.sh`; generic `Rscript` bypasses the image's required
HPC library-path initialization. The prepared launch script demonstrates the
complete invocation, resource limits, and read-only mounts.

After replay, rebuild its report and production contract locally with:

```bash
python3 oxygen/code/O2_supply_demand_MAP/analysis/structural_identifiability/summarize_culture_replay.py \
  oxygen/results/structural_identifiability/20261109_seed144_wp0
```

The summarizer requires passed replay checks and unchanged frozen inputs. It
does not run the model. Raw replay tables/RDS, independent scores, full
distribution comparisons, observation/conditioning manifests, logs, and
runtime fingerprints remain inspectable. Run fingerprints describe this
replay, not a statistical likelihood or a confidence region.

## WP1 culture adapter

The shared helpers are:

- `util/o2_supply_demand_map_profile.R`: active/fixed/inactive mapping,
  natural/free conversion, actual start populations, sample contracts,
  candidate acceptance, and fingerprint/overwrite guards.
- `util/o2_supply_demand_map_profile_culture.R`: load the immutable, verified
  production v2 source into an isolated environment and adapt its specification,
  start-population generator, and observation validation. The existing frozen
  backend still performs DEoptim, worker/DLL preflight, L-BFGS-B refinement,
  protocol penalties, simulation, and scoring. No second optimization stack or
  simulator is introduced.

The ordinary legacy culture CLI in this checkout is unchanged. New diagnostics
must use the explicit v2 adapter, not assume the legacy CLI has been upgraded.
Neither production CSVs nor the frozen WP0 source tree are edited.

`zero_extended` admits natural-scale zero for both `p_misseg` and `buffer_beta`,
with fixed transforms `log1p(p_misseg/0.015)` and `log1p(buffer_beta/1)`.
Other production transforms/bounds remain unchanged. Fixed zeros bypass the
transform entirely. The unrestricted production and extended maps both have
20 active coordinates. At `p_misseg=0`, `k_o_mis` is inactive for this data-only
culture objective; at `buffer_beta=0`, `buffer_n_exp` is inactive while
`buffer_smax` stays free. Each single null therefore has 18 active coordinates.
The generic mapper retains a data-inactive parameter if it still enters a
complete-objective penalty. The culture adapter currently supports only the
verified production data-only objective, not new priors or joint penalties.

Starts contain all 20 fitted natural coordinates. Projection onto exact
constraints/inactive values is saved explicitly. A 64-machine-epsilon boundary
tolerance handles production log-transform roundoff; any adjustment is labeled
`roundoff_to_bound`. It is not an epsilon replacement for zero.

## Validation and single benchmark

The prepared engineering study is:

```text
oxygen/results/structural_identifiability/20261109_culture_v2_wp1/
```

`test_culture_profile.R <project> <study>` runs pure mapping/orchestration tests,
the existing kernel/event-accounting oracles, integrated production and
extended-domain replays, exact-zero/PSOCK tests, and a native DEoptim entry spy
with a cheap analytic objective. Test fixtures get fresh directories. The test
runner disables automatic legacy helper loading so that only the frozen v2
model is used.

`run_culture_constrained_fit.R <project> <study>` refuses to run without passed
tests for the current helper fingerprints and the expected SIF SHA256. It
launches exactly one constrained fit, injects declared complete starts into the
native DE population, records the population actually received by DEoptim, and
logs evaluator calls by process/phase. It independently replays and recomputes
the saved endpoint score. Its output remains a candidate until scheduler
completion is independently checked.

After retrieving the outputs and pipe-separated `sacct` accounting into
`manifests/benchmark_scheduler.tsv`, finalize without loading a simulator:

```bash
scripts/agentRrunner.sh \
  oxygen/code/O2_supply_demand_MAP/analysis/structural_identifiability/finalize_culture_benchmark.R \
  /absolute/project/root /absolute/engineering/study/root JOB_ID
```

The finalizer checks external exit status, source/config/start/result hashes,
exact target values, full nuisance feasibility, retained sample identities,
and independent score agreement. It refuses completed outputs. The guards
support safe reuse of compatible saved endpoints as new starts. This benchmark
does not implement interrupted-DE continuation or batch scheduling; the later
pilot and review have separate orchestration.

The benchmark fixes `p_misseg=0`, uses the production 256-member population,
at most 100 DE iterations and 50 local iterations, and eight cores. It measures
that declared bounded search, not a converged profile minimum. At the benchmark
stage, an extended-domain unrestricted baseline and multiple starts per profile
point remained future work; they are implemented in the subsequent pilot.
No confidence thresholds, intervals, COMBOS calculations,
manuscript regeneration, or production-fit modifications are performed.

The completed 2026-11-09 benchmark is externally validated: job 21798460,
75m25s total job wall time, 3728.97 MiB Slurm MaxRSS, and 29,242 evaluator
invocations. The local optimizer returned code 0, while DE reached its
100-generation budget; the endpoint is provisional rather than a validated
profile minimum. Eight nuisance coordinates reached a bound. Full evidence is
in the study's `reports/culture_v2_engineering_report.md` and
`reports/engineering_summary.json`. `manifests/code_manifest.tsv` points to
archived implementation snapshots, and `manifests/artifact_manifest.tsv`
fingerprints the final study artifacts, excluding runtime caches.

## Culture pilot and optimization review

The parent pilot is `oxygen/results/structural_identifiability/20261109_culture_v2_pilot`.
Its 125 accepted attempts and 20 original point candidates remain immutable.
See `CULTURE_PROFILE_PILOT.md` for its implementation. Accepted feasible
candidates alone do not establish independent-start recovery or a converged
profile minimum.

The approved follow-up is
`oxygen/results/structural_identifiability/20261110_culture_v2_optimization_review`.
Its scope, priorities, tolerances and stopping rules are recorded in
`docs/culture_profile_optimization_review_plan.md`. It retains the production-v2
culture simulator, scoring and selected-day protocol. No tumor fitting is
included. Exact zero constraints, fitted initial states and observation-error
coordinates retain the verified WP1 mapping.

`util/o2_supply_demand_map_review.R` adds separately allocated global and local
execution. The global adapter intercepts the frozen backend's local entry,
without calling L-BFGS-B. The inherited native trace is retained verbatim and
is not evidence that a local optimizer ran. Actual-call tests, requested and
received populations, minimum-generation traces and fresh-process score
checks are saved. Derivative-free local searches and one-sided boundary
probes use the full objective with normal passage-day selection.

`culture_review.py` is a budgeted, idempotent dependency controller. Resolve
the deployed immutable version through `manifests/orchestrator_versions.json`;
do not execute an old controller or copy a stale local state to HPC. It
distinguishes fresh independent, transferred, parent and local origins.
Independent initialization banks record every proposal and use first-feasible
members without ranking by objective. Screening changes never narrow model
bounds or insert accepted endpoints into pure-independent populations.

Every job uses `xxlarge`. Its verified single-job wall limit is 12 hours;
globals request eight cores and locals one. No job-concurrency cap or array
percentage throttle is imposed. Actual closed allocations plus worst-case
live commitments are checked against the new 10000-core-hour ceiling before
submission. Attempts and reserve rounds have separate hard caps. Parent-pilot
costs are reported separately.

P2 and the declared D_base/D_sat sentinels depend on scientific recovery
gates. Each expanded domain needs its own unrestricted reference. A failed
gate withholds dependent work, retains already submitted results and remains
unresolved at the declared stopping rule. Warm agreement cannot substitute
for two completed fresh populations near the current minimum. Physical
probability limits and arbitrary search ceilings remain distinct; the
chromosome grid is unchanged.

`watch_culture_review.py` reads state only. `sync_culture_review_closed.py`
captures HPC state, pulls closed artifacts, verifies accepted payload hashes
and installs the captured local state last; it never pushes state. The final
reporting snapshot replays accepted candidates inside the fingerprinted
runtime and separates original curves from domain sentinels. It preserves
raw negative deltas, complete distributions, component/function traces and
recording-scope event accounting. `audit_culture_review_search.py` summarizes
actual evaluation logs by phase; worker-call seconds are not allocated
core-hours. Reports remain diagnostic profile objectives, without formal
confidence intervals or likelihood-ratio thresholds.

Evaluation-error exports group protocol failure reasons while retaining cohort labels.
`check_culture_search_audit_format.py` runs beside a frozen source manifest and
reconciles the bounded category export against the archived v03 metrics. Raw
worker errors retain candidate-specific numeric details; summaries keep those
details out of oversized spreadsheet cells. This metadata-only repair does not
change fitting or prediction output.
