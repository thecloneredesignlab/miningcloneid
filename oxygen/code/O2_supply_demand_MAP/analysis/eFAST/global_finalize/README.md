# Global eFAST completion and Sobol comparison

This stage uses only the separate in-vivo **global** model. Sobol N131072 is
read-only input; all endpoint neighborhoods remain archived.

The first approved batch is N=2049,4097,8193,16385, M=4, five independent phase
seeds, and all 201 oxygen points from 0 to 5% in 0.025% increments. There is no
upper-N restriction. Later batches can use `--root NEW_ROOT --resolutions
32769 65537 ...`; the current immutable batch is never rewritten. No additional
batch is silently submitted just because a diagnostic failed.

## Approved sequential refinement (2026-11-12)

The follow-up batch evaluates **32769, 65537, 131073**, retaining M=4, the same
five phase seeds, 14 inputs and all 201 oxygen points. Its separate output root
is `oxygen/results/eFAST/global_refinement_20261112`. Invoke `submit.py` with
`--previous-root .../global_finalize_20261112 --resolutions 32769 65537 131073`.
Use the committed source revision and an explicit `--root` for the new batch.

All three new levels are registered at launch, with model N65537 depending on
the completed N32769 report and model N131073 on the N65537 report. Reporting
success does not require diagnostic success. All three levels run regardless
of the diagnostic outcome. The completed Sobol report and previous FAST index
exports are copied and checksum-verified; **no Sobol model, estimator or report
job is submitted**. N32769 compares against N16385 from the certified batch.
All seven certified FAST resolutions and the separate historical N513 appear
in the final comparison. The sampler and estimator definitions remain fixed.

New N designs are not nested prefixes. Complete old designs are not rerun.
For a new trajectory's first point, reuse is allowed only after byte equality
of all 14 sampled values, matching oxygen axes, original model/SIF identity,
and the old model block's identity and content certificate are verified.
Other points are evaluated under their own new design. Each reused point gets
an `exact_reference_reuse.json` record; model receipts distinguish design
positions, newly evaluated positions and reused positions. Checkpoints retain
the same identity and do not re-evaluate certified oxygen chunks or blocks.

There are 15,890 model shards and 210 estimator shards, covering 3,227,362,530
parameter-oxygen design positions across the three new levels (both outputs
computed together). `comparison/design_budget.tsv` separates this design
budget from Sobol's different definition of N. Exact reused positions reduce
actual model calls and are counted separately. All stages and their scoped
monitors retain xxlarge, 12h initially, no node constraint, no array throttle,
and six-hour-threshold/six-hour-increment extension with monitor renewal.

`test_refinement_exports.py` uses explicitly synthetic temporary indices to
exercise the cross-batch comparison, independent validation and dynamic
synthesis. It never creates scientific results or invokes a model/estimator.
`finish_report.py` derives conclusions from actual completed source tables;
it no longer carries the first batch's fixed N, diagnostic counts or FAST
oxygen-preference conclusions into later runs.

## Execution

1. Commit and push the explicit new source files. Sync those Git objects to RED,
   retaining the archived checkout and source results. Verify file digests.
2. On RED, run `python3 .../global_finalize/submit.py --revision COMMIT`.
3. Submission freezes controllers and copies certified small source exports.
   It submits preparation, a reusable full-201-O2 model smoke block, independent
   Sobol reporting, four model arrays, four estimator arrays, four reports, and
   the final comparison. All are xxlarge, initially 12 hours, one CPU, no node
   constraint, and no array throttle.
4. Each registered stage gets a six-hour-threshold / six-hour-increment monitor.
   Monitors also use xxlarge / 12h, renew before their own end, and call verified
   copies of the user's `/home/4482173/submit_hourly_add_hours.sh` and
   `add_hours_to_running_array_jobs.sh`. Targets are restricted by an immutable
   registry and read back after extension.

`--export=ALL` is passed without a variable list, avoiding the get-user-env
retrieval path that held a previous Sobol task. Root/role/N are positional
arguments, not export-list entries.

## Numerical identity and recovery

The R evaluator differs from the frozen Sobol evaluator only in its method
label and bootstrap source location. All fixed-O2 matrix construction,
configuration checks, dominant eigenpair calculation and 50/100 digit fallback
remain identical. Preparation verifies the SIF, frozen model sources, bounds,
and global seed25 inputs. It compares new-evaluator outputs against cached
Sobol outputs, including any recovered samples in that cached block.

The sample design is SALib 1.5.2 `fast_sampler.sample` in encoded space. Log10
inputs are exponentiated afterward. Direct `fast.compute_orders` avoids unused
bootstrap intervals; equivalence to `fast.analyze` S1/ST and an analytic
additive reference are checked before production. No custom sensitivity
estimator is substituted.

Each model task owns a nonoverlapping <=1024-sample part of one parameter
trajectory and one phase. It checkpoints every eight oxygen points and keeps
sample order. Complete output files require matching identity and SHA256
receipts. Retries skip certified blocks and resume certified oxygen chunks.
Different eFAST N/M designs are not presumed nested. Existing Sobol outputs
are not treated as eFAST samples.

## Diagnostic-only convergence

Per output / parameter / index / oxygen, report five-repeat max-min and
adjacent-new-N mean absolute difference. References are 0.05 and 0.02.
All entries remain in the sources, figures and interpretation. First-new-N
entries have `no_adjacent_new_design`, not a manufactured pass from a zero
delta. Legacy N513 is a separate historical reference: its old metadata did
not fully certify best_params/fit_config identities. Five-repeat ranges are
not calibrated 95% confidence intervals.

## Reports and interpretation

The Sobol report includes all 91 pairs, the 36 focal pairs, ST-S1, per-band
rankings, N65536/N131072 rank stability, heatmaps and strongest-pair curves.
Process means are descriptive means of parameter indices, never group variance
fractions or additive ST contributions. Noemi's low-death / high-missegregation
/ high-buffering statements are evaluated separately for both outputs.

Comparison sources retain every eFAST N, repeat spread, resolution change,
pointwise S1/ST difference from Sobol, and per-oxygen ranks. Comparison heatmaps
use Sobol ploidy row order and shared color limits within output/index.
P/C/S annotations, logarithmic positive-O2 display with zero preserved, and
1%/3% boundaries follow the approved figures. S is ploidy-based on both outputs.

Existing global response surfaces remain conditional on their fixed background;
S2 does not identify interaction direction. WGD's fitted-bound observation is
not evidence of biological absence or identifiability. Neighborhood computation,
redraw and refinement are not part of this stage.

## Validation and collection

Run `test_workflow.py` locally and in the original SIF. Preparation adds SALib
and frozen-model numerical checks. Check model/statistical receipts, source
hashes, array dimensions, ordered samples, all 5 repeats, all 201 oxygen points,
diagnostic denominators, and final tables. Collect only exported summaries,
figures, reports and receipts; raw model blocks stay on RED. Render final PDFs
locally and inspect figures before claiming final visual acceptance.

## Reporting import recovery

The first original-SIF preparation passed, including exact agreement at 603
cached Sobol model points. The first standalone Sobol report exposed a Python
module-name collision: the canonical helper directory also contains report.py.
The controller now loads its frozen reporting file by explicit path; a
regression test creates a conflicting canonical report module.

`recover_reporting.py --revision FIX_COMMIT` preserves all four original model
arrays, estimator arrays, the immutable plan and outputs. It submits only six
replacement reporting jobs using a SHA256-verified overlay, cancels only the
five superseded pending reports/final stage, and adds six scoped extension
monitors. Original failure evidence and both scheduling ledgers are retained.
The overlay cannot enter model/analyze modes. Its model/estimator/report bodies
must match the originally frozen scientific sources; only the entry import is
repaired. Nine local contract checks are repeated in the original SIF at each
replacement reporting entry.

Use the newer read-only collector with `--remote-script` when collecting the
recovery metadata. `verify_exports.py --root LOCAL_ROOT` independently validates
all repeated-index summaries, diagnostics, row coverage, comparison values and
matched color limits, and emits output-variance and method-alignment tables.

## Local closeout (Python 3.12; no production environment changes)

The local verification/plot environment is recorded in
`oxygen/results/eFAST/global_finalize_20261112/validation/local_environment.json`.
`requirements-local-python312.txt` describes only that local environment; the
HPC computation continues to use the original SIF and its original library lock.

After the final HPC report completes, run these scripts in order:

1. `collect.py --root LOCAL_ROOT --remote-script HPC_REPO/oxygen/code/O2_supply_demand_MAP/analysis/eFAST/global_finalize/collect.py`
2. `verify_exports.py --root LOCAL_ROOT`
3. `reuse_global_responses.py --source oxygen/results/eFAST/sobol/full_N4096 --root LOCAL_ROOT`
4. `publication.py --root LOCAL_ROOT`
5. `finish_report.py --root LOCAL_ROOT`
6. `validate_pdfs.py --root LOCAL_ROOT --pdftoppm /path/to/pdftoppm`

The final synthesis is in `interpretation_zh.md`; publication figures and English
Methods/Results are in `publication/`. The original checksum-certified exports
are not replaced. Publication redraws correct title/oxygen-label overlap in the
strongest-pair and comparison curves, without new models or estimators.

Independent checks also test the known zero effect of O2_crit/n_O at O2=0,
verify the original global-reference input identities, and validate band means
with a separate trapezoidal calculation. Final diagnostic band counts are
disjoint; original HPC inclusive-window counts are retained only as provenance.
The final narrative is scoped to this approved M4 batch and is reviewed against
the source tables; it is not a scientific conclusion template for future runs.
