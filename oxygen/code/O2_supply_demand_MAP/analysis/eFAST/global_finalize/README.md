# Global eFAST completion and Sobol comparison

This stage uses only the separate in-vivo **global** model. Sobol N131072 is
read-only input; all endpoint neighborhoods remain archived.

The first approved batch is N=2049,4097,8193,16385, M=4, five independent phase
seeds, and all 201 oxygen points from 0 to 5% in 0.025% increments. There is no
upper-N restriction. Later batches can use `--root NEW_ROOT --resolutions
32769 65537 ...`; the current immutable batch is never rewritten. No additional
batch is silently submitted just because a diagnostic failed.

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
