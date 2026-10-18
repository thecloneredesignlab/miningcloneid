# In-vivo parameter map and separate second-order Sobol pilots

The fixed-O2 matrix and fitted semantics are unchanged. This workflow generates
new SALib 1.5.2 Sobol/Saltelli designs with second-order calculation enabled;
FAST outputs are not relabeled or supplied to the Sobol analyzer.

## Contracts

- All 14 active fixed-O2 parameters and all 91 unordered pairs.
- Global: original independent in-vivo uniform/log-uniform bounds; one design
  per scramble, using seed25 for fixed configuration/reference checks.
- Neighborhood: optimum ±0.10 times the full **natural** parameter span,
  clipped to the original bounds; preserve original transform/distribution.
- Eight representative endpoints: seed25,464,165,2,214,138,66,395, from the
  existing audited best/worst/most-clipped/farthest-point selection.
- Oxygen 0,0.5,2.5,5%, all on Figure 4's 201-point reference grid.
- N=256,512,1024; five independent randomized Sobol scrambles. A scramble's
  unit-cube design is paired across scopes for comparisons, while the five
  repeats are independent. Nested prefixes are verified against direct sampling.
- Exactly 30N parameter vectors per second-order design. Blocks A,AB_i,BA_i,B
  are evaluated independently and restored to SALib's interleaved row order.
- No convergence-based endpoint exclusion. Preserve negative estimates and
  near-constant-output undefined values; never delete failed model sample rows.
- Eigenpair residual/nonnegativity checks; if needed certify the unchanged
  matrix at 50/100 decimal digits and retain the matrix and correction record.
- Both outputs from each matrix evaluation: dominant mean ploidy and leading
  eigenvalue (net live-cell growth rate, day^-1).

The method map is a mechanistic audit. It is not a formal COMBOS/structural
identifiability result. Numerical S2 is a pure pairwise variance component;
ST-S1 can additionally include interactions of order three or higher.

## Run on RED

After committing/pushing and pulling this branch on RED:

```bash
bash oxygen/code/O2_supply_demand_MAP/analysis/eFAST/sobol/submit_pilot.sh
```

Defaults use the existing independent in-vivo fitting directory, iteration4
Figure 4 data, completed eFAST input audits/private Rcpp template and
`o2_supply_demand_map_r442_hpc_exact_salib152_20260923.sif`. Environment overrides
are SOBOL_REPO_ROOT, SOBOL_FIT_ROOT, SOBOL_EFAST_ROOT, SOBOL_FIGURE4_DIR,
SOBOL_SIF and SOBOL_OUT_DIR. Output defaults to `oxygen/results/eFAST/sobol`.

Slurm: xxlarge,12h, no specified node or array concurrency cap. Jobs comprise
one runtime validation job (seven numerical tests and nine reference/cache
checks), 1350 evaluation blocks, nine case summaries, an index/selection report,
54 conditional response surfaces and a final report. Each worker has a private
node-local Rcpp cache initialized from the audited existing template. All
submitted model/helper/code hashes are checked; editing them during a run stops
that run rather than silently mixing source versions.

Evaluation uses 5,529,600 operator calls; both outputs are obtained together.
N256/N512 reuse checked prefixes. Case summaries average five scrambles;
neighborhood scope figures then take the median of eight endpoint means, and
remain explicitly pilot results. There is no automatic 500-endpoint extension.

## Results and diagnostics

`model_structure/` has 18 parameter rows (14 active, four absent), 91 candidate
pair links, equations, source locations and hashes. `global/pilot/global/`
and `neighborhood10pct/pilot/seed*/` contain bounds, reference optima, ordered
samples, A/B/mixed-block outputs, per-repeat S1/ST/S2, per-case summaries and
numerical QC. Scope reports and 91-pair heatmaps live under each scope's
`pilot/`. Two outputs have separate color scales. Raw negative S2 estimates
are displayed with a diverging scale, and unstable cells have an × mark.

SALib bootstrap half-widths use 200 resamples,95% nominal confidence. They are
reported as estimator diagnostics, not calibrated randomized-QMC confidence.
Independent repeat and resolution checks use operational tolerances: repeat
range <=0.05 and mean difference N512 to N1024 <=0.02. These thresholds do not
establish mathematical convergence. All checks apply to individual parameter
or parameter-pair/output/oxygen cells; a pooled percentile cannot mask a
specific unstable interaction.

Six response pairs are chosen transparently after all 91 pairs are analyzed,
covering both scopes and outputs before filling by estimated strength; stable
candidates rank first. If estimates are unstable, they remain labeled
candidates. Conditional surfaces use nominal 41x41 grids with each optimum
inserted exactly, hold the other 12 parameters at that optimum, and report
Y(i,j)-Y(i,j0)-Y(i0,j)+Y(i0,j0). They are deterministic model diagnostics, so
the five randomized repeats apply to the Sobol designs rather than five
identical evaluations of each surface. Ploidy/growth and their contrasts have
separate scales; each output's scale is aligned across oxygen conditions.

`status.json` is complete only when all nine cases and all 54 response surfaces
are complete. Scheduler terminal state/exit code and numerical/analysis status
are separate records. Dependencies diagnose incomplete arrays; there is no
silent partial aggregation. `full_scale_projection.json` extrapolates measured
pilot CPU time/storage for 501 domains (500 neighborhoods plus one global) and
201 oxygen values at N1024, explicitly excluding queue/initialization/analysis
overhead. Convergence diagnostics determine whether that resolution is adequate.

Local meaningful numerical checks:

```bash
python3 oxygen/code/O2_supply_demand_MAP/analysis/eFAST/sobol/test_sobol.py
```

Raw matrix blocks and large ordered arrays remain on HPC. Tracked deliverables
include input/provenance manifests, source tables, diagnostic summaries,
response-surface source tables, figures and interpretation. Reproduction
requires the recorded external fitting inputs and immutable SIF.

Collect only after the final report job is COMPLETED with exit 0:0:

```bash
python3 oxygen/code/O2_supply_demand_MAP/analysis/eFAST/sobol/collect_pilot.py
```

The collector checks completion, verifies copied checksums against the HPC
inventory, and records exclusions and the final Slurm accounting state.
Surface anchors tolerate 64 machine epsilons of endpoint conversion roundoff;
the actual optimum is retained and the Sobol sampling bounds are unchanged.
For the initial surface-boundary correction, `resume_postprocessing.sh`
archives the original plan, verifies all nine original index receipts, records
the changed plotting/test sources, and submits only the remaining downstream
stages. Existing evaluated response surfaces require a separate audit.

## Resolution refinement: N2048 and N4096

The separate follow-up run keeps the same eight endpoints, four oxygen values,
14 parameters, 91 pairs and five scrambled repeats. It analyzes N=1024,2048,4096,
using N2048→N4096 for the resolution diagnostic and retaining every endpoint.
The earlier pilot remains under its original output root.

## Full N4096 run: all 500 fitted neighborhoods and 201 oxygen points

```bash
bash oxygen/code/O2_supply_demand_MAP/analysis/eFAST/sobol/submit_full.sh
```

Output: `oxygen/results/eFAST/sobol/full_N4096`. One global domain plus
seed1 through seed500, all 14 active factors/all 91 pairs, five independently
scrambled repeats per domain, oxygen 0:0.025:5 from the hashed Figure 4 grid.
N4096 has 122880 parameter vectors per design. All 75150 matrix-block tasks
evaluate their 4096 vectors at the same 201 oxygen points; both outputs come
from one canonical operator call, totaling 61870694400 calls.

Global bounds/distributions match the independent in-vivo fit. Neighborhoods
remain optimum ±10% of the full natural span, clipped; log-uniform parameters
retain their transform. Unit-cube designs are paired across domains, while the
five scrambles are independent. Shared-unit transport is checked against the
canonical SALib design exactly. N2048 prefixes provide diagnostics without
extra model evaluation. Every seed is included regardless of diagnostics.

The submission freezes the committed workflow in a read-only code snapshot.
The canonical evaluate.R and mathematical model are unchanged. Each task
uses private node-local Rcpp caches, initially warmed once for the snapshot.
The workflow is bootstrap/tests → 501 configuration/storage checks → validation
report → all 75150 evaluation tasks → 501 summaries → ensemble report →
54 representative response surfaces → final report. All stages use xxlarge,
12h; evaluation 1 CPU/4G, case and ensemble summary 1 CPU/16G, final 1 CPU/8G.
There is no specified node or `%` array concurrency cap. The actual cluster
MaxArraySize must accommodate 75150; verify it before submission.

Outputs are losslessly compressed float64 arrays (NPZ), with per-point solver
flags, per-oxygen minimum gap/maximum residual, and full exception matrices.
The canonical temporary CSV is parsed on node-local scratch; it is not kept
as a second multi-terabyte copy. An atomic NPZ checkpoint is saved every eight
oxygen conditions. On resubmission, its identity/axes/CRC/completed points are
checked and already completed oxygen chunks are retained. Analysis similarly
checkpoints each eight-oxygen batch. Completion receipts checksum every block.
No bad sample row is deleted. The storage gate requires 1.35 TiB free, allowing
for roughly 980 GiB uncompressed values+flags plus reserve; compression is
measured during the 501 configuration validations. This is a correctness and
resource gate, never an index-convergence gate.

Per-seed `indices_summary.tsv.gz` contains N4096 five-repeat means, repeat range,
N2048→N4096 difference, bootstrap half-width, variance, group labels and
defined/passed/not-passed status. The full two-resolution repeated arrays remain
on HPC as `indices_repeats.npz`. Every completed domain gets S1/ST/S2 figures.
Ensemble reports first average five repeats per seed, then summarize the 500
endpoint means. Medians/quantiles use defined values and explicitly retain
undefined counts; no convergence filtering is applied. A companion S2 heatmap
shows the fraction of all 500 endpoints passing diagnostics. Global and
neighborhood outputs and ploidy/growth have separate color scales.

Low/high windows remain 0–1% and 3–5%. Process-group summaries are parameter
index means, not formal group variance decompositions; ST is not summed into
a variance fraction. Hypothesis-pattern fractions describe optimizer-endpoint
robustness, not biological confidence. All sources preserve raw negative S2
estimates. Selected response surfaces are conditional four-oxygen diagnostics
at the existing eight representative endpoints plus the global reference;
five scrambles apply to the Sobol designs, not duplicate deterministic surfaces.

Do not rerun submit_full.sh after jobs are recorded. Failed workers can be
resubmitted with the same frozen script and task IDs; it will inspect receipts
and checkpoints. Record retry IDs and correct dependencies for downstream
stages rather than silently aggregating missing blocks. An incomplete model
rectangle stops that domain's summary; diagnostic failure does not.

## Independent global post-processing

After all 150 global model blocks have completed, the full-bound domain can be
reported before the 500-neighborhood run finishes:

```bash
bash oxygen/code/O2_supply_demand_MAP/analysis/eFAST/sobol/submit_global_postprocess.sh
```

This submits one independent xxlarge/12h/1-CPU/16G job using the original
SALib 1.5.2 SIF, without a node restriction or a dependency on the neighborhood
array. The scheduler's user CPU quota still applies. Original model outputs,
fit inputs and the existing calculation snapshot are mounted read-only.
The existing evaluation and downstream job dependencies are preserved.

The result directory is:
`oxygen/results/eFAST/sobol/full_N4096/global/standalone_N4096/`.
It contains three PDF/PNG heatmaps (S1, ST, S2), N4096 five-repeat summaries,
all repeated estimates at N2048/N4096, 200-resample 95% bootstrap intervals,
numerical and repeat/resolution diagnostics, low/high-oxygen summaries,
interpretation, provenance and a checksum inventory. Raw model arrays stay
in the original HPC global model directory. Each output has an independent
color scale; negative sampling estimates and undefined values are preserved.

The global job uses a new frozen post-processing snapshot. The analysis
identity is distinct from each model-block identity, so interrupted analysis
can resume and completed analysis can be reused safely. Results are isolated
from the original full-run summaries. The job runs test_global.py in the
actual SIF before analyzing the submitted output.

After the global job is COMPLETED with exit 0:0 and its own status.json says
complete, collect its allowlisted deliverables:

```bash
python3 oxygen/code/O2_supply_demand_MAP/analysis/eFAST/sobol/collect_global.py
```

Collection verifies source inventory stability and every copied artifact's
SHA256. The local output path is the same repository-relative directory.

```bash
bash oxygen/code/O2_supply_demand_MAP/analysis/eFAST/sobol/submit_resolution_pilot.sh
```

Output defaults to `oxygen/results/eFAST/sobol/resolution_pilot/N2048_N4096`.
`SOBOL_REUSE_FROM` points to the completed N1024 pilot. Preparation checks all
model/helper and SIF hashes, case/input receipts and the exact nested parameter
prefix for each case and repeat. Workers recheck each source output, receipt
and parameter prefix, evaluate only the new 3072 vectors per matrix block,
then restore the complete N4096 output with explicit sample-ID offsets.
Recovery matrices for incremental points keep local IDs; add the recorded
`incremental_sample_id_offset` (1024) to map them to the complete design.

The design contains 22,118,400 operator evaluations, of which 5,529,600 are
verified N1024 results and 16,588,800 are new. Both outputs come from each
operator. There are still 1350 evaluation tasks; each evaluates 12,288 new
sample×oxygen points. xxlarge,12h, no specified node or array concurrency cap.
Nine case summaries, updated 91-pair heatmaps, 54 newly selected conditional
surfaces and a final report run after the evaluation array.

`resolution_comparison.tsv` aligns S1/ST/S2 values and repeat ranges across all
three resolutions. `reused_indices_validation.json` requires exact reproduction
of the old N1024 S1/ST/S2 means before response surfaces are launched. Per-cell
diagnostics use the actual N2048/N4096 ladder; raw negative values and undefined
cells remain visible. This is a representative resolution pilot, not the full
500-neighborhood or 201-oxygen analysis.

After the final report reaches COMPLETED/0:0, collect with the explicit roots:

```bash
python3 oxygen/code/O2_supply_demand_MAP/analysis/eFAST/sobol/collect_pilot.py \
  --remote-dir /share/lab_crd/taoli/Project/soft_couping_eFAST/oxygen/results/eFAST/sobol/resolution_pilot/N2048_N4096 \
  --local-dir oxygen/results/eFAST/sobol/resolution_pilot/N2048_N4096
```
# Resume the 355 unfinished N4096 neighborhoods

The October 18 recovery uses `recover_neighborhood.py`,
`run_neighborhood_recovery.sbatch`, and `submit_neighborhood_recovery.sh`.
It preserves the original full sampling plan and frozen model snapshot.
The 145 completed neighborhoods and completed global summaries are reused.

The authorized missing original model task IDs are 25318, 27551, 27552,
33636, 33637, 41044, 45450, 48408, 53810, 63791, 68684, 69826, 71559,
and 74956. All failed during Apptainer user lookup before model execution.
Only these 14 blocks are submitted again; existing valid model receipts
also make the original worker return without model evaluation.

Legacy statistical checkpoints store a B-block identity because the original
summary code overwrote its analysis-identity variable. The recovery verifies
the original plan, frozen code, runtime, unit designs, case bounds and B-block
identity before importing any certified cells. A complete temporary NPZ is
used only if it extends the durable checkpoint without conflicting values.
Original checkpoints remain untouched. Undefined indices at zero variance
remain certified results; convergence never excludes a seed or cell.

Remaining work is split into one seed/repeat and at most 16 unfinished oxygen
points, retaining both outputs, N2048/N4096 and the original 200-bootstrap
seed formula using the **global** oxygen index. Each analyzer call is saved
atomically. Workers lock their own shards, never a shared legacy checkpoint.
An oxygen-major float64 cache is constructed once per needed seed/repeat;
each source block is checksum/identity checked, and each consumed oxygen
slice is SHA256 checked. Cache creation only rearranges existing outputs.

Run the submission driver with `nohup` on RED after syncing the committed
branch. It first runs all tests and preparation in the original SIF on an
xxlarge compute allocation. After successful preparation, it submits only
the missing model blocks, ready statistical shards, affected shards waiting
for model repair, 355 case finalizers, a complete-input audit, aggregation,
54 response surfaces and final reporting. All jobs request xxlarge/12 hours;
no node or array concurrency cap is specified. Only the two obsolete pending
downstream jobs 20839426 and 20839427 are cancelled after replacement submission.

Default recovery directory:
`oxygen/results/eFAST/sobol/full_N4096/recovery_20261018/`.
Use `plan.json`, `missing_model_blocks.tsv`, `reused_cases.tsv`, `migrations/`,
`slurm/submission_jobs.tsv` and `slurm/launch_receipt.json` for provenance.
Raw model outputs, caches, repeated NPZ arrays and checkpoints remain on HPC.
Completed case tables/figures use their original `neighborhood10pct/full/seed*/`
locations. A final audit requires all 501 domains and unchanged reused cases,
not convergence. Volatile work after the last old checkpoint cannot be recovered.
