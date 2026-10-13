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
