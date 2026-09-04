# Figure 7 segment-selection revision

Scope: Figure 7 and its supplementary figures only. Figure 6 and manuscript
text remain unchanged. The external soft_couping_org model is read-only.

## Computation contract

In vivo remains continuous natural growth. In vitro continuous and passage
are distinct propagations. Within each canonical culture segment, the passage
branch evaluates every positive integer day, selects the state whose live-cell
count has the smallest absolute distance to the target, and uses its normalized
composition to seed the next segment. Ties select the earliest candidate day.
The experimental clock advances by the segment duration, not the selected day.
At a segment boundary the stored mean is the reseeded selected-state mean
(right-continuous); any difference from the unselected segment-end mean is
recorded, not interpolated away.

The user additionally confirmed the current HPC inoculum-supply gate. A target-
selected population smaller than the required inoculum is replaced by the
smallest eligible positive-day population, if available; otherwise the protocol
is infeasible and that endpoint stops. Every subsequent ensemble cell requires
the full original weight of 50; it is NA/gray when any endpoint is infeasible.
Feasible endpoint weights and first failure days are retained separately. No
missing cell is silently interpreted as ploidy zero, and no survivor-only mean
is substituted. The HPC selector differs from the older local model helper;
its source checksum is recorded and the full current HPC implementation is
the validation authority.

Canonical conditions preserve the approved six-lineage equal-weight summaries:
geometric-mean inoculum and expansion ratio; median representative duration
(currently 5 days). This is the fitting **selection algorithm** under unified
conditions, not replay of the original six experimental schedules. Fixed O2,
fixed p_misseg, pure initial 2N-6N populations and the expm propagator retain
the existing Figure 7 intervention definition. The unmodified external
`ivt_extract_passage_end_state` is the selection oracle in the independent R
reference implementation. This validation does not claim expm and the fitting
Euler simulator are bitwise identical numerical solvers.

Both branches retain all 50 q10 optimizer endpoints via parameter-identity
deduplication and multiplicity-weighted arithmetic means. Data grids are daily
0-10000; in vivo oxygen 0-5 by .025; in vitro oxygen 0-20 by .1. Fixed p_misseg
values are .005, .01, .1, .2 and .3. Initial ploidy is 2,3,4,5,6N.

Continuous caches may be reused only after model-source, endpoint signature,
multiplicity, task-grid and finite-output checks. The old scalar-reset passage
arrays are never reused. New run ids and a versioned profile isolate them.

## Validation

`segment_selector_validation.tsv` compares actual C++ output with independently
propagated R states passed through the external fitting selector. It checks
every selected day, all daily means and final state vectors, using synthetic
edge cases and both actual in-vitro clusters at oxygen 0,.5,2,20 and p .005,.3.
`passage_vs_continuous_validation.tsv` measures actual array differences; it
does not require a nonzero difference and never hardcodes equality.

The main A/B layout uses 44 identical 38-mm square data panels with four aligned
rows, a shared 1-7N log-color scale and unchanged arithmetic averaging. Panel A
shows oxygen 0-2 on y and log10 p_misseg on x. Panel B shows daily 0-1000-day,
oxygen 0-2 slices, with initial 4N above 2N. The full-range A moves to Supp7-8;
inverse response and three full-range time diagnostics occupy Supp7-9 to 12.
The entrypoint map is `Code/config/manifests/figure7_entrypoints.tsv`; the shared
manifest with unrelated uncommitted Figure 6 changes is deliberately untouched.

Headless HPC runner: `Code/hpc/run_figure7_full_range_hpc.sh`, fixed node
hpctpa3pc0028 and the supplied R 4.4.2 SIF. Use 60 workers and one BLAS thread
per worker. Compute and draw stages can be run separately using `--compute-only`
and `--draw-only`; do not modify a running shell script. Archive previous
publication files before renumbering. Verify checksums when copying outputs
back to the local iteration4.
