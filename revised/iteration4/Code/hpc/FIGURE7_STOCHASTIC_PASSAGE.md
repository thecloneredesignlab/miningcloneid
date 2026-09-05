# Figure7 v4: reproducible threshold-triggered stochastic passage

Scope: iteration4 Figure7 and related supplements only. Figure6 and manuscript
text are frozen. This replaces the historical v3 fixed-window selector; it is
not a claim of exact reproduction of that fitting protocol. All growth/death/
missegregation kinetics still come exclusively from the external read-only
soft_couping_org/O2_supply_demand_MAP model and the supplied latest fitting roots.

## Clock and sampling contract

The canonical six-lineage summaries supply only inoculum and target population.
Use round(inoculum)=569799 and ceiling(target)=6180783 cells. The old median
five-day duration is source metadata, not a scheduling constraint.

At each day d, propagate the previous post-passage population by expm(M*1 day).
The first integer day with total live cells >= target triggers passage. Use
that day's entire population, including any threshold overshoot. Never rewind
to a fractional crossing time, an earlier day or a nearest-target state. Draw
exactly the inoculum without replacement, resume from that sampled composition,
and continue the real clock to day 10000. Store post-passage mean ploidy at d.
No threshold crossing means continued culture, not protocol failure. Numerical
errors fail validation; there is no survivor-only ensemble mean.

expm populations are real-valued expected counts. At a passage let T=round(total).
For normalized p, integerize with one uniform U and cumulative systematic
rounding: count_i=floor(T*cumsum(p)_i+U)-floor(T*cumsum(p)_(i-1)+U), fixing the last
cumulative count to T. Counts are nonnegative and sum to T; their expectations
are T*p, including rare states. Sequential conditional rhyper draws implement
the exact multivariate hypergeometric distribution conditional on those counts.
Integerization is an explicit approximation of expected populations, not a
full stochastic birth/death simulation. No per-day rounding occurs in culture.

## RNG and averaging

Master seed 20260904; R L'Ecuyer-CMRG/Inversion/Rejection, pinned R 4.4.2 SIF.
The sorted condition key is (invitro, cluster, original optimizer seed number,
oxygen, p_misseg, initial ploidy). Every key owns a separate nextRNGStream stream;
replicate r owns its (r-1)th nextRNGSubStream. The catalog and RNG settings are
saved with results. Keys never use worker id, task scheduling, chunk size or
run id. The first R replicas are unchanged when increasing R. Parameter-identity
deduplication shares the operator only: every original optimizer endpoint keeps
independent sampling streams, including the 35-fold identical group.

Average replicas within each endpoint, then average the same 50 q10 endpoints
equally. Save within-endpoint variance, Monte Carlo SE of the 50-endpoint mean,
and between-endpoint-mean variance separately. A 20/50/100-repeat full-time pilot
selects a fixed R; all final grid means must have MCSE <=0.01N before publication.
Pilot convergence is evidence at selected conditions, not a proof for the grid.
The pilot screening SE is the largest observed endpoint SE divided by sqrt(50),
assuming equally large variance in all 50 independent endpoints. It does not
apply the 0.01N target to an individual endpoint. If 100 repeats are insufficient,
the campaign stops for a compute-strategy decision before a larger full run.

Daily state, log count, RNG states, accumulators, current endpoint and oxygen are
checkpointed at 1000-day boundaries. Completed operators are atomically reduced
in deterministic order; their temporary checkpoint is removed only afterwards.
Resume requires identical source/model/config fingerprints. The no-crossing
optimization uses one deterministic scan: no random event implies all replicas
are exactly equal and no random numbers are consumed. It does not discard seeds.

## Data and figures

Daily 0:10000; initial 2:6N; fixed p_misseg .005,.01,.1,.2,.3; in-vitro oxygen
0:20 by .1. In-vivo continuous data remain oxygen 0:5 by .025; no in-vivo passage
variant exists. Main Figure7B uses the existing display subset and layout;
Supp7-12 shows full-range passage; Supp7-13 diagnoses trajectories, integer-day
cell counts and stochastic variability. A and continuous arrays are preserved
after model, parameter, multiplicity and grid validation. Old passage caches are
never reused. Historical runs/publication files are retained in versioned roots.

## Execution

Use run_figure7_stochastic_campaign_hpc.sh RUN_ID BASELINE_RUN on hpctpa3pc0028.
It runs the pilot, locks R, computes all new passage data, renders headlessly,
validates whole PDF words/bounds and seals a positive-list SHA256 return manifest.
60 worker processes, BLAS one thread; checkpoint batching stays within the user
allocation of 64 CPUs/512GB. Temporary/cache/output writes stay in iteration4.
The immutable source version is committed/pushed before HPC pull and launch.

Validation covers an independent daily R implementation; integer-day ceiling
and overshoot; slow/no growth; count conservation; hypergeometric mean/variance;
rare-state rounding; identical reruns; worker permutation; parallel equality;
bitwise checkpoint continuation; and deterministic scalar-dilution equivalence.
Actual external model matrices are tested on HPC. Final visual inspection and
checksum-verified return to local iteration4 are required before handoff.
