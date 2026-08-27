# Model-fix assessment

## Root-cause finding

The implementation findings agree with both code-review documents.

1. The previous extraction function chose the simulated day whose live-cell
   count was closest to the observed harvest, and that selected state fed
   growth, flow/karyotype predictions, and the next passage.
2. The previous propagation reconstructed the next state from a selected
   chromosome fraction and unconditionally scaled it to the observed next
   inoculum, including upward scaling when predicted cells were insufficient.
3. The previous job adapter collapsed multiple O1/O2 records in a segment with
   medians and one shared segment state, so parallel experimental lineages were
   not independent state chains.

## Implemented state contract

For passage \(p\) with recorded duration \(T_p\), all fitted modalities now use
the same fixed endpoint:

\[
X_{\mathrm{end},p}=X_p(T_p), \qquad
N_{\mathrm{end},p}=\sum_i X_{\mathrm{end},p,i}
\]

\[
g_{\mathrm{pred},p}
=
\frac{\log N_{\mathrm{end},p}-\log N_p(0)}{T_p}
\]

`selected_day` is retained only as a compatibility alias equal to
`endpoint_day = passage_duration`. The old argmin is exported separately as
`closest_day_diagnostic` and cannot select a state.

At a nonterminal boundary with next recorded inoculum \(N_{\mathrm{required}}\):

\[
X_{p+1}(0)=
\begin{cases}
X_{\mathrm{end},p}
\dfrac{N_{\mathrm{required}}}{N_{\mathrm{end},p}},
&N_{\mathrm{end},p}\ge N_{\mathrm{required}}\\[8pt]
X_{\mathrm{end},p},
&N_{\mathrm{end},p}<N_{\mathrm{required}}
\end{cases}
\]

Thus `boundary_scale` is either
\(N_{\mathrm{required}}/N_{\mathrm{end},p}\le1\) or exactly 1. No
cell-insufficiency objective penalty was added.

The adapter now constructs six deterministic chains:

- 2N-C, 2N-O1, 2N-O2
- 4N-C, 4N-O1, 4N-O2

Every passage keeps its own duration, initial count, final count, observation
ID, parent, and state history. The objective passes one shared `run_params`
object to both cohorts and all six scenarios; no lineage parameter, random
effect, or stochastic clone parameter was introduced.

## Likelihood unit and weighting

Each modality asserts that `passage_id` is unique before aggregation.
Aggregation is:

1. single-cell mean (karyotype) or single-passage/sample likelihood;
2. mean within each lineage/scenario;
3. equal mean across available lineages within each 2N/4N cohort;
4. equal mean across the two cohorts.

This prevents a lineage with more passages or single-cell observations from
being treated as more independent cultures. Only the true depth-zero
karyotype anchors are evaluated once under an explicit `INITIAL` group rather
than copied into O1 and O2. Non-root landmark records are evaluated once
against the unique formal source segment state: `SUM-159_NLS_2N_A6M_seed`
against `2N-C-A1`, and `SUM-159_NLS_4N_A4M_seed` against `4N-C-A2`.

Custom observation grids are filtered to finite days in
\([0,\mathrm{passage\_duration}]\), with zero and the exact endpoint always
present. The segment runner independently rejects any grid that could extend
the C++ simulation beyond the recorded passage endpoint.

## Verification evidence

The targeted structural/migration tests and the entire oxygen test suite pass.
The full suite has two pre-existing data-availability skips: one local seed1
report mirror and one materialized multi-warmup regression result.

The current non-overwriting replay is stored under
`results_review_fixes/`. The seed10 and seed340 parameter replays both produce:

- 114 unique passage predictions across six scenarios;
- exact expected cumulative passage/endpoint times;
- 114/12/20 unique growth/karyotype/flow likelihood units;
- no boundary scale greater than 1;
- exact shared parameters across all scenarios.

Seed340 exercises the insufficient-state branch at 33 boundaries. All 33 have
`cell_number_after == cell_number_before`, `boundary_scale == 1`, and zero
component-wise difference between the 133-dimensional endpoint and next-passage
state vectors. Its other 75 nonterminal boundaries are exact proportional
downsamples.

## Remaining risks and decision

Code-structure decision: **GO for a separately authorized full refit**.

Scientific-result decision: **NO-GO for interpreting the old best-fit
parameters as final results**. The structural change alters simulated histories,
likelihood weighting, and objectives; seed10/seed340 here are diagnostic
replays, not reoptimized fits. A fresh multi-seed optimization and post-fit
quality review are required before replacing published or manuscript-facing
results.

Additional residual risks:

- the lineage identity parser intentionally relies on the current canonical
  `_<2N|4N>_<C|O1|O2>_A<number>_seed` data-ID contract;
- `selected_day` remains in output solely for consumer compatibility, so future
  code should prefer `endpoint_day`;
- severe predicted depletion is now exposed rather than hidden by artificial
  replenishment, which may produce long low-cell trajectories and materially
  different optima;
- this deterministic mean-field replay does not establish biological
  replication or causal interpretation of fitted latent mechanisms;
- the inaccessible `/share/.../fit_invitro_O2_buffering_500seed` baseline was
  not directly compared or modified.
