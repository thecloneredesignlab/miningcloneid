# Figure7: native-fit passage replay and 1000-day extrapolation

This Figure7-only implementation is at the small-acceptance stage. It does not
change the formal Figure7 data, rendering entry, Figure8, historical `_ver2`
products, shared passage helper, or locked biological model. No full campaign or
redrawing is authorized by this entry.

## Scientific contract

Read all 33 in-vitro fitted parameters from the approved C01 seed94 / C02 seed80
joint winners. Read culture configuration from `fit_result.rds$ctx$invitro_cfg`,
not the in-vivo `fit_config.rds`. Change scenario O2 and mechanistic `p_misseg`
(0.1/0.3), chromosome-grid upper bound, blocking beta, and the explicitly pure
44/88-chromosome initial condition only. Effective missegregation probability is
still computed by the original oxygen/ploidy-dependent model. Pure initial
conditions are used only at the true root; subsequent composition is inherited.

The recorded experimental adapter DAG is replayed with its declared parentage,
native previous-segment fallback, actual inocula/targets, observation windows,
and native target fallback. All original branches, including short reoxygenation
branches, are retained in prefix diagnostics. Scenario O2 is held constant from
day zero: the observed oxygen sequence is not replayed.

For each initial cohort, identify the longest nominal-duration terminal path in
each of the control, O1 and O2 groups (lexical terminal ID breaks ties). Each
long-term route then cycles its own last five **raw complete measured** segment
tuples `(initial_cells, final_cells, obs_days_local, duration_days)` in their
original order. Missing measurements do not become complete because a fallback
exists; target greater than inoculum is **not** an eligibility requirement.
Three- and ten-segment cycles are sensitivity scenarios, not extra replicates.
Short reoxygenation offshoots are not long-term cycling templates.

The extension is an explicit protocol assumption beyond the experiments, not an
observed 1000-day culture schedule.

## Passage and elapsed time

Simulate each tuple over its full native candidate window using the original
Euler step and logistic crowding. The unchanged native selector chooses the
positive observation day whose expected live-cell count is closest to the
tuple's target (earlier candidate wins ties; absent target uses the last positive
candidate). This is **not** threshold-triggered stochastic passage.

Given selected composition `f_k` at adopted time `tau_k`, the next segment starts
from `I_(k+1) * f_k`. There is no stock-feasibility gate, no `N < 1` stop, no
`target not reached` stop, no added random sampling, and no standardized common
569799-cell inoculum. Record the actual scaling factor, including values above
one. This is conditional expected reseeding, not inventory-conserving evidence
that such a culture could be maintained physically.

Total elapsed time is the sum of **selected** segment durations, including the
experimental prefix; discarded candidate tails do not advance time. Continue to
1000 days from pure day zero, not 1000 days after the experiment. For the final
segment, perform selection using the full candidate window **before** truncating
the retained path at day 1000. Do not force early passage. At an exact final
passage boundary, keep the pre-passage endpoint and do not apply an unused reseed.
Blocking conditions select their own passage times using the same algorithm;
differences in passage calendars must be reported.

## Accounting and origin diagnostics

Daily window `[d-1,d)` counts each new viable low daughter once by the mother's
chromosome state immediately before division (`low < threshold`). A is
`100*C_HL/(C_HL+C_LL)`, B is `100*C_LL/(C_HL+C_LL)`. Ordinary/WGD daughter
multiplicities come from the original kernel. Neither deaths nor maternal
removal are subtracted. Passage is not a birth. Both shares are NA at day zero
or when the combined production is zero.

First-generation direct H-to-L labels D are not inherited: a labeled cell exits
this group at its own first division or death. Never-H-to-L lineage U follows
`U' = [c*(G-K_HL)-diag(mu)] U`. D follows
`D' = c*K_HL*X - (c*lambda+mu)*D`. The current low groups are X, D, X-D, X-U and U;
rates are evaluated at their current chromosome states and contemporaneous O2
and crowding, not at the former high-ploidy mother's state. X is not altered by
these passive labels. Low-group D/U distributions are retained for the current
origin diagnostic, not promoted to historical `_ver2` products.

Normalization of numerical states does not discard absolute population size:
`log_cells` tracks that size and daily expected daughter counts are recovered
from the corresponding day's **post-reseeding** start size. Small expected
counts remain legal. Numerically vanishing source-group conditional means need
separate caution; meaningful-source CDF comparisons use a reported mass threshold
of `1e-12`, not a population-stop or exclusion rule. Native uniform fallback at
zero selected mass is retained, flagged, and invalidates unidentified ancestry
metrics rather than inventing labels. Its occurrence fails numerical acceptance
but does not introduce a stock-feasibility stop.

## Summaries and acceptance scope

Keep control/O1/O2 outputs separate. Summary curves use equal-route medians and
type-7 25th/75th percentiles. This IQR describes **protocol-scenario spread**, not
random passage replicates, optimizer-seed uncertainty or a confidence interval.
No deterministic trajectory is duplicated to manufacture 200 repeats.

The dedicated entry is:

```sh
Rscript --vanilla Code/Figures/run_Figure7_passaged_acceptance.R \
  --audit-root=/absolute/path/to/a/new/backed-up/iteration5/audit/directory
```

It requires a verified backup receipt and refuses to overwrite a previous
execution. The same small suite can continue with `--resume=yes` after its local
safety timeout. Reuse requires exact core/config/test hashes, unchanged raw inputs
and model lock, all three complete 1000-day routes, and matching registered RDS
SHA-256 values. Cached 1000-day paths are not recalculated; finite prefix and
independent checks are replayed. Unregistered interrupted writes are moved to a
recoverable audit subdirectory. Source snapshots and stage receipts are retained.
It uses current raw fits and the private locked model, never old audit
results or archived predictions. It performs a bounded small suite, not the
Figure7 full task manifest: four representative conditions at 137/139/141N,
coverage at 141N, blocking/oxygen/threshold probes, additional unblocked oxygen/p
crossing probes (including the low-O2/high-p blocking baseline), and 3/10-segment
template sensitivity. Results and compiled caches go only into that audit.

Grid acceptance requires two adjacent passes on complete 1000-day paths for A/B,
mean ploidy, log population, full-state CDF, upper-bound loss and identical passage
calendars. Thresholds are in `Code/config/figure7_passaged.json`. Independent
absolute-count Euler checks test source counters and labels; original-kernel
checkpoints test actual root and late-extension states. A finite-experiment grid
certificate does not imply a 1000-day certificate or acceptance of every full
Figure7 condition. Template sensitivity is reported separately, not called grid
error. No figures are created by these tests.
