# Figure7 and supplementary Figure7: companion figure notes

## Authoritative update: 2026-10-24

The revised A/B share one daily viable-low-daughter denominator. For the
interval [d-1,d), A=100*C_HL/(C_HL+C_LL), B=100*C_LL/(C_HL+C_LL), classified
by maternal chromosome state immediately before division. Valid paired
shares sum to100%. Day0 and zero-output intervals are NA. Birth accounting
does not subtract death or maternal removal and does not count passaging.

The previous first-generation proliferation metric is now characterized in
supp_fig7_8, separately from permanent historical ancestry. The complement
of surviving first-generation direct-origin cells is "other low", which
includes descendants, NOT "never flux". At fixed current chromosome state,
origin does not confer an additional division-rate parameter. See the
fully revised Figure7_calculation_methods.md and _EN.md for all formulas,
simulation, distribution weights, daily counters and validation requirements.
Current culture replays the native fitting protocol and independently cycles
the last five complete measured control/O1/O2 tuples. Selection uses the full
positive-observation-day window nearest the target; composition is rescaled
to the next actual inoculum without a stock/N<1/unreached-target stop. The
1000-day adopted clock includes the recorded prefix. Summaries use three
deterministic protocol routes, not200 stochastic repeats. Conditional reseeding
is an extrapolation assumption, not proof of physical culture sustainability.
The historical _ver2 definitions, data and original protocol are unchanged. Production completion is
determined from the new dataset's completion receipt, not this note.

Proposed caption: "Expected viable low-ploidy daughter production in daily
intervals, partitioned by the mother's chromosome state at division. A,
high-ploidy mothers; B, low-ploidy mothers. Both counts use their combined
total as denominator. In vitro curves are equal-route medians and ribbons
are type-7 IQRs across three protocol routes, not confidence intervals;
in vivo retains the single continuous trajectory. Undefined intervals are omitted.
β=0 shows unblocked source accounting;
β>0 shows the specified pathway intervention, not an attribution of the
unperturbed model's growth."

## Superseded notes retained for the previously published data only

All references to "current B" below describe the pre-2026-10-23 publication,
not the newly implemented daily-source definition.

Current 2026-10-14 outputs use separate instantaneous and historical-lineage
versions. Figure7/supp_fig7_x A retains immediate high-to-low inflow, while B
now measures the current expected low-cell proliferation contribution of
currently-live first-generation daughters born directly from high mothers.
The direct-origin label is lost when the original cell divides or dies and
is never inherited by offspring. No age window is imposed and no death is
subtracted from proliferation. figure7_ver2/supp_fig7_x_ver2 A measures current
permanently-tagged low-cell number share, and its historical-lineage B is
unchanged. The current B is no longer equated to ver2 B. See
Figure7_versions_notes.md for the complete current definitions and diagnostics.
The B formula below also applies to the current publication. The culture
protocol notes describe the preserved October7 source simulations: the current
display-only amendment reuses those trajectories and draws no new passages.

Explanatory subtitles and footnotes are intentionally outside the exported
figures. Plot titles, axes, legends and categorical row/column labels remain.

## Blockade notation and source inputs

The Greek letter β denotes the blockade-strength coefficient. β=0 is no
blockade, β=0.5 redirects half of viable ordinary high-to-low daughters, and
β=1 redirects all such daughters to the mother's chromosome state. Viable
daughter yield is conserved. WGD, division frequency and direct death are
not directly altered by the intervention. In the repaired culture model,
population-dependent crowding responds to each trajectory and can therefore
change actual division activity indirectly. The implementation's frozen `q`
field is this same β coefficient;
only its visible notation changed, not the calculations.

C01 seed94 and C02 seed80 are the complete-parameter joint-fitting winners.
Both native environment scopes are retained. Tumor dynamics are continuous;
culture uses 200 canonical threshold-triggered stochastic passage repeats.
Culture lines/ribbons represent repeat medians/IQRs, not confidence intervals.
The horizon is 1,000 days.

The October7 culture repair reads each winner's authoritative
`fit_result.rds$ctx$invitro_cfg` and retains its complete 33-parameter culture
vector. Only the specified oxygen and native p_misseg counterfactuals change.
The saved DT=0.05-day Euler dynamics multiply division-linked birth/removal
and WGD by the configured crowding factor, while direct death is unscaled.
Culture passage times and samples are generated anew on the actual elapsed
clock under each blockade trajectory; historical passage events are not
replayed. β conditions use paired RNG streams, not enforced identical event
times. In-vivo rows are value-identical to the preceding numerical table.

## Figure7 and supp_fig7_5, _6, _7

Each A/B figure has four environment/start rows (in vivo 4N, in vivo 2N,
in vitro 4N, in vitro 2N) and C01/C02 columns. The high/low threshold is 66
chromosomes: low N<66; high N>=66.

| Figure | O2 (%) | p_misseg | A | B |
| --- | --- | --- | --- | --- |
| Figure7 | 0.5 | 0.1 | Flux percentage | Growth percentage |
| supp_fig7_5 | 0.5 | 0.3 | Flux percentage | Growth percentage |
| supp_fig7_6 | 5 | 0.1 | Flux percentage | Growth percentage |
| supp_fig7_7 | 5 | 0.3 | Flux percentage | Growth percentage |

A: `100 I_HL / (I_HL + I_LL)`, the fraction of newly produced viable low cells
from high mothers, not the fraction of the existing low-cell population with
high-ploidy ancestry. Missing complete-blockade curves can indicate no low
births: an undefined denominator is not zero.

B (current first-generation proliferation revision):
`100 sum(J[L,L] %*% D[L]) / sum(J[L,L] %*% x[L])`, where J is the
ordinary-division plus WGD viable-daughter production operator. Only mothers
currently in L contribute. The numerator counts viable low daughters produced
by direct-origin mothers currently in L; the denominator counts the same production
by all mothers currently in L. D labels directly born high-to-low daughters,
but never their low-mother-born offspring. The original cell exits D at its
first division or death. An earlier entrant that is still alive and has not
divided remains eligible, without an age window. The two production rates share the same low-population
size normalization, so their ratio is the displayed contribution.

Death, dividing-mother removal, high-mother inflow and passage dilution are
not subtracted or added. This is a proliferation contribution, not a signed
net-growth contribution, stock ancestry share, or blockade-induced growth
reduction. The natural finite values lie within 0-100%; they are not clipped.
Numerically negligible denominators are undefined, not assigned a percentage.
The preserved in-vivo calculation masks low mass <=1e-14 or low birth rate
per low cell <=1e-8/day. Repaired culture flux/birth denominators are masked
at <=1e-14 per live cell per day.
Main B and B in supp_fig7_5/6/7, plus every supp_fig7_4 cell, use this definition
and a 0-100% axis. A retains its original source-flux definition unchanged.

When a pure 4N start has no initial low cells at the 66-chromosome threshold,
the first low cells are direct entrants. Their later low-born descendants
are not D, so current B is not automatically 100%. The unchanged ver2 B
retains permanent ancestry and can still be 100%. Complete blockade can
leave no low compartment, making B undefined rather than zero; where low
proliferation remains defined but no direct entrants exist, current B is zero.

## supp_fig7_1 and _2

Both families and environments use O2 anchors 0/5%, p_misseg=0.3, 4N start,
day1000, β=0. supp_fig7_1 displays the source/daughter production ledger for
ordinary division and WGD, per 1,000 live cells per day. supp_fig7_2 displays
signed high-to-low births, low-to-low births, dividing-mother removal and
direct-death terms for N<66, normalized per 1,000 low-ploidy cells per day
(not per 1,000 total live cells). Dashed lines show biological net change; actual
stochastic passage losses are saved separately in the passage-event data.
This signed budget is a separate diagnostic, not the denominator of the new
proliferation percentage. The October7 repair recomputes the culture source
ledger using each repeat's actual day1000 crowding factor, then averages
per-live-cell production. The low-cell-normalized culture budget is averaged
after calculation within each repeat. In-vivo ledger/budget rows are retained.

## supp_fig7_3 and _4

Each supplement has eight environment/family/start matrix pages. Threshold
rows are 2, 2.5, 3, 3.5 and 4N (44,55,66,77,88 chromosomes); oxygen columns
are 0,1,2,3,4,5% in vivo, additionally 20% in vitro. p_misseg=0.1. Each cell
contains β=0,0.5,1 time courses. _3 uses the A flux definition, _4 the B growth
definition. Proliferation cells share the 0-100% domain and never clip data.

At 20% culture oxygen, low-threshold flux percentages can oscillate with
passage. A C01 seed94 / 4N / threshold44 / β=0 spot check has 123 passages
per repeat, eight days apart, with zero sampled N<44 cells after every
passage. Between passages the canonical deterministic live-cell dynamics
produce fractional expected abundances. Percentages of such tiny compartments
must not be described as measurement-resolved flux estimates. The raw
time courses are retained without smoothing or a new cell-count filter.

## Interpretation limitations

The historical October5 net ratio could become large or negative when low
births, death and mother removal cancelled; it also counted tagged high
mothers' production into L. It is preserved as historical evidence, not the
new proliferation definition. Original v4/v5/v6 results remain immutable.
Current bookkeeping amendment: `figure7_direct_source_growth_20261014_v10`,
reusing the original total trajectories with a non-inherited direct-origin
label and no new stochastic draws. The v9 permanent-lineage B is superseded.
Source repaired run: `figure7_invitro_repair_20261007_v7`, with 170 culture
tasks and all historical in-vivo rows retained. No refitting or modification
of the canonical source model occurred. The isolated propagator is checked
against canonical live counts, saved experimental lineage replay, both
crowding modes, stochastic passage accounting and timestep refinement.
All 948 pilot checks and 10 complete-run numerical checks passed. Culture
inoculum 569799 and passage threshold 6180783 remain the approved data-derived
protocol; these are not seed-fitted parameters. Previous chromosome-grid and
interpretation limitations remain, and no direction of biological agreement
is imposed by the repair. Receipts and pre-repair exports are under
`audit/figure78_invitro_repair_20261007/`.
