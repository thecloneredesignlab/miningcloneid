# Figure7 instantaneous and historical-lineage versions, 2026-10-14

## Supersession notice, 2026-10-23

The section below describes the previous publication, NOT the revised main
Figure7 B. New A/B are complementary daily viable-low-daughter shares by
high versus low maternal state, with one combined denominator. The former
first-generation metric is characterized in supp_fig7_8 alongside historical
lineage and current-cell composition. The _ver2 definitions do not change.
For current formulas, parameter sources, simulation, weights, NA handling,
and caption text, use Figure7_calculation_methods.md or its English _EN.md
counterpart. Actual full-run completion requires a validated data receipt.

## Prior-publication reference only

## Current publication: first-generation direct-source proliferation

Figure7 B, supp_fig7_4 and the B panels of supp_fig7_5/6/7 now use currently
live first-generation low-ploidy daughters born directly from high-ploidy
mothers. These original entrant cells retain their direct-origin label only
until their own first division or death. Their offspring never inherit it.
Their current viable low-daughter production counts as their proliferation
contribution, but those daughters' subsequent proliferation does not belong
to the direct-origin source group. No entrant-age window is imposed: an
earlier entrant that is still alive and has not divided remains eligible.
All initial cells are unlabelled.

The current B ratio is 100 sum(J[L,L] D[L])/sum(J[L,L] X[L]). D has only
high-to-low birth input and original-cell division/death loss, with no
J[L,L] D[L] daughter-addition term. This is a non-inherited birth-origin
label, unlike the permanent historical ancestry label T retained by ver2.
Neither numerator nor denominator subtracts direct death, maternal removal
or passage dilution from viable low-cell proliferation. These processes
still affect the current live source population normally.

The original total-population trajectories and saved stochastic passages
are reused. At each saved culture passage, direct-cell retention is its
conditional expectation given the saved within-chromosome-state total-cell
sample. No new passage or origin-label RNG draws are made. Culture medians
and IQRs summarize per-repeat conditional expected contributions over the
existing phenotype/passage trajectories; they are not confidence intervals.
Current B need not equal A or ver2 B, nor is it constrained to be far below
100%. Undefined denominators remain missing rather than zero.

Run: figure7_direct_source_growth_20261014_v10. Receipts and recoverable
pre-amendment exports are under audit/figure7_direct_source_growth_20261014.
A and all ver2/other figures are protected. No fits, native parameters,
biological kernels, oxygen/blockade conditions or total-cell passage samples
are changed. The earlier definitions below are archived, not current B.

## Archived v9 amendment: current historical-lineage proliferation

The superseded v9 instantaneous-version Figure7 B, supp_fig7_4, and the B panels of
supp_fig7_5/6/7 use the contribution of currently-low cells with historical
high-to-low ancestry and their descendants to current expected low-cell
proliferation. The ratio is 100 sum(J[L,L] T[L])/sum(J[L,L] X[L]). Only
currently-low mothers and viable low daughters contribute; no direct death,
maternal removal, passage dilution, first high-to-low inflow, or cumulative
past proliferation is subtracted or added. No entrant-age window is imposed.
The current B values and IQRs equal the corresponding ver2 B values by design.
The historical label defines the current source population, not a cumulative
growth statistic. A remains the original instantaneous inflow ratio. All
ver2 figures, supplementary diagnostics 1/2, flux matrices 3, and Figure8
remain unchanged. Undefined denominators remain missing, not zero.

This display-only amendment uses run figure7_current_growth_20261014_v9 and
audit/figure7_current_growth_20261014. The original v8 replay, parameter
vectors, stochastic passages, and numerical receipts are immutable; no new
simulations or stochastic draws are needed. The following text records the
superseded v8 instantaneous B definition, not the current publication B.

## Archived v8 definition

The user confirmed two separate versions and explicitly accepted identical
instantaneous A/B formulas. Figure7 and supp_fig7_3/4/5/6/7 now use immediate
viable high-to-low inflow divided by all viable low-cell production (H-to-L
plus L-to-L). B is an instantaneous inflow contribution to positive low-cell
production, not proliferation of an entrant-age cohort. A/B curves and IQRs
therefore coincide exactly; this is not independent evidence.

figure7_ver2 and corresponding supp_fig7_x_ver2 use a historical transition
label. All initial cells are untagged. An ordinary high-to-low daughter and
all its descendants are tagged permanently, including descendants of an
initial 2N population that went high and subsequently returned low. This
label is not initial 4N identity. Low/high membership always refers to the
current chromosome state, using N<threshold and N>=threshold.

Historical A is 100 times the tagged cell count currently in L divided by
the total current L count. It is a surviving-population ancestry share, not
a transition rate, integrated flux, or the tag fraction of the whole live
population. Historical B is 100 sum(J[L,L] T[L])/sum(J[L,L] X[L]); only
currently-low mothers and viable low daughters contribute. J includes
ordinary division and WGD daughter production, without death, maternal
removal, high-mother first-tag inflow, or passage dilution in this ratio.
This B reproduces the October6/7 proliferation definition. Undefined
denominators are missing, not zero. True percentages are not clipped.

Both versions retain complete C01 seed94 and C02 seed80 context-specific
parameters, pure 2N/4N starts, beta=0/0.5/1, all existing imposed O2 and
p_misseg counterfactuals, and the original chromosome grid. Tumor trajectories
are continuous. Culture uses fit_result.rds$ctx$invitro_cfg, its DT/crowding,
and the existing stochastic threshold passages. No fits, model kernels, or
new stochastic passage samples are generated. Replay uses saved X/U passage
anchors and must reproduce old daily flux, ploidy, tag fraction, final states,
and historical B before figures are published. Quantiles are computed within
each condition across the executed culture repeats, with type-7 median/IQR;
nominal 200 repeats are retained. A proven deterministic no-passage case can
have one executed repeat, as in the authoritative source run.

Instantaneous supp_fig7_1/2 remain the original production ledger and signed
low-population budget, with their original units. They are not percentage
curves or denominators for B. Their ver2 diagnostics partition production and
budgets by daughter lineage. Tagged ordinary production is J*T+K_HL*U;
never-tagged production is (J-K_HL)*U. First high-to-low daughters of never-tagged
high mothers therefore belong to the tagged child group. WGD daughters inherit
their mother's label. Both parts sum to the original production ledger.
In the low budget all H-to-L inflow goes to T, while L-to-L births, mother
removal and direct death are split by their respective current T/U populations.
Both groups share the original total-low-cell normalization, so their budget
components sum to the original diagnostic. Direct death is not crowding-scaled.

The 66-chromosome threshold and O2/p_misseg combinations remain 0.5/0.1 (main),
0.5/0.3 (supp5), 5/0.1 (supp6), and 5/0.3 (supp7). Supplementary matrices retain
the five threshold rows, six tumor or seven culture oxygen columns, and eight
environment/family/start pages. The blockade label is Greek beta, not q.
Explanatory methods remain outside the exported graphs.

Run: figure7_instant_history_20261014_v8. Numerical inputs are the immutable
figure7_flux_full_20261001_v4 tumor trajectories and corrected
figure7_invitro_repair_20261007_v7 culture trajectories. Backups and QA are
under audit/figure7_instant_history_20261014. Figure8 and all non-Figure7
artifacts are protected by SHA256. This extension does not establish
biological lineage replacement or individual-cell death from a share alone.
