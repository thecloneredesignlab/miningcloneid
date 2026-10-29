# Figure7: daily daughter-source accounting and origin-population analysis

Updated 2026-10-24. This is the English counterpart of `Figure7_calculation_methods.md`. Current in-vitro prediction uses native-fit passage replay and last-five-segment cycling. Historical `_ver2` data and their original protocol remain unchanged. Numerical completion requires the dataset's `provenance/completion.tsv`; publication additionally requires the figure acceptance/promotion record. This document alone does not establish that production images have been updated. Previous methods and outputs were backed up first.

## 1. Three distinct questions

1. Figure7 A/B partition newly produced viable low-ploidy daughters within the same daily interval by their mother's high versus low chromosome state immediately before division. These are complementary production shares.
2. `supp_fig7_8` characterizes cells currently in the low compartment, separating first-generation direct origin from permanent historical ancestry. Cell stock, proliferation output, per-cell rates, and ploidy distributions are different estimands.
3. `_ver2` retains its historical-lineage current-cell share and current proliferation share, not a daily birth-source partition.

None of these is a net-growth contribution, the causal growth loss caused by blockade, or simply a fraction of maternal division events. The categories are parental source, not “division versus missegregation”: both parental compartments can divide and missegregate.

## 2. Parameters, locked model, and experimental conditions

Each family uses one complete joint-fitting winner, without refitting or pooling optimizer seeds:

| Family | Joint-fit directory | Winner seed | in vivo / in vitro parameters |
| --- | --- | ---: | --- |
| C01 | fit_joint_tsne_vi_seed366_C01_vt_seed144 | 94 | 26 / 33 |
| C02 | fit_joint_tsne_vi_seed25_C02_vt_seed144 | 80 | 26 / 33 |

The fitting root is configured in `Code/config/figure78_inputs.json`. Complete natural-scale parameters come from `joint_best_params_long.tsv`: `shared_invivo` and `invitro_effective`, respectively, not only 14 coupling parameters. In vivo uses `fit_config.rds`; culture uses the authoritative `fit_result.rds$ctx$invitro_cfg`, preserving DT, Crowding, crowding mode and K.

Only the SHA256-locked private model under `Code/vendor/figure78_model` is used. The complete entry is `Code/Figures/run_Figure7.R`, dispatched by `figure7_current_pipeline.R`; current-only rendering is `draw_Figure7.R --analysis=instantaneous`. Retired workers, archived predictions and Figure8 predictions are not simulation inputs.

Main Figure7: O₂=0.5%, native p_misseg=0.1. Supplements5/6/7 use .5%/.3, 5%/.1, and 5%/.3. Each includes both environments, both families and pure 2N/4N starts. N-unit22; the fitting grid and unchanged continuous in-vivo grid are22–154 chromosomes. Current culture prediction compares upper bounds137/139/141N and accepts141N (22–3102 chromosomes) only after all convergence gates pass. This is an explicit prediction-grid setting, not a replacement of fitted parameters. L means N<66, H means N≥66. β=0, .5, 1.

Changing native p_misseg does not fix the effective state/oxygen-dependent missegregation probability to that same number. Fitted p_mis_base, k_o_mis and the canonical transformation remain active. Complete parameter passthrough does not imply every parameter affects an imposed-oxygen, pure-start diagnostic.

Matrices supp_fig7_3/_4 use thresholds44,55,66,77,88, in-vivo oxygen0,1,2,3,4,5%, and an additional20% in vitro, with native p_misseg=.1. L/H are redefined by each subplot's threshold.

## 3. Canonical viable-daughter operator

Let x_j(t) be the live count of maternal chromosome state j; λ_j and μ_j are its current division and direct-death rates. B_ij already contains the expected viable ordinary daughters in state i per event, including native multiplicities, boundary losses and viability filters. W is the original WGD output operator.

\[
J_{ij}=\lambda_j[(1-p_{wgd})B_{ij}+p_{wgd}W_{ij}],\quad
G=J-\operatorname{diag}(\lambda),\quad M=G-\operatorname{diag}(\mu).
\]

Do not multiply B by an additional factor of2. Preserve the canonical WGD multiplicity. Count each expected viable low daughter once, irrespective of subsequent death. Daughter output is not maternal event count.

Blockade redirects the fraction β of ordinary viable H→L output to the mother's chromosome state, conserving each mother's total viable yield. β=0 is unblocked; β=1 completely blocks this pathway. Native WGD, λ and μ do not change directly; culture density feedback can change actual activity indirectly. The compatible field `q` means β, not Figure8's mixture ratio Q.

## 4. Revised A/B: identical daily denominator

Let c(t)=1 in vivo or the canonical culture crowding factor in vitro. A curve point labeled day d refers to [d−1,d):

\[
C_{HL}(d)=\int_{d-1}^{d}c(t)\sum_{j\in H}\sum_{i\in L}J_{ij}x_j(t)\,dt,
\quad
C_{LL}(d)=\int_{d-1}^{d}c(t)\sum_{j\in L}\sum_{i\in L}J_{ij}x_j(t)\,dt.
\]

\[
A(d)=100C_{HL}/(C_{HL}+C_{LL}),\qquad
B(d)=100C_{LL}/(C_{HL}+C_{LL}).
\]

Classify by the mother's chromosome state immediately before division, not ancestry, starting cohort or a future state. A high-mother-born low daughter that divides later in the same day produces LL daughters: she is now a low mother. Every birth is counted once. For nonzero denominator, each matched repeat satisfies A+B=100%.

Do not subtract death, maternal removal, exit from L or passage loss. Those affect population stocks, not the partition of births. If no viable low daughters are produced, both shares are NA. Day0 has no completed interval and is NA, not0 or100.

In vivo, an augmented matrix exponential exactly integrates the original continuous model's daily production. It is not an endpoint rate fraction. Culture accumulates `DT × c × J × x` at the same DT=.05 Euler substep starts used by the fitted model. Counters use daily-start live population as a common scale, saving its log count; each substep is weighted by its actual population expansion, not equally by normalized fractions.

A culture passage occurs on the positive observation day selected by the original native algorithm, after that day's growth and births. The selected composition is multiplied by the next tuple's actual inoculum. Reseeding contributes neither HL nor LL births. A boundary is reported once as the pre-passage endpoint; the next day's birth count starts from the reseeded inoculum. In vivo has no passage.

### Small numerical example

An independent two-state test starts entirely high. HL output rate is .4/day and LL output rate2/day, with DT=.5. Unnormalized `(L,H)` states are `(0,1)→(.2,1.3)→(.56,1.69)`. Daily counts per initial live cell are HL=.46 and LL=.20, denominator .66; A=69.69697%, B=30.30303%. A passage at the end of the day changes neither birth count. This is an accounting test, not a replacement for the actual model parameters or numerical results.

## 5. Current culture protocol, extrapolation and route summaries

Preserve the fitted culture configuration: DT=.05 days, crowding enabled, logistic mode and K=10¹². Read actual inocula, targets, complete observation windows and parentage from each family's winning-fit 2N/4N experimental adapter. Current culture does not use a common569,799-cell inoculum, a6,180,783-cell trigger or stochastic passage. Those are older-result rules; the unchanged continuous in-vivo adapter retains its old numerical count scale without passing it to culture prediction.

First propagate the full original candidate window. The unchanged native selector chooses the positive observation day with expected live count nearest the target, preserving its tie/fallback rules. Normalize the selected composition and rescale it to the next tuple's actual inoculum. Insufficient stock, expected population below1 or an unreached target adds no stop or reseeding cap. Record scaling factors and shortages. Numerical normalization retains absolute expected size in log_cells. This is conditional reseeding, not inventory-conserving proof of physical culture sustainability.

Replay the recorded DAG from truly pure44/88-chromosome roots; short reoxygenation branches remain prefix diagnostics only. For control/O1/O2 separately, choose the longest nominal-duration terminal path and cycle its last five raw complete measured tuples in original order, preserving paired inoculum, target and window. Target greater than inoculum is not an eligibility requirement. Scenario oxygen is fixed from day0 rather than replaying measured oxygen changes. The1000-day clock includes the prefix and adds selected passage durations only. Select the final segment using its full candidate window before truncation; do not reseed at day1000 without a subsequent segment. Each β selects its own passage calendar; record differences.

Current culture has three deterministic protocol routes, not200 random repeats. Calculate daily shares and current-population measures within each route, then take equal-route medians and type-7 IQRs; do not pool absolute daughter counts across routes to change weighting. IQR denotes protocol-scenario spread, not optimizer-seed uncertainty, random-replicate uncertainty or a confidence interval. Last-three/ten cycles are separate small-acceptance sensitivity evidence, not extra summary replicates. D/U are passive expected labels rescaled with the same composition, without new RNG.

## 6. Two independent origin partitions in supp_fig7_8

Restrict all groups to cells currently in L. X is the total live state, D is surviving direct H→L daughters that have not yet divided, U is never-H→L ancestry, and T=X−U.

First-generation view: `all_low=X_L`, `direct_first_generation=D_L`, `other_low=X_L−D_L`. Other low includes descendants of D; it must not be called “originally low” or “never flux.”

Historical view: `all_low=X_L`, `historical_lineage=T_L`, `never_HL_lineage=U_L`. Ancestry is inherited permanently. With a pure4N start, all later low cells may have H→L ancestry and U_L can remain absent.

The direct cohort obeys:

\[
\dot D=cK_{HL}X-[c\operatorname{diag}(\lambda)+\operatorname{diag}(\mu)]D.
\]

New H→L daughters enter D; their first division or death removes them. Offspring never inherit D. There is no arbitrary age window. Use the daughter's current low-state λ and J, not her high mother's rates. The permanent untagged state obeys `U′=c(G−K_HL)U−μU`, with inherited tagged state T=X−U.

For any current-low group Z, calculate separately:

- cell share `100 ΣZ/ΣX_L`;
- low-daughter production share `100 Σ_j∈L[Σ_i∈L cJ_ij]Z_j / Σ_j∈L[Σ_i∈L cJ_ij]X_j`;
- per-cell expected low-daughter output `Σ_j∈L[Σ_i∈L cJ_ij]Z_j/ΣZ`, daughters per current cell per day;
- per-cell division rate `Σ_j∈L cλ_j Z_j/ΣZ`, day⁻¹;
- mean ploidy `Σ_j∈L(N_j/22)Z_j/ΣZ`.

Empty-group means and distributions are NA, not zero. Percentage denominators require nonempty total low population or positive total low-parent output, respectively. Within each view, the selected source group and its complement partition the corresponding total stock/output.

### 6.1 Reporting validity of conditional means: postprocessing only

The three conditional-mean time courses in `supp_fig7_8` (per-cell low-daughter output, division rate and ploidy) share one validity rule. For every saved route, day and current-low origin group g, let f_g be its mass fraction of **all currently viable cells**, rather than its fraction of the low compartment. Retain the means only when f_g≥10⁻⁸, ancestry is identifiable, and f_g and all three means are finite. **There is no absolute expected-count cutoff.** An expected count below one does not enter the mean-validity mask or route eligibility for the median/IQR.

The relative threshold is a numerical/reporting criterion, not an extinction criterion, passage-stop rule or experimental detection limit. The same group-level mask applies to all three means. Direct/historical source labels are initially empty and remain empty under β=1. A pure4N root has no never-H→L low lineage; all its low groups are also empty initially or under complete blocking. Means of these structurally empty groups are NA even when saved weights contain tiny positive roundoff residuals. A genuine zero division/output rate in an otherwise valid group remains zero; values are not selected by rate magnitude, curve shape or the presence of a spike.

Expected counts are reported separately using the saved log(N_total)+log(f_g), with `log_expected_cells`, `expected_cells_lt1` and count-status fields. Extreme totals are not exponentiated and no inoculum or counting scale is changed. Unavailable counts receive a separate flag and do not exclude otherwise reliable composition means; possible exponentiation underflow likewise introduces no cutoff. A below-one mean is a deterministic expected-mass-weighted ratio, not a measured average over fractional cells, and does not establish the existence or extinction of a physical source population.

Apply validity before route aggregation. All three control/O1/O2 routes must qualify for an in-vitro median and type-7 IQR. If fewer than three qualify, all three mean summaries and their q25/q75 are NA, with the actual eligible-route count recorded separately. The single continuous in-vivo route uses the same group rule. Missing intervals are not filled with zero, interpolated, smoothed or bridged by lines/ribbons.

Reporting-only sensitivity checks compare relative thresholds 10⁻⁹/10⁻⁸/10⁻⁷, with no absolute-count filtering branch. The publication threshold remains10⁻⁸; sensitivity reports changes in valid support. Original simulation tables and receipts remain immutable. The derived table, route validity, expected-count flags, exclusion reasons, temporal coverage and sensitivity outputs are retained in `data/Figures/Figure7/postprocess/population_validity/`, tied to the original scientific fingerprint and input SHA256 values and authenticated before plotting. Subsequent postprocessing/redraw uses the retained data in that current directory.

This rule affects only these three conditional-mean time courses in `supp_fig7_8`. Cell shares, low-daughter production shares, original CDF data, Figure7 A/B, Figure8 and historical `_ver2` products are unchanged. CDFs retain the definitions in section7; do not describe them as having received this mean-validity filter.

At the same family, environment, oxygen and chromosome state there is no origin-specific growth parameter. Mean rate differences reflect composition, crowding, state-dependent division and viable-daughter yields. Previous B>A does not establish faster source-cell growth because previous A was not D's current stock share. Compare D's output share with D's stock share, alongside its per-cell rates.

## 7. Distributions, summaries and presentation

At days0,10,100,500,1000, save cell-weighted ploidy, actual division-rate and expected low-daughter-output distributions. Normalize within each nonempty group/repeat, then average valid-repeat CDFs equally. Ribbons are the25th/75th percentiles of repeat CDFs, not the within-cell IQR and not a confidence interval. Ploidy uses exact chromosome lattice points; rate CDFs use101 common grid points from0 to the condition's original maximum rate, including each repeat's contemporaneous crowding factor. No KDE introduces fictitious states.

Time metrics are computed per protocol route and summarized by median and type7 q25/q75. The CDF repeats described above correspond to the three current culture routes, not200 draws. Do not pool all cells before calculating rates or confuse protocol IQR, within-cell distributions and optimizer-seed differences.

`supp_fig7_8.pdf` contains48 pages: two origin definitions × two oxygen anchors(.5,5) × two native p_misseg(.1,.3), each with three β=0/.5/1 time-profile pages and three unblocked β=0 distribution pages. Each page contains both environments/families/initializations. Unblocked accounting is explicitly separate from blockade comparisons. PNG previews page1; the PDF/page manifest is the complete product.

## 8. Unchanged historical and diagnostic definitions

Historical A=`100 ΣT_L/ΣX_L`; historical B=`100 LL(T_L)/LL(X_L)`, without subtracting death. Stock and production shares are not complementary and need not sum to100%. Existing historical formulas and passive replays remain unchanged.

supp_fig7_1 shows ordinary/WGD mother-daughter production matrices. supp_fig7_2 shows signed HL inflow, LL output, maternal removal and death budgets. These are not revised daily A/B percentage curves. Historical diagnostic supplements retain their lineage partitions and denominators.

## 9. Data, verification and publication boundaries

The single canonical directory is `data/Figures/Figure7`. Plot tables include `daily_source_counts_summary.tsv`, `low_population_time_summary.tsv` and `low_population_distribution_summary.tsv`. The isolated campaign's `routes/` RDS files retain daily states, HL/LL per-start-live-cell counts, log starting size, absolute expected births, paired shares and passage events, with SHA receipts and task fingerprints. Rendering does not repropagate trajectories. Independent counter-column medians must not be summed and treated as a matched-route ratio.

Compatible keys `flux_percent`/`growth_percent` in `instant_time_summary.tsv` now explicitly mean daily high-parent/low-parent daughter production, with `source_accounting=daily_viable_low_daughters`. Rendering refuses to relabel old first-generation-B data. `history_time_summary.tsv` retains historical definitions.

Back up first, run independent analytic/Euler/native-selection tests and formal-entry1000-day representative acceptance, then compute the full campaign. There are320 conditions:170 culture conditions at three grids plus150 unchanged continuous conditions, giving660 workers. Every condition uses β=0/.5/1 and culture uses three routes per β. HPC uses xxlarge, array1–660 without a `%` concurrency limit. No refitting. Every culture condition/β/route must pass both137→139 and139→141 comparisons for daily shares, mean ploidy, log size, full-state CDF, boundary loss and calendar identity. A failure keeps results isolated and blocks promotion; no automatic further grid expansion. New data and figures are staged until all acceptance gates pass, then promoted by file allowlist.

Acceptance covers route-level complementarity and NA masks, median/reversed-IQR complementarity, empty day0, X/U/D and native passage-anchor agreement, no origin-specific rate, empty-group NA and legal CDFs. Report floating-point state/label correction counts and maximum magnitudes; do not claim zero numerical corrections or silently clip material negative states. Conditional means of numerically vanishing origin groups do not support mechanism claims. Figure8 and historical `_ver2` data, original receipts and images remain byte-identical. The ancestry diagnostics inside supp_fig7_8 follow the new current trajectories and are not a republication of `_ver2`.
