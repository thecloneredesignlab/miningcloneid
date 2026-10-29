# Additional Figure7 and supp_fig7_8 condition using environment-specific fitted p_misseg

This condition adds Figure7_fitted_p_misseg and supp_fig7_8_fitted_p_misseg. It does not replace fixed-p_misseg=0.1/0.3 products, Figure8, or historical_ver2. Publication status is established by this condition's completion.json and figure acceptance records, not by this document alone.

## Parameters and scenarios

Use the complete current winning joint fits: C01 seed94 and C02 seed80, with 26 shared_invivo parameters or 33 invitro_effective parameters. Culture configuration comes from fit_result.rds$ctx$invitro_cfg. Native p_misseg is 0.117778681235383 (C01 in vivo), 0.164105992427092 (C02 in vivo), and 0.0104720786599741 (both in vitro). It is not overridden; p_mis_base, k_o_mis and all other parameters remain context specific.

The effective probability remains p_eff(k,O2)=p_mis_base+p_misseg*mu_eff(k,O2)/(mu_eff(k,O2)+k_o_mis), rather than a forced constant effective probability. Use the dedicated SHA256-locked Figure7/8 model. Fixed oxygen (0.5% or 5%) and pure 2N/4N starts remain explicit diagnostic overrides; this does not restore fitted Gaussian starts or observed oxygen histories. One N corresponds to 22 chromosomes; low L is 22≤k<66, high H is k≥66. Beta=0,0.5,1 redirects the corresponding ordinary H-to-L daughter yield to the maternal state, preserving viable column yield. WGD retains the canonical multiplicity.

## Finite horizon and culture protocol

Each trajectory runs to 1000 cumulative days from its actual root, not an asymptotic equilibrium. In vivo uses continuous fixed-oxygen dynamics. Culture retains the original fitted DT, density feedback, inocula, targets, observation windows and positive-day nearest-target selector. Beyond the recorded prefix, control/O1/O2 separately cycle their last five complete measured passage tuples. Full-window selection precedes final truncation to day 1000. No inoculum-feasibility stop is added; reseeding is not birth, and no stochastic passage rule is introduced.

Grid convergence is established anew for these fitted values:137/139/141N in vitro and 7/9/11N in vivo. Adjacent grids must pass the existing limits for daily shares, overall mean ploidy, log count, state CDF, upper-bound loss and culture passage calendars. The retained grids are 141N and 7N respectively; failed conditions cannot be published.

## Daily complementary source accounting

At plotted day d, the interval is [d-1,d). C_HL and C_LL count viable low daughters produced by mothers currently high or low immediately before division. A=100*C_HL/(C_HL+C_LL), B=100*C_LL/(C_HL+C_LL). Count each viable low daughter once; do not subtract death or maternal removal. A former H-to-L daughter that divides while low contributes her new low daughters to B regardless of ancestry. Day0 and zero-production windows are NA in both panels.

## Current origin groups and distributions

D contains living direct H-to-L daughters that have not divided once. The label exits at first division or death, is not inherited, and has no arbitrary age window. The direct view compares X_L,D_L,X_L-D_L; the complement includes D descendants and must not be called never-flux. The permanent-lineage view compares X_L,T_L,U_L, where U has never undergone H-to-L and T=X-U; this label is inherited.

Cell shares divide by all current low cells. Low-daughter output shares divide by all current low mothers' expected viable low-daughter output. Mean division rate, mean viable low-daughter output and mean ploidy weight each group's current states, using current low-state lambda, viable daughter operator J and contemporaneous crowding. They do not inherit an ancestral high-state rate or subtract death. Daughter viability and WGD multiplicity remain canonical.

Culture metrics are calculated within each protocol route, then summarized by equal-route medians and type-7 IQRs. These describe three deterministic protocols, not 500 optimizer seeds, random repeats, biological uncertainty or confidence intervals. CDFs normalize cell weights within each nonempty route, then use equal-route mean CDFs and route IQRs, matching current supp_fig7_8.

## Conditional-mean validity and products

All three conditional means require identifiable origin, finite means, and group mass≥1e-8 of all current living cells. Theoretical empty groups are NA. All three culture routes must qualify before showing medians/IQRs; no interpolation or bridging across missing intervals is permitted. Expected count<1 is a reporting flag only. Such means are deterministic expected-mass ratios, not observed means of fractional living cells. Cell shares, output shares and CDFs do not receive this conditional-mean mask. Relative-threshold coverage sensitivity is reported at 1e-9,1e-8,1e-7.

The added main figure is the 0.5% O2 daily A/B page. The added supplement contains 12 time-profile pages and 12 unblocked CDF pages across both oxygens and origin views (24 total). Its PNG previews page 1 only. Retained data reside in data/Figures/Figure7/fitted_p_misseg; drawing needs its current prepared, tables, postprocess and provenance directories, without old versions or archives.
