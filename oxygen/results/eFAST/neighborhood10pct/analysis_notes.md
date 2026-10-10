# Interpretation and figure guide: all 500 in-vivo neighborhoods

## Mechanistic interpretation

This analysis provides limited support for the complete proposed transition
from death dominance at low oxygen to greater missegregation and buffering
influence at high oxygen. The declared comparison uses normalized oxygen-window
AUCs over 0–1% and 3–5% oxygen, averaged across parameters within each mechanism.
Death-specific parameters are `mu_hp` and `gamma_mu`; the shared stress-response
parameters `O2_crit` and `n_O` form a separate group.

All four comparisons have 500 evaluable endpoints. The following percentages
come from `summaries/mechanism_support.tsv`:

| Criterion | Ploidy S1 | Ploidy ST | Growth S1 | Growth ST |
| --- | ---: | ---: | ---: | ---: |
| Death has the largest mechanism mean at low O2 | 12.8% | 11.2% | 22.8% | 21.0% |
| Missegregation has a higher mechanism mean at high O2 | 54.0% | 50.8% | 63.8% | 66.4% |
| Buffering has a higher mechanism mean at high O2 | 58.2% | 55.4% | 36.0% | 40.2% |
| All three conditions hold within the same endpoint | 7.6% | 6.6% | 5.2% | 5.6% |

The median ploidy S1 mechanism means are consistent with death influence
decreasing (0.0364 to 0.0041) and missegregation increasing (0.0239 to 0.0441).
Buffering is already influential at low oxygen (0.1187), with a modest increase
to 0.1310 at high oxygen. For growth, buffering decreases from 0.1136 to 0.0291,
while missegregation increases from 0.0437 to 0.0709. Thus, individual parts of
the proposed shift occur, but low-oxygen death dominance and the complete
three-part shift have weak support across the fitted endpoints. These
mechanism means summarize parameter indices; ST means overlap through shared
interactions and are not additive group variance fractions.

These conclusions are conditional on independent parameter sampling within
each fitted optimum ±10% of the full natural parameter span, clipped to the
original in-vivo bounds. Five phase estimates are averaged within each
endpoint, and all endpoints receive equal weight. The separate global-range
analysis is in the parent folder and uses a different sampling domain.

## Reading `figures/efast_four_panel.pdf`

- **A: S1, white to deep purple.** Each value estimates the parameter's main
  contribution to output variance within the neighborhood design.
- **B: ST, white to orange.** Each value includes the parameter's main effect
  and interactions involving it. ST−S1 does not identify a particular interacting
  parameter pair. Oxygen is fixed during each sensitivity analysis.
- **Left: dominant mean ploidy; right: asymptotic net live-cell growth rate.**
  The growth output is the fixed-oxygen operator's leading real eigenvalue.
  The growth output has units day⁻¹; the displayed sensitivity indices are dimensionless.
  Both heatmaps show the median across 500 endpoints of their five-phase means.
  The source table also gives means, quartiles and valid-endpoint counts.
- **Four independent colorbars.** Each starts at zero and ends at the maximum
  of its own displayed data. Compare numerical colorbar values when comparing
  outputs or S1 with ST. Darker color indicates greater sensitivity, without
  indicating whether increasing a parameter increases or decreases the output.
  The retained correlation panel is `../figures/figure4b_spearman_direction.pdf`.
- **Three row annotations:** 1 is the Figure 4 process annotation; 2 preserves
  the iteration5 Figure 4B O2 correlation group; 3 is the new sensitivity group.
  Each panel's row order follows its own ploidy S1/ST sensitivity group and then
  descending global peak. Growth follows that order to keep the two outputs aligned.
- **Oxygen axis:** 201 measured points from 0 to 5%, using a symmetric-log
  display axis with 0.025% linear threshold. Dashed lines mark 1% and 3%, the
  boundaries used for the low- and high-oxygen windows. AUCs use physical oxygen
  coordinates, rather than distances on the transformed display axis.

The new group follows the Figure 4B decision rule: a median-curve global peak
strictly above 0.3 and a low-versus-high contrast with BH-adjusted q<0.05.
Complete endpoint curves are resampled after retaining their five-phase means.
All ploidy S1 peaks are below 0.3, so every S1 row falls in the operational
`O2-independent` category. This label means that the curve did not meet this
classification rule; parameters can still have substantial S1 contributions
and oxygen-dependent curves. For ST, `buffer_smax` and `p_mis_base` are High O2;
`buffer_beta`, `mu_hp` and `O2_crit` are Low O2. Interpret these ST groups with the
repeat diagnostics below. Group assignments and exact statistics are in
`efast_o2_sensitivity_classification.tsv`.

## Convergence and numerical validity

All 500 endpoints were retained. Of the endpoint-level combined checks, two
passed, 87 have stable repeats but untested resolution, and 411 have at least
one failed check. The five-phase p90 range threshold is 0.10; the N257-to-N513
p90 change threshold is 0.05. These are practical index-unit tolerances.

| Output/index | Repeat checks passed | Resolution checks passed / tested |
| --- | ---: | ---: |
| Ploidy S1 | 500/500 | 8/8 |
| Ploidy ST | 90/500 | 2/8 |
| Growth S1 | 495/500 | 7/8 |
| Growth ST | 474/500 | 6/8 |

Ploidy ST is the principal numerical convergence limitation. Its interaction
and grouping interpretation remains provisional. Repeat stability of S1 does
not establish resolution convergence for the 492 untested endpoints. The
diagnostic tables provide each endpoint's separate statuses.

The 590 recovered model evaluations used their unchanged canonical matrices,
with nonnegative normalized eigenvectors and agreement at 50 and 100 decimal
digits. The maximum original-versus-recovered change was 5.30×10⁻⁶ in ploidy
and 6.49×10⁻¹¹ day⁻¹ in growth. The maximum S1/ST change was 1.47×10⁻⁸, and
no pattern of defined/undefined indices changed. Matrix proof artifacts and
the valid original rows remain on HPC; compact corrections and index-impact
tables are under `slurm/`. Numerical recovery addresses invalid evaluations;
the phase and resolution convergence diagnostics remain separate.
