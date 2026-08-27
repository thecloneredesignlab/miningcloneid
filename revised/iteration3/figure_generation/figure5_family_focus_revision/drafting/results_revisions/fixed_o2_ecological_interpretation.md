---
section_id: fixed_o2_ecological_interpretation
title: Sustained missegregation can reverse the oxygen-ploidy response at high oxygen
primary_claims: []
main_figures: [Figure 6]
supplemental_figures: [Supplementary Figures 6.2 and 6.3]
input_versions:
  figure_description: Figure 6 legend and fixed-input leading-eigenvector Methods in the iteration 3 manuscript
  claim_graph: not provided
---

## Results Text

We asked whether oxygen shifted long-run chromosome-state composition predictably when its fitted coupling to chromosome missegregation was removed. Within each of the six joint warm-start pairs, we retained the 50 endpoints with the lowest full joint objectives and independently imposed oxygen and effective per-chromosome missegregation across the fitted chromosome-state operator. Normalizing its leading eigenvector yielded the model-derived dominant mean ploidy at each condition; a small leading spectral gap marked conditions in which convergence toward that rank-1 composition would be slow. Figure 6 displays one prespecified primary pair from each landscape family.

The fixed-input response was consistent through the low-oxygen range but was not monotonic across the complete 0--5% O2 interval. In all three displayed ensembles, dominant ploidy rose as oxygen increased from anoxia through the hypoxic range (Figure 6B), directionally consistent with the reported positive association between tissue oxygen and cancer ploidy. At lower imposed missegregation probabilities, ploidy changed gradually or continued to rise. At higher probabilities, dominant ploidy instead reached a maximum and then declined as oxygen increased further, producing a high-O2 reversal in each displayed ensemble. All six warm-start pairs nevertheless retained positive 5%-minus-0% contrasts at the prespecified low, intermediate, and high missegregation inputs, showing why endpoint contrasts alone do not capture the intervening nonmonotonic response.

The chromosome-state operator provides a model-consistent route for this reversal. At a fixed per-chromosome missegregation probability, higher-chromosome mothers have more chromosome copies at risk and therefore a larger expected number of missegregated copies per division. Each modeled event generates complementary chromosome-gain and chromosome-loss daughters, and the post-missegregation survival filter retains a fraction of those altered daughters. Recurrent production and retention of chromosome-loss descendants can therefore replenish lower-chromosome states, whose smaller total chromosome-error burden per division and state-dependent fitness can increase their contribution to the dominant composition. This encoded high-to-low influx is consistent with the decline in dominant ploidy at high oxygen and sustained missegregation. Figure 6 reports the leading-eigenvector composition of the complete operator, however, rather than directional transition fluxes; it therefore does not quantify that influx or establish it as the unique cause of the reversal.

Restoring each endpoint's fitted oxygen-dependent missegregation function traced different paths through the fixed-input surfaces (Figure 6A). The median coupled trajectory decreased by 1.112 and 1.172 ploidy units in the displayed C01 and C02 pairs, respectively, but increased by 1.175 units in C03; across all six pairs, four trajectories decreased and two increased. Repeating the prespecified comparisons with the lowest 5%, 10%, and 20% of endpoints and with seed-weighted or unique-endpoint summaries did not change the modal results across the six pairs (Supplementary Figure 6.2). Small-gap regions remained locally sensitive in exact ploidy magnitude, but at least 90% of endpoints agreed on low, intermediate, or high dominant-ploidy class at every weak-gap cell in the three displayed ensembles (Supplementary Figure 6.3).

Thus, the fixed-input analysis predicts that oxygen favors higher ploidy through the hypoxic range, whereas sustained missegregation can shift the long-run composition back toward lower chromosome number at higher oxygen. Recurrent production and retention of chromosome-loss descendants is an encoded route to that result, not a separately quantified flux mechanism. Because the fitted oxygen-dependent missegregation functions traverse different parts of the response surface, the direction of the fully coupled oxygen--missegregation--ploidy response remains unresolved.

## Revision Notes

1. Figure 6 directly supports the nonmonotonic dominant-ploidy response of the complete fixed-input transition-growth operator, including the high-O2 decline at elevated imposed per-chromosome missegregation.
2. The operator explicitly generates complementary chromosome-gain and chromosome-loss daughters and applies post-missegregation viability filtering, so recurrent generation and retention of lower-chromosome descendants is a model-consistent route to the observed reversal.
3. Directional transition fluxes were not calculated, and no ablation separated chromosome-loss influx from high-ploidy reproductive loss, boundary effects, WGD, or other operator terms. The proposed route is therefore mechanistically supported by model structure but not demonstrated as the unique or quantitatively dominant cause.
4. Fixed O2 values are standardized model inputs in a post-fit asymptotic diagnostic, not measured oxygen interventions or finite-time fitted trajectories.
