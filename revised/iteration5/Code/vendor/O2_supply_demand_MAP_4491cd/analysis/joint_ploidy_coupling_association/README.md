# Joint ploidy-category and coupling association

This folder classifies materialized 0–1000 day mean chromosome trajectories
and relates those classes to joint-fit parameter values and ratio classes. It
does not run the numerical model.

## Prespecified categories

- CatA: both 2N and 4N trajectories stay at or above `44 - 0.5` chromosomes.
  The tolerance is necessary because the 2N trajectory starts at exactly 44.
- CatB: both cohorts end at or below 30 chromosomes and each shows one
  effective collapse without a qualifying middle plateau.
- CatC: a cohort-level CatC curve contains at least two effective drops
  separated by a plateau of at least 60 days and favors the two-transition
  representation by `BIC_two - BIC_one <= -10`. At seed level, CatC includes
  CatC+CatC and the complementary CatB(2N)+CatC(4N) pattern: 2N starts at 44,
  which is already the middle plateau of the 4N 88→44→22 trajectory.
- CatU: missing, threshold-edge, unclassified, or biologically incompatible
  2N/4N results; the complementary CatB+CatC pattern is not treated as a
  disagreement.

Every numeric rule is exported in `ploidy_category_definition.tsv`; an 81-cell
sensitivity grid varies high tolerance, low endpoint, plateau duration, and
BIC support.

## Files

| File | Exact responsibility | Main outputs |
|---|---|---|
| `classify_joint_ploidy_trajectories.R` | Reads materialized `predict1000_ploidy_seed_day_mean.tsv`, calculates terminal levels, crossings, slopes, drop episodes, plateau diagnostics and one/two-transition BIC, then makes cohort and concordant seed categories. It also materializes pair assignments, pair-balanced category archetypes, deterministic representative trajectories, signed cutoff margins, and sensitivity-agreement summaries. | Cohort features, seed/pair categories, category/sensitivity/confidence/input-quality tables, plot-ready trajectory/archetype summaries and manifest. |
| `analyze_within_pair_ploidy_coupling.R` | Joins Cat labels to the soft-coupling master table and compares in-vivo values, in-vitro values, ratios, log2 ratios, deltas, penalties and objectives within each pair. | Cat×parameter summaries, Cat×Class proportions, Kruskal effect tables and joined master table. |
| `analyze_between_pair_ploidy_coupling.R` | Gives pairs equal weight, bootstraps pair means and performs leave-one-pair-out summaries of Cat×parameter effects. | Pair-balanced Cat×parameter summaries/effects and leave-one-pair-out table. |
| `analyze_ploidy_ratio_class_association.R` | Builds seed-level audit tables and pair-level descriptive Cat×Class summaries for every parameter, including pair-balanced class cells, natural-value/ratio summaries, dominant-class Cramér's V, normalized mutual information, pair-stratified tests only when estimable, and in-vivo/in-vitro/ratio representation comparisons. | Pair-level Cat×parameter/Class tables, global/cell association audits, explicit estimability table, metric-comparison table, report summary and manifest. |

Association outputs are model-derived and descriptive. They do not establish
causality or external biological validation. CatU is always retained in audit
tables and excluded only from the prespecified CatA/B/C association tests.
When a pair contains only one ploidy category, within-pair Cat comparisons and
pair-stratified Cat×Class tests are marked non-estimable rather than assigned a
misleading p-value; cross-pair Cramér's V remains a descriptive, pair-confounded
pattern.

Tests: `test-joint-ploidy-coupling-association.R` and
`test-joint-coupling-stage-split.R`.
