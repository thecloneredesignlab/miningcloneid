# Figure 5D DE-initial/optimizer revision

Status: active as panel D of the iteration2 manuscript Figure 5. This record
retains its legacy filename for build-history compatibility. The earlier
generalized-posterior workflow is not an active input and no HPC output is
required.

## Scientific scope

Figure 5D compares two numerical-search distributions for the *in vivo* and
*in vitro* values of each of the 14 paired parameters:

1. the pair-specific DE initial population reconstructed by replaying the
   actual `joint_deoptim_initial_population()` routine for 500 joint seeds and
   400 paired population members per seed; and
2. 500 feasible, unprojected optimizer endpoints from one selected pair in
   each of C01, C02, and C03.

Within each family, the selected pair is the one whose in-vivo member has the
lower original standalone separate-in-vivo objective. This gives C01 seed366,
C02 seed25, and C03 seed311, each paired with the common in-vitro seed10
anchor.

## Changes from the preceding active version

- Removed every generalized-posterior distribution, credible interval,
  convergence diagnostic, temperature-sensitivity result, and family-basin
  statement from active Figure 5D, its caption, methods, results, and SI table.
- Corrected the meaning of the reference distribution: it is the actual DE
  initial population, not a distribution induced by fitting intervals or
  weighted by the Welsch objective penalty.
- Replayed seed1--seed500 with `NP_used=400` and
  `joint_warmup_sigmaN=0.1216`, producing 200,000 initial values per family and
  parameter.
- Made the gray distribution pair-specific in C01/C02/C03 because the
  warm-start centers, deltas, and feasible delta bounds differ among pairs.
- Reconstructed the two context values from each sampled center and delta and
  display them separately on the parameter's natural scale. The upper blue
  half is *in vivo* and the lower pink half is *in vitro*; neither context is
  divided by the other in the active panel.
- Used a gray dashed outline for the actual pair-specific DE initial
  population and a solid C01/C02/C03 outline for the optimizer endpoints.
- Used separate peak normalization for each density and a shared x grid within
  each parameter; the graph therefore compares horizontal location, width, and
  shape rather than peak height.
- Added medians and 5th--95th percentile spans for both distributions.
- Consolidated the repeated C01/C02/C03 strips into one shared top header,
  placed the upper-*in vivo*/lower-*in vitro* orientation labels at the left of
  every parameter row, and moved each parameter name and compact functional
  description to a two-line label at the right.
- Added a Supplementary Table with the exact source paths and calculation rules
  in TeX comments, per-family summaries, endpoint/initial width ratios, active
  bound fractions, and cross-family directions.
- Replaced identifiability claims with the narrower statement that endpoint
  concentration measures numerical fitted-solution stability under the tested
  search design.

## Evidence boundary

The initial populations and optimizer endpoints are numerical-search
distributions. Neither is a Bayesian posterior, confidence distribution, or
biological replicate distribution. Endpoint concentration relative to
initialization shows repeated numerical convergence under the tested design
but does not prove structural identifiability. An exactly zero 5th--95th
percentile endpoint span can mean that at least 90% of starts reached the same
reported parameter value. Active-bound occupancy is reported separately by
context because convergence at a bound may reflect optimizer pressure rather
than a resolved biological optimum. The paired ratio is retained only as a
secondary directional/search audit.

## Active products

- `Code/Figures/prepare_Figure5F_de_initial_population.R`
- `Code/Figures/audit_Figure5F_prior_optimizer_inputs.R`
- `Code/Figures/build_Figure5F_prior_optimizer_products.R`
- `Code/Figures/build_Figure5F_supplementary_table.R`
- `Code/Figures/draw_Figure5.R`
- `data/Figures/Figure5/figure5f_prior_optimizer_density.tsv`
- `data/Figures/Figure5/figure5f_prior_optimizer_summary.tsv`
- `data/Figures/Figure5/figure5f_prior_optimizer_cross_family.tsv`
- `data/Figures/Figure5/figure5f_prior_optimizer_readiness.tsv`
- `data/Figures/Figure5/figure5f_de_initial_population_context_values.rds`
- `data/Figures/Figure5/figure5f_context_initial_optimizer_density.tsv`
- `data/Figures/Figure5/figure5f_context_initial_optimizer_summary.tsv`
- `manuscript/tables/supp_figure5f_prior_optimizer.tex`
- `manuscript/tables/data/supp_figure5f_prior_optimizer_values.tsv`
- `data/Figures/Figure5/panels/figure5d_context_initial_optimizer_endpoints.png`
- `Figures/assembled_fig5.png`

## Verified numerical summary

- 67,368 direct-context density rows
  (`14 x 3 families x 2 contexts x 2 roles x 401`).
- 84 direct-context family-by-parameter summaries
  (`14 x 3 families x 2 contexts`).
- Cross-family endpoint direction agreement for 12 of 14 parameters.
- Direction exceptions: `p_wgd` and `n_O`.
- The width-ratio and direction checks are paired-ratio audits rather than the
  coordinate of the active panel. All 14 paired ratios had endpoint 90% widths
  below the corresponding DE-initial width in all three selected families; 13
  were below half the initial width in all three families. `p_wgd` reached
  0.690 of the initial width in one family, and 11 parameters had a maximum
  family-specific width ratio of zero.
- Six parameters had at least 10% active-bound occupancy in one or more
  selected families: `alpha_o2`, `gamma_growth`, `mu_hp`, `gamma_mu`,
  `buffer_smax`, and `buffer_n_exp`.
- `data/Figures/Figure5/figure5f_de_initial_population_log2_ratios.rds`
- `data/Figures/Figure5/figure5f_de_initial_population_config.tsv`
- `data/Figures/Figure5/figure5f_de_initial_population_readiness.tsv`
