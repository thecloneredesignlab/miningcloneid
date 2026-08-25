# Figure 6 multi-seed robustness revision audit

Status: completed and validated on 2026-08-06. All work and outputs in this audit are confined to `revised/iteration2`.

## Reviewer-request checklist

### 1. Justify the selected number of clusters and subclusters

- The primary warm-start partition was audited over `k = 2,...,8` on the frozen pooled 14-parameter t-SNE coordinates. The average silhouettes were 0.815694 for `k = 3` and 0.816651 for `k = 4`; the improvement was only 0.000957, and `k = 4` created a singleton. The retained design is therefore the smaller non-singleton `k = 3` partition.
- In 100 independent 80% endpoint subsamples drawn without replacement, re-selection favored `k = 3` 58 times and `k = 4` 42 times. This is reported as near-tied search-coverage structure, not evidence for three biological groups.
- A `k = 3` partition estimated directly in standardized 14-parameter endpoint space agreed poorly with the t-SNE partition (ARI 0.0523). This supports explicitly treating C01-C03 as embedding-dependent numerical-search regions.
- Within-region silhouette maximization selected `k = 2` for C01, C02, and C03, with silhouettes 0.179, 0.096, and 0.330. Because separation is weak to moderate, these divisions are called warm-start strata rather than subtypes.
- The frozen archive contains the 228,000 saved t-SNE coordinates but not their corresponding standardized 14-parameter input vectors. t-SNE seed/perplexity sensitivity cannot be re-estimated from this archive; that limitation is stated in the manuscript.

Primary evidence: `data/Figures/Figure6/cluster_selection_audit.tsv`, `cluster_k_selection.tsv`, `cluster_k_resample_selection_frequency.tsv`, and `cluster_stability.tsv`.

### 2. Define an objective/distribution-based method for excluding poor fits

- Hard eligibility requires a finite complete full joint MAP objective, valid unique seed, exactly one finite value for each of 14 in-vivo response parameters, feasibility before projection and at the reported endpoint, no applied projection, and an available local configuration. All 3,000 endpoints passed these checks.
- Within each warm-start pair, endpoints are ranked by the complete retained full joint MAP objective, including data-fit and soft-coupling terms.
- The primary set is the pair-specific lowest objective decile: 50 of 500 endpoints per pair. Nested lowest-5% and lowest-20% sets, containing 25 and 100 endpoints per pair, are used for threshold sensitivity.
- Inclusion is unweighted after selection. Objective values define eligibility but are not converted into likelihood, posterior, confidence, or biological weights.
- Every retained operator result must have status `ok`, finite dominant-mode outputs, and forced-input error no larger than `1e-8`. All 600 endpoints in the widest 20% set passed operator QC.

Primary evidence: `joint_seed_fit_qc.tsv`, `joint_seed_acceptance.tsv`, `joint_seed_acceptance_summary.tsv`, and `joint_multiseed_operator_qc.tsv`.

### 3. Evaluate multiple joint-fit seeds rather than only the best seed

- Figure 6D calculations use 50 accepted endpoints from each of six warm-start pairs, for 300 primary endpoints rather than six pair minima. The main panel displays C01Sc01, C02Sc01, and C03Sc02 as C01, C02, and C03 in a one-row layout; all six pairs remain in the analytical and robustness tables.
- Each endpoint is evaluated over 201 oxygen concentrations and 60 fixed effective per-chromosome missegregation probabilities, yielding 12,060 operator evaluations per endpoint and 3,618,000 primary evaluations.
- The unmodified fitted trajectory is also evaluated at all 201 oxygen values for every retained endpoint.

Primary evidence: `joint_multiseed_surface_summary.tsv` and `joint_multiseed_trajectory_summary.tsv`.

### 4. Average acceptable seeds and display their distribution

- Figure 6D displays the pointwise median dominant mean ploidy across the 50 primary endpoints in each warm-start pair. The legend reports ploidy in its original units while the color spacing is logarithmic.
- The unmodified trajectory is shown as a pointwise median with a 10th-90th percentile band.
- Figure 6E uses a denser direct-evaluation grid for C01Sc01, C02Sc01, and C03Sc02: each panel contains 496 curves at fixed effective per-chromosome missegregation probabilities from 0.005 through 0.500 in increments of 0.001; each curve contains 201 oxygen concentrations, and the vertical coordinate is the pointwise median dominant mean ploidy across 50 accepted numerical endpoints.
- Four black line-type overlays identify the already evaluated trajectories at fixed effective missegregation probabilities 0.01, 0.10, 0.20, and 0.30. The line types are solid, dash--dot--dot, dot--dash, and dotted, respectively. They add no new calculation or color encoding; the complete extracted rows are recorded in `figure6d_highlighted_trajectories.tsv`.
- Figure 6D marks regions where 50% of retained endpoints have spectral gap below 0.005 with clipped purple (`#9B59B6`) diagonal hatching and a dotted boundary. The hatching adds no fill and therefore preserves the underlying dominant-ploidy color scale; the overlay legend identifies the weak-gap region and the black model-implied median trajectory. The corresponding mapped boundary remains recorded in `figure6d_spectral_gap_boundary.tsv` for provenance but is no longer drawn in Figure 6E.
- Figure 6E does not interpolate the original 60-point Figure 6D grid. It performs direct operator calculations for all 496 fixed inputs. Exact duplicate parameter endpoints are evaluated once and restored with their original optimizer-seed multiplicity before calculating the pointwise median; the three displayed pairs contain 50, 15, and 50 unique parameter endpoints, representing 50 endpoints each.
- The two panels answer separate questions. Figure 6D reports final ploidy over the two-dimensional fixed-O2/fixed-missegregation grid and overlays the unmodified model-implied final effective missegregation trajectory. Figure 6E reports the O2-ploidy response conditional on each fixed missegregation input.
- White crosses mark grid locations where fewer than 80% of endpoints agree on the `dominant_mean_ploidy >= 2` versus `< 2` regime.
- Exact duplicate 14-parameter endpoints are additionally collapsed to one equal-weight endpoint in a sensitivity analysis. In the primary set, each pair contains 15-50 unique endpoints among its 50 optimizer seeds.

Primary evidence: `joint_unique_parameter_endpoint_counts.tsv`, `joint_seed_vs_unique_parameter_surface_comparison.tsv`, and `joint_seed_vs_unique_parameter_trajectory_comparison.tsv`.

### 5. Confirm whether oxygen-response conclusions remain consistent across seeds

- At fixed effective missegregation, every retained primary endpoint in every warm-start pair predicts higher dominant mean ploidy at 5% than at 0% O2 at low, intermediate, and high diagnostic missegregation. Pair-median changes span 0.466-1.678, 0.434-1.693, and 0.953-4.888 ploidy units, respectively.
- The unmodified coupled trajectory is not universal: it decreases in C01Sc01, C01Sc02, C02Sc01, and C03Sc01, but increases in C02Sc02 and C03Sc02. Pair medians span -1.172 to +1.175.
- All 48 pair-by-diagnostic modal results are invariant across the 5%, 10%, and 20% objective sets, with minimum within-set modal support of 100%.
- All 144 seed-weighted versus unique-endpoint-weighted qualitative comparisons retain the same modal result.
- The primary `q10` median high/low-ploidy classification is unchanged at every grid point after exact endpoint deduplication. C02Sc01 is the most duplicated primary pair, with 50 seeds representing 15 unique endpoints; its maximum pointwise median-ploidy difference after deduplication is 0.0149.
- A narrow local threshold sensitivity remains in the smaller C02Sc01 `q05` set: 25 seeds represent nine unique endpoints, grid-regime agreement is 0.999917, and the maximum local median-ploidy difference is 0.862. This does not change any prespecified modal conclusion, but it prevents claiming pointwise identity of all surfaces.
- The scientifically supported conclusion is conditional: oxygen raises the dominant ploidy at fixed effective CIN across the sampled accepted endpoints, whereas the emergent coupled oxygen-CIN-ploidy topology remains warm-start-pair dependent. The result is not presented as a causal intervention or biological confidence interval.

Primary evidence: `joint_seed_response_delta_summary.tsv`, `joint_seed_claim_robustness.tsv`, `joint_seed_cutoff_consistency.tsv`, and `joint_seed_robustness_audit_summary.tsv`.

### 6. Clarify the clustering-coordinate terminology

- The figure and manuscript use `Pooled 14-parameter t-SNE coordinate 1` and `Pooled 14-parameter t-SNE coordinate 2`.
- The coordinates are explicitly described as unitless visualization coordinates that approximate local neighborhoods; they are not parameters, latent biological variables, or mechanistic distances.
- The informal legacy label has been removed from the Figure 6 code and manuscript.
- The former combined panel B is now separated into Figure 6B, `Response classes in parameter space`, and Figure 6C, `Full MAP fit quality across classes`. The original response-surface and fixed-input curve-family panels are consequently labeled Figure 6D and Figure 6E. Every panel label includes a period and is immediately followed by its panel title.

## Delivered files

- Main figure: `Figures/assembled_fig6.png` and `Figures/assembled_fig6.pdf`
- Figure 6E panel (legacy internal `figure6d_*` file prefix retained for provenance): `data/Figures/Figure6/panels/pair_fixed_p_miss_eff_o2_ploidy_curve_family_three_pair_grid.png` and `.pdf`
- Figure 6E source tables (legacy internal `figure6d_*` file prefix retained): `data/Figures/Figure6/figure6d_fixed_p_curve_family.tsv`, `figure6d_highlighted_trajectories.tsv`, `figure6d_spectral_gap_boundary.tsv`, `figure6d_dense_endpoint_manifest.tsv`, `figure6d_dense_endpoint_qc.tsv`, and `figure6d_displayed_pairs.tsv`
- Response-class supplement: `Figures/supp_fig6-1_response_class_diagnostics.png` and `Figures/supp_fig6-1_response_class_diagnostics.pdf`
- Joint-ensemble robustness supplement: `Figures/supp_fig6-2_joint_ensemble_robustness.png` and `Figures/supp_fig6-2_joint_ensemble_robustness.pdf`
- Standalone data and drawing entry points: `Code/Figures/data_Figure6.R`, `draw_Figure6.R`, `data_Supp_Figure6_1.R`, `draw_Supp_Figure6_1.R`, `data_Supp_Figure6_2.R`, and `draw_Supp_Figure6_2.R`
- Analysis implementation: `Code/Figures/util/analysis/figure6_robustness.R`
- Manuscript integration: `manuscript/ltee_hypoxia_model.tex`
- Main validation: `data/Figures/Figure6/figure6_multiseed_validation.tsv` — 16/16 checks pass after dense-grid rendering
- Figure 6E model validation (legacy internal file prefix): `data/Figures/Figure6/figure6d_dense_model_validation.tsv`
- Figure 6E plotting validation (legacy internal file prefix): `data/Figures/Figure6/figure6d_validation.tsv` — 17/17 checks pass, including the four complete highlighted trajectories and retained mapped-boundary provenance
- Supplementary Figure 6-2 validation: `data/Figures/Supp_Figure6_2/supp_figure6-2_validation.tsv` — 14/14 checks pass
- Output manifests: `data/Figures/Figure6/figure6_output_manifest.tsv` and `data/Figures/Supp_Figure6_2/supp_figure6-2_output_manifest.tsv` — every listed artifact exists and has a recorded MD5 checksum
- Weak-gap regime-robustness figure: `Figures/supp_fig6-3_weak_gap_regime_robustness.png` and `Figures/supp_fig6-3_weak_gap_regime_robustness.pdf`
- Supplementary Figure 6-3 standalone entry points: `Code/Figures/data_Supp_Figure6_3.R` and `Code/Figures/draw_Supp_Figure6_3.R`
- Weak-gap methods and interpretation audit: `Supplementary_Figure6-3_weak_gap_regime_robustness.md`
- Supplementary Figure 6-3 validation: `data/Figures/Supp_Figure6_3/supp_figure6-3_data_validation.tsv` (27/27) and `data/Figures/Supp_Figure6_3/supp_figure6-3_figure_validation.tsv` (9/9)
- The earlier top-10 eigenspectrum outputs and `Archived_Figure6_eigenmode_audit.md` are retained as mathematical provenance but are not a manuscript-facing supplementary figure.

## Interpretation boundary

Numerical seeds are optimizer endpoints, not biological replicates or posterior draws. C01-C03 and Sc01-Sc02 are warm-start search strata, not biological subtypes. The fixed-O2 surfaces are post-fit operator diagnostics and do not substitute for oxygen-by-CIN perturbation experiments.
