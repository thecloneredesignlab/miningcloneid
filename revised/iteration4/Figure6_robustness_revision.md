# Figure 6 paired-context, multi-seed revision audit

Status: completed with paired *in vivo* and *in vitro* calculations on 2026-08-25. All writes in this revision are confined to `revised/iteration2`.

## Scope and provenance

- Figure 6 now evaluates both parameter vectors from the same retained joint endpoint. The *in vitro* joint rows never substitute the separate-fit seed-10 vector.
- All model calculations load the current-branch implementation under `oxygen/code/O2_supply_demand_MAP`. The exact files and MD5 checksums used are recorded in `data/Figures/Figure6/invitro_model_code_provenance.tsv`.
- The user reported synchronization with commit `83953a874401e42cd176432786f889a896adc959`; reproducibility is nevertheless defined by the recorded current source paths and MD5 values because the live files are not asserted to be byte-identical to that commit.
- Exact context-specific parameter duplicates are evaluated once and then restored to their optimizer-seed multiplicities before ensemble summaries.

## Reviewer-request checklist

### 1. Justify the selected number of clusters and subclusters

- Primary warm-start regions were audited over `k = 2,...,8` on the frozen pooled 14-parameter t-SNE coordinates. Average silhouettes were 0.815694 at `k = 3` and 0.816651 at `k = 4`; the difference was 0.000957, and `k = 4` created a singleton. The smaller non-singleton `k = 3` partition was retained.
- Across 100 independent 80% endpoint subsamples, `k = 3` was selected 58 times and `k = 4` 42 times. This is reported as near-tied numerical search structure, not three biological groups.
- Within-region silhouette maximization selected `k = 2` for C01, C02, and C03, with silhouettes 0.179, 0.096, and 0.330. These divisions are called warm-start strata rather than subtypes.
- The coordinate labels are now “Pooled 14-parameter t-SNE coordinate 1/2.” They are unitless visualization coordinates with no direct mechanistic interpretation; the legacy informal terminology is absent.

Evidence: `cluster_selection_audit.tsv`, `cluster_k_selection.tsv`, `cluster_k_resample_selection_frequency.tsv`, and `cluster_stability.tsv` under `data/Figures/Figure6`.

### 2. Define an objective-based method for excluding poor fits

- Hard eligibility requires a finite complete full joint MAP objective, valid unique seed, complete finite parameter records, feasibility before projection and at the reported endpoint, no applied projection, and an available configuration.
- Endpoints are ranked within each warm-start pair by the complete retained joint objective. The primary set is the lowest 10% (50/500); nested 5% and 20% sets contain 25 and 100 endpoints.
- Objective values determine eligibility only. They are not converted to likelihood, posterior, confidence, or biological weights.
- Every retained post-fit operator result must be finite, have status `ok`, and reproduce a forced effective missegregation input within `1e-8`.

Evidence: `joint_seed_fit_qc.tsv`, `joint_seed_acceptance.tsv`, `joint_seed_acceptance_summary.tsv`, `joint_multiseed_operator_qc.tsv`, and `joint_multiseed_operator_qc_invitro.tsv`.

### 3. Evaluate multiple joint-fit endpoints in both contexts

- Figure 6A evaluates 50 retained endpoints per pair, 201 oxygen values, and 60 fixed effective per-chromosome missegregation values in each context.
- All six warm-start pairs are analyzed; C01Sc01, C02Sc01, and C03Sc02 are displayed as C01, C02, and C03.
- This gives 7,236,000 complete-surface evaluations across the two contexts and six pairs.
- Figure 6B directly evaluates 496 fixed probabilities from 0.005 to 0.500 at 0.001 intervals and 201 oxygen values for each displayed pair. There are 11,465,040 distinct operator evaluations per context after exact-endpoint deduplication and 22,930,080 across the paired display.

Evidence: `joint_multiseed_surface_summary.tsv`, `joint_multiseed_surface_summary_invitro.tsv`, `figure6d_dense_model_validation.tsv`, and `figure6_invitro_dense_validation.tsv`.

### 4. Show ensemble summaries rather than only best seeds

- Figure 6A uses the seed-weighted pointwise median dominant mean ploidy. The unmodified fitted trajectory is shown as a pointwise median with a 10th--90th percentile band.
- Figure 6B shows the endpoint-level numerical inverse of the dense forward family. Four representative forward curves at fixed `p_miss,eff = 0.01, 0.10, 0.20, 0.30` and the equal-weight mean across all 496 fixed-input curves are overlaid.
- Supplementary Figure 6-2D reports the modal result and exact endpoint support for every context--pair--diagnostic combination, including cutoff and exact-endpoint-deduplication sensitivity.
- Supplementary Figure 6-3 reports regime agreement, local class sensitivity, across-endpoint ploidy spread, and local-jump distributions within weak-gap regions for both contexts.

### 5. Test whether conclusions remain consistent across seeds and contexts

- In the *in vivo* vectors, all six pair ensembles increase in dominant mean ploidy from 0% to 5% oxygen at low, intermediate, and high fixed missegregation inputs. Pair-median changes span `+0.466`--`+1.678`, `+0.434`--`+1.693`, and `+0.953`--`+4.888`.
- In the paired *in vitro* vectors, the corresponding response is approximately flat at the low input (`+0.00103`--`+0.00168`), increases modestly at the intermediate input (`+0.0807` after rounding), and decreases at the high input (`-0.823` after rounding).
- Restoring each fitted missegregation function gives four decreasing and two increasing *in vivo* trajectories, but approximately flat *in vitro* trajectories in every pair.
- All 96 context--pair--diagnostic modal results are unchanged across the 5%, 10%, and 20% objective sets, with 100% within-set support. All 288 seed-weighted versus unique-context-endpoint comparisons retain the same modal result.
- The inverse map supports one stable value in 27.67%--69.84% of the displayed *in vivo* domains and 84.28% of every displayed *in vitro* domain. At target ploidy 4, supported *in vivo* values rise between 1% and 5% oxygen, whereas *in vitro* values fall from 0.02174 at 0% to 0.006725 at 5%.
- At least 99.93% of weak-gap cells in every displayed context--pair panel retain 90% endpoint agreement on low/intermediate/high dominant ploidy. Narrow boundary cells can still show large local changes, so this is regime robustness rather than exact pointwise certainty.

Evidence: `joint_seed_claim_robustness*.tsv`, `joint_seed_cutoff_consistency*.tsv`, `joint_seed_vs_unique_parameter_robustness*.tsv`, `figure6*_inverse_*summary.tsv`, and `data/Figures/Supp_Figure6_3/supp_figure6-3_weak_gap_pair_summary.tsv`.

## Delivered figure structure

- Main Figure 6A: paired *in vivo*/*in vitro* oxygen--CIN--ploidy response surfaces, three displayed pair columns.
- Main Figure 6B: paired context-specific inverse maps with fixed-input reference curves and the equal-weight dense-grid mean.
- Supplementary Figure 6-1A/B: separate-fit response classes for 500 *in vivo* and 500 *in vitro* endpoints in aligned 8-by-1 class layouts; absent classes are blank.
- Supplementary Figure 6-1C: pooled parameter-space overlay with circles for *in vivo*, triangles for *in vitro*, and outlined context--class objective winners.
- Supplementary Figure 6-1D: split context-specific relative-objective distributions by response class.
- Supplementary Figure 6-2A--C: shared cluster and objective-selection diagnostics; panel D: context-specific multi-seed/cutoff robustness.
- Supplementary Figure 6-3A--D: paired-context weak-gap regime robustness.

## Validation

- Main context layout: `figure6_context_validation.tsv`, 6/6 checks passed.
- Main multi-seed workflow: `figure6_multiseed_validation.tsv`, 22/22 checks passed.
- *In vitro* dense grid: `figure6_invitro_dense_validation.tsv`, 4/4 checks passed.
- *In vitro* inverse: `figure6_invitro_inverse_validation.tsv`, 14/14 checks passed.
- Supplementary Figure 6-1: `supp_fig6-1_validation.tsv`, 7/7 checks passed.
- Supplementary Figure 6-2: `data/Figures/Supp_Figure6_2/figure_validation.tsv`, 4/4 checks passed; analytical validation 14/14.
- Supplementary Figure 6-3: context validation 7/7, data validation 27/27, and figure validation 11/11.

## Interpretation boundary

Optimizer endpoints are neither biological replicates nor posterior draws. C01--C03 and Sc01--Sc02 are numerical search strata, not biological subtypes. Fixed-oxygen surfaces and inverse maps are post-fit asymptotic operator diagnostics; they do not establish a causal oxygen-by-CIN effect or a biological requirement for one missegregation probability.
