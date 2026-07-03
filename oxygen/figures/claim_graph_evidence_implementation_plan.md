# Claim Graph Evidence Implementation Plan

## Purpose

This plan translates the panel-first claim graph into figure-content work. The claims are treated as fixed for this version. The purpose here is not to decide whether each claim is sufficiently supported, but to identify the existing reports, figures, tables, and code paths that can already supply figure content for each claim node.

No new analysis code is proposed in this pass. The first implementation should reuse, extract, annotate, or recompose existing outputs and reports wherever possible.

## Working Rule

For each claim node, record:

- the intended figure/panel;
- the existing evidence assets already available;
- the existing code or report workflow that generated or can regenerate those assets;
- the first no-code action needed to turn the asset into a manuscript panel;
- deferred evidence questions for later audit.

The first panel pass should be a coherent visual storyboard, not a final evidence audit.

## Existing Evidence Assets

### Data Stream Figures

Primary existing assets:

- `docs/optimization_data_streams/optimization_data_streams_overview.png`
- `docs/optimization_data_streams/optimization_data_streams_overview.pdf`
- `docs/optimization_data_streams/in_vivo_optimization_data_streams.png`
- `docs/optimization_data_streams/in_vivo_optimization_data_streams.pdf`
- `docs/optimization_data_streams/in_vitro_optimization_data_streams.png`
- `docs/optimization_data_streams/in_vitro_optimization_data_streams.pdf`
- `figures/in_vivo_optimization_data_streams.pdf`
- `figures/in_vitro_optimization_data_streams.pdf`
- `figures/optimization_data_streams_overview.pdf`

Supporting tables:

- `docs/optimization_data_streams/measurement_events.tsv`
- `docs/optimization_data_streams/invivo_burden_long.tsv`
- `docs/optimization_data_streams/invivo_ploidy_cells.tsv`
- `docs/optimization_data_streams/invivo_necrosis_endpoint.tsv`
- `docs/optimization_data_streams/invitro_lineage_timeline.tsv`
- `docs/optimization_data_streams/invitro_passage_observations.tsv`
- `docs/optimization_data_streams/invitro_kary_cells.tsv`
- `docs/optimization_data_streams/invitro_flow_density_grid.tsv`
- `docs/optimization_data_streams/optimization_data_streams_manifest.tsv`
- `docs/optimization_data_streams/optimization_data_streams_qc.tsv`

Context reports:

- `docs/in_vivo_optimization_data_streams_report.md`
- `docs/in_vitro_optimization_data_streams_report.md`

### Fit Reports

Packaged report archive:

- `/Users/4470246/Downloads/reports.zip`

Relevant contents:

- `in_vitro/fit_report_seed10.html`
- `in_vivo_top_10/rank01_fit_report_seed146.html`
- `in_vivo_top_10/rank02_fit_report_seed22.html`
- `in_vivo_top_10/rank03_fit_report_seed454.html`
- `in_vivo_top_10/rank04_fit_report_seed195.html`
- `in_vivo_top_10/rank05_fit_report_seed336.html`
- `in_vivo_top_10/rank06_fit_report_seed370.html`
- `in_vivo_top_10/rank07_fit_report_seed447.html`
- `in_vivo_top_10/rank08_fit_report_seed292.html`
- `in_vivo_top_10/rank09_fit_report_seed51.html`
- `in_vivo_top_10/rank10_fit_report_seed438.html`
- `multi-warm-up_results.html`
- `multi_warmup_joint/v01_i01_fit_report_seed440.html`
- `multi_warmup_joint/v02_i01_fit_report_seed392.html`
- `multi_warmup_joint/v03_i01_fit_report_seed330.html`
- `multi_warmup_joint/v04_i01_fit_report_seed129.html`
- `multi_warmup_joint/v08_i01_fit_report_seed118.html`

Existing report/render code:

- `oxygen/code/O2_supply_demand_MAP/report/render_fit_report.R`
- `oxygen/code/O2_supply_demand_MAP/report/render_invitro_fit_report.R`
- `oxygen/code/O2_supply_demand_MAP/util/o2_supply_demand_map_joint_parameter_plot.R`
- `oxygen/code/O2_supply_demand_MAP/analysis/multi_warmup/multi_warmup_results_report.R`
- `oxygen/code/O2_supply_demand_MAP/analysis/fit_results/extra_results_report.R`

### Joint Soft-Coupling Outputs

Existing manuscript-facing assets:

- `docs/fit_report_seed1_ratio_vivo_to_vitro_all_soft.png`
- `docs/fit_report_seed1_ratio_vivo_to_vitro_all_soft.pdf`
- `docs/fit_report_seed1_ratio_vivo_to_vitro_all_soft.tsv`
- `docs/fit_joint_sigma1e6_fit_report_seed38_ratio_vivo_to_vitro.png`
- `docs/fit_joint_sigma1e6_fit_report_seed38_ratio_vivo_to_vitro.pdf`
- `docs/fit_joint_sigma1e6_fit_report_seed38_ratio_vivo_to_vitro.tsv`

Existing code:

- `oxygen/code/O2_supply_demand_MAP/util/o2_supply_demand_map_joint_parameter_plot.R`
- `oxygen/code/O2_supply_demand_MAP/report/render_fit_report.R`
- `oxygen/code/O2_supply_demand_MAP/analysis/fit_results/compare_joint_sigma_soft_coupled_paired_seeds.R`

### Fixed-O2 And Resource-Regime Reports

Existing report assets:

- `/Users/4470246/Downloads/mode_parameter_contribution_report.html`
- `/Users/4470246/Downloads/mode_parameter_contribution_report (1).html`
- `/Users/4470246/Downloads/fixo2_eigen_attractor_embedding_curve_class_report.html`
- `/Users/4470246/Downloads/pooled_embedding_curve_class_report.html`

Existing code:

- `oxygen/code/O2_supply_demand_MAP/analysis/parameter_landscape_clustering/mode_parameter_contribution_analysis.R`
- `oxygen/code/O2_supply_demand_MAP/analysis/parameter_landscape_clustering/mode_parameter_contribution_report.R`
- `oxygen/code/O2_supply_demand_MAP/analysis/parameter_landscape_clustering/dominant_ploidy_parameter_contribution_analysis.R`
- `oxygen/code/O2_supply_demand_MAP/analysis/parameter_landscape_clustering/dominant_ploidy_parameter_contribution_report.R`
- `oxygen/code/O2_supply_demand_MAP/analysis/best_fit_parameter_feature/01_fixed_o2/FixO2_invivo.R`
- `oxygen/code/O2_supply_demand_MAP/analysis/best_fit_parameter_feature/05_FixO2_eigen_attractor_based_clustering/fixo2_eigen_attractor_feature_builder.R`
- `oxygen/code/O2_supply_demand_MAP/analysis/best_fit_parameter_feature/05_FixO2_eigen_attractor_based_clustering/fixo2_eigen_attractor_report.R`
- `oxygen/code/O2_supply_demand_MAP/analysis/best_fit_parameter_feature/06_combine_FixO2_eigen_attractor/plot_fixo2_eigen_attractor_embedding_curve_class.R`
- `oxygen/code/O2_supply_demand_MAP/analysis/best_fit_parameter_feature/06_combine_FixO2_eigen_attractor/render_fixo2_eigen_attractor_embedding_curve_class_report.R`

Planning docs:

- `docs/fixed_o2_ploidy_monotonicity_implementation_plan.md`

### Parameter Landscape And Candidate Regime Outputs

Primary archive:

- `/Users/4470246/Downloads/landscape_subcluster.zip`

Key contents:

- `landscape_subcluster/seed_cluster_original_fit_prior_positions.csv`
- `landscape_subcluster/seed_cluster_original_fit_prior_positions_README.md`
- `landscape_subcluster/SeedParameterTables/invivo_best_params_by_seed.csv`
- `landscape_subcluster/SeedParameterTables/invitro_best_params_by_seed.csv`
- `landscape_subcluster/pooled_invivo_invitro/full_data_in_vivo_clustring/Figures/pooled_invivo_invitro_initial_vs_best_umap_best_clusters.png`
- `landscape_subcluster/pooled_invivo_invitro/full_data_in_vivo_clustring/Figures/pooled_invivo_invitro_initial_vs_best_tsne_best_clusters.png`
- `landscape_subcluster/pooled_invivo_invitro/full_data_in_vivo_clustring/Tables/pooled_invivo_invitro_best_seed_groups_by_method.csv`
- `landscape_subcluster/pooled_invivo_invitro/full_data_in_vivo_clustring/Tables/pooled_invivo_invitro_best_subclusters_by_method.csv`
- `landscape_subcluster/pooled_invivo_invitro/full_data_in_vivo_clustring/Tables/pooled_invivo_invitro_best_subcluster_summary_by_method.csv`

Additional local derived plots:

- `/Users/4470246/Downloads/landscape_subcluster_tsne_invivo_driver_violins.png`
- `/Users/4470246/Downloads/landscape_subcluster_tsne_invivo_driver_violins.pdf`
- `/Users/4470246/Downloads/landscape_subcluster_tsne_invivo_pmisbase_nO_violins.png`
- `/Users/4470246/Downloads/landscape_subcluster_tsne_invivo_pmisbase_nO_violins.pdf`

Existing reports:

- `/Users/4470246/Downloads/parameter_landscape_clustering_umap_cluster_report (2).html`
- `/Users/4470246/Downloads/parameter_landscape_clustering_pca_cluster_report (2).html`
- `/Users/4470246/Downloads/parameter_landscape_pooled_pca_interpretation.html`
- `/Users/4470246/Downloads/pooled_embedding_curve_class_report.html`

Existing code:

- `oxygen/code/O2_supply_demand_MAP/util/build_multi_warmup_pairs_from_landscape_subclusters.R`
- `oxygen/code/O2_supply_demand_MAP/analysis/parameter_landscape_clustering/clustering_runner.R`
- `oxygen/code/O2_supply_demand_MAP/analysis/parameter_landscape_clustering/clustering_analysis.R`
- `oxygen/code/O2_supply_demand_MAP/analysis/parameter_landscape_clustering/clustering_report.R`
- `oxygen/code/O2_supply_demand_MAP/analysis/parameter_landscape_clustering/full_data_in_vivo_clustring.R`
- `oxygen/code/O2_supply_demand_MAP/analysis/parameter_landscape_clustering/parameter_landscape_utils.R`
- `oxygen/code/O2_supply_demand_MAP/analysis/best_fit_parameter_feature/04_combine_parameter_landscape/plot_pooled_embedding_curve_class.R`
- `oxygen/code/O2_supply_demand_MAP/analysis/best_fit_parameter_feature/04_combine_parameter_landscape/render_pooled_embedding_curve_class_report.R`

Planning docs:

- `docs/parameter_landscape_prior_aware_preprocessing_plan.md`

## Figure-by-Figure Evidence Plan

### Figure 1: Experimental System And Data Streams

Claim nodes:

- C1: Matched SUM159 near-2N and near-4N lineages reveal environment-dependent ploidy selection across in vitro and in vivo settings.
- C3: The tumor setting is not simply lower oxygen; it reflects a broader in vivo stress regime.

Existing evidence to reuse:

- `docs/optimization_data_streams/optimization_data_streams_overview.pdf`
- `docs/optimization_data_streams/in_vivo_optimization_data_streams.pdf`
- `docs/optimization_data_streams/in_vitro_optimization_data_streams.pdf`
- `docs/optimization_data_streams/measurement_events.tsv`

First no-code panel actions:

1. Use the existing overview figure as the first draft of Fig. 1.
2. Split or crop the in vivo and in vitro stream figures into subpanels if the overview is too dense.
3. Add panel captions emphasizing:
   - same isogenic 2N/4N system;
   - absolute time-frame difference;
   - terminal in vivo ploidy/necrosis versus longitudinal in vitro oxygen-known passage history;
   - oxygen is known in vitro and latent/effective in vivo.

Deferred evidence questions:

- Whether Figure 1 should include all data streams or reserve some streams for supplement.
- Whether the current `figures/` copies or `docs/optimization_data_streams/` copies should be the canonical source.

### Figure 2: Model Mechanism

Claim nodes:

- R0: Resource limitation rewires ploidy evolution through opposing costs and buffering advantages.
- C2: Resource limitation can both oppose and promote high ploidy.
- C8: In vitro ploidy reduction can arise from CIN-generated variation plus missegregation survival buffering.

Existing evidence/code to reuse:

- Model implementation:
  - `oxygen/code/O2_supply_demand_MAP/model/model_O2_supply_demand_MAP.R`
  - `oxygen/code/O2_supply_demand_MAP/model/model_O2_supply_demand_MAP.cpp`
  - `oxygen/code/O2_supply_demand_MAP/util/o2_supply_demand_map_common_semantics.R`
- Existing fit-report schematic/report sections:
  - `oxygen/code/O2_supply_demand_MAP/report/render_fit_report.R`
  - `oxygen/code/O2_supply_demand_MAP/report/render_invitro_fit_report.R`
- Current claim explanation:
  - `/Users/4470246/Downloads/notes_hypoxiaLTEE.txt`

First no-code panel actions:

1. Draft a schematic panel from the model components already implemented:
   - O2/resource stress;
   - proliferation damping;
   - stress-associated death;
   - chromosome missegregation;
   - WGD;
   - post-missegregation survival;
   - high-ploidy buffering;
   - viable chromosome-loss daughters.
2. Use model terms and parameter names from fit reports/code to avoid inventing mechanisms.
3. Place C8 as a subpanel or inset showing how high ploidy can both persist and move downward through viable chromosome-loss products.

Deferred evidence questions:

- Whether the schematic should be pure conceptual art or include parameter/function annotations.
- Whether direct model equations belong in main Figure 2 or supplement.

### Figure 3: Separate In Vitro And In Vivo Fits

Claim nodes:

- C7: Separate in vitro fits support a culture regime where oxygen stress reshapes ploidy without strong elimination of high-ploidy cells.
- C1/C3 support: separate fits connect the same data-stream setup to distinct context behavior.

Existing evidence to reuse:

- `/Users/4470246/Downloads/reports.zip`
  - `in_vitro/fit_report_seed10.html`
  - `in_vivo_top_10/*.html`
- Current extracted/available reports in `/Users/4470246/Downloads/reports/` if present:
  - `in_vitro/fit_report_seed10.html`
  - `in_vivo_top_10/rank01_fit_report_seed146.html`
- Existing report code:
  - `oxygen/code/O2_supply_demand_MAP/report/render_fit_report.R`
  - `oxygen/code/O2_supply_demand_MAP/report/render_invitro_fit_report.R`

First no-code panel actions:

1. Use the in vitro seed10 report as the current in vitro milestone panel source.
2. Extract or screenshot the existing growth/ploidy fit panels that show:
   - low oxygen slows proliferation;
   - ploidy trajectories are reshaped;
   - high-ploidy states remain present/permissive.
3. Extract parameter values or report table snippets showing weak lower-bound high-ploidy growth-penalty terms:
   - `alpha_o2 = 0.500`
   - `gamma_growth = 1.500`
4. Use the in vivo top-10 reports as context, but do not overload Fig. 3 with all top seeds. Pick one representative in vivo fit or show a compact top-fit montage.

Deferred evidence questions:

- Whether to show one best in vivo seed, multiple representative in vivo seeds, or only in vitro for the main Fig. 3.
- Whether the in vitro-only fit should be re-exported as standalone PNG/PDF panels rather than screenshots from HTML.

### Figure 4: Joint In Vivo/In Vitro Context Differences

Claim nodes:

- C4: Joint soft-coupled fits suggest that in vivo differs from in vitro in how stress is translated into proliferation, chromosome instability, and survival filtering.

Existing evidence to reuse:

- `/Users/4470246/Downloads/reports.zip`
  - `multi-warm-up_results.html`
  - `multi_warmup_joint/*.html`
- Existing ratio figures:
  - `docs/fit_report_seed1_ratio_vivo_to_vitro_all_soft.png`
  - `docs/fit_report_seed1_ratio_vivo_to_vitro_all_soft.pdf`
  - `docs/fit_report_seed1_ratio_vivo_to_vitro_all_soft.tsv`
  - `docs/fit_joint_sigma1e6_fit_report_seed38_ratio_vivo_to_vitro.png`
  - `docs/fit_joint_sigma1e6_fit_report_seed38_ratio_vivo_to_vitro.pdf`
  - `docs/fit_joint_sigma1e6_fit_report_seed38_ratio_vivo_to_vitro.tsv`
- Existing code:
  - `oxygen/code/O2_supply_demand_MAP/util/o2_supply_demand_map_joint_parameter_plot.R`
  - `oxygen/code/O2_supply_demand_MAP/analysis/multi_warmup/multi_warmup_results_report.R`
  - `oxygen/code/O2_supply_demand_MAP/analysis/fit_results/compare_joint_sigma_soft_coupled_paired_seeds.R`

First no-code panel actions:

1. Use existing in vivo/in vitro ratio forest plots as the first Fig. 4 draft.
2. Add a compact table or annotation for the three context-specific axes:
   - lower effective proliferative ceiling in vivo;
   - stronger stress-linked chromosome missegregation in vivo;
   - more ploidy-dependent post-missegregation survival filter in vivo.
3. Use multi-warmup joint reports to show whether these axes recur across representative starts.

Deferred evidence questions:

- Which joint run is the canonical manuscript run.
- Whether soft-coupling projection/objective semantics need to be settled before this becomes a main result.
- Whether Fig. 4 should show one best joint fit or a family/representative summary.

### Figure 5: Fixed-O2 Resource-Regime Behavior

Claim nodes:

- C5: Dominant ploidy state is resource-regime dependent.
- S1: Steady-state ploidy mode is predictable from fitted parameters, but predictability depends on O2/resource level.
- S2: Low-resource ploidy behavior is mainly death dominated.
- S3: Intermediate-resource behavior is least reducible and most tradeoff-dependent.
- S4: High-resource ploidy behavior is governed by baseline missegregation with buffering modifying outcome.
- S5: Continuous dominant-ploidy analyses support discrete mode analyses.
- S6: Full fixed-O2 ploidy-response curve classes provide a richer phenotype than binary mode labels.

Existing evidence to reuse:

- `/Users/4470246/Downloads/mode_parameter_contribution_report.html`
- `/Users/4470246/Downloads/mode_parameter_contribution_report (1).html`
- `/Users/4470246/Downloads/fixo2_eigen_attractor_embedding_curve_class_report.html`
- `/Users/4470246/Downloads/pooled_embedding_curve_class_report.html`
- Existing code:
  - `oxygen/code/O2_supply_demand_MAP/analysis/parameter_landscape_clustering/mode_parameter_contribution_report.R`
  - `oxygen/code/O2_supply_demand_MAP/analysis/parameter_landscape_clustering/dominant_ploidy_parameter_contribution_report.R`
  - `oxygen/code/O2_supply_demand_MAP/analysis/best_fit_parameter_feature/01_fixed_o2/FixO2_invivo.R`
  - `oxygen/code/O2_supply_demand_MAP/analysis/best_fit_parameter_feature/05_FixO2_eigen_attractor_based_clustering/fixo2_eigen_attractor_report.R`
  - `oxygen/code/O2_supply_demand_MAP/analysis/best_fit_parameter_feature/06_combine_FixO2_eigen_attractor/render_fixo2_eigen_attractor_embedding_curve_class_report.R`

First no-code panel actions:

1. Use the mode-parameter contribution report to extract:
   - AUC-by-O2 summaries;
   - top feature rankings by O2 level;
   - representative low/intermediate/high O2 panels.
2. Use the fixed-O2/eigen attractor report to show dominant ploidy as a function of O2.
3. Use curve-class report outputs to show full O2-response class shapes, not only binary mode.
4. Arrange Fig. 5 around the resource axis:
   - low O2;
   - intermediate O2;
   - high O2;
   - full curve-class summary.

Deferred evidence questions:

- Whether the main figure should emphasize binary mode, continuous dominant ploidy, or full curve classes.
- Whether fixed-O2 attractor results should be paired with finite-time trajectory caveats in the main figure or supplement.

### Figure 6: Parameter Landscape And Candidate Biological Regimes

Claim nodes:

- C6: Multiple fitted solution regions may represent distinct biological regimes rather than purely optimizer noise.
- S6: Full fixed-O2 ploidy-response curve classes provide richer phenotype overlays.
- Methods guardrail: parameter landscapes need prior-aware preprocessing and original/derived feature characterization.

Existing evidence to reuse:

- `/Users/4470246/Downloads/landscape_subcluster.zip`
- `/Users/4470246/Downloads/parameter_landscape_clustering_umap_cluster_report (2).html`
- `/Users/4470246/Downloads/parameter_landscape_clustering_pca_cluster_report (2).html`
- `/Users/4470246/Downloads/parameter_landscape_pooled_pca_interpretation.html`
- `/Users/4470246/Downloads/pooled_embedding_curve_class_report.html`
- `/Users/4470246/Downloads/landscape_subcluster_tsne_invivo_pmisbase_nO_violins.png`
- `/Users/4470246/Downloads/landscape_subcluster_tsne_invivo_pmisbase_nO_violins.pdf`
- Existing code:
  - `oxygen/code/O2_supply_demand_MAP/util/build_multi_warmup_pairs_from_landscape_subclusters.R`
  - `oxygen/code/O2_supply_demand_MAP/analysis/parameter_landscape_clustering/clustering_runner.R`
  - `oxygen/code/O2_supply_demand_MAP/analysis/parameter_landscape_clustering/clustering_report.R`
  - `oxygen/code/O2_supply_demand_MAP/analysis/parameter_landscape_clustering/full_data_in_vivo_clustring.R`
  - `oxygen/code/O2_supply_demand_MAP/analysis/parameter_landscape_clustering/parameter_landscape_utils.R`

First no-code panel actions:

1. Use the landscape-subcluster UMAP/t-SNE cluster figures as the first Fig. 6 landscape panels.
2. Use `seed_cluster_original_fit_prior_positions.csv` to summarize original-parameter differences for candidate regions.
3. Include the existing `p_mis_base` and `n_O` violin plot as the first focused in vivo cluster-driver panel.
4. Use the pooled embedding curve-class report to overlay or reference fixed-O2 response classes on the landscape.

Deferred evidence questions:

- Whether UMAP, t-SNE, or PCA should be the main landscape view.
- Whether in vivo cluster structure is stable enough for a main claim.
- Whether model-implied biological feature characterization needs to be regenerated before the final panel.

## Supplemental And Guardrail Panels

### Fixed-O2 Attractor Versus Finite-Time Trajectories

Existing evidence/code:

- `docs/fixed_o2_ploidy_monotonicity_implementation_plan.md`
- `oxygen/code/O2_supply_demand_MAP/analysis/dense-grid_monotonicity_classification/compare_eigen_initial_ploidy_predictions.R`
- `oxygen/code/O2_supply_demand_MAP/analysis/best_fit_parameter_feature/01_fixed_o2/FixO2_invivo.R`

Panel goal:

Show that fixed-O2 attractors and finite-time experimental trajectories are related but not identical. This should qualify interpretation of mode/attractor analyses.

### Prior-Aware Preprocessing

Existing evidence/code:

- `docs/parameter_landscape_prior_aware_preprocessing_plan.md`
- `oxygen/code/O2_supply_demand_MAP/analysis/parameter_landscape_clustering/parameter_landscape_utils.R`

Panel goal:

Document how parameter transformations, bounds, and prior-unit scaling affect landscape interpretation.

### Ploidy-Independent Division Ablation

Existing evidence:

- Current in vitro fit report and the weak lower-bound `alpha_o2`/`gamma_growth` interpretation.

Panel goal:

This is currently an ablation-design panel, not an evidence panel. It belongs in supplement or future-work unless the ablation is run.

## Priority Order For First Panel Draft

1. Figure 1: use existing optimization-data-stream figures directly.
2. Figure 3: extract in vitro and in vivo fit panels from existing reports.
3. Figure 4: use existing joint ratio plots and multi-warmup reports.
4. Figure 5: extract fixed-O2 mode/AUC/curve-class panels from existing reports.
5. Figure 6: use landscape-subcluster figures, summary table, and existing violin plots.
6. Figure 2: draft mechanism schematic from model components and current claim notes.

This order prioritizes panels where reports or figures already exist, then fills in the conceptual mechanism schematic.

## Later Evidence Audit

After the first panel draft exists, revisit each node with an evidence status:

- `supported_main`
- `supported_supplement`
- `provisional`
- `needs_analysis`
- `drop_or_merge`

This status should not be assigned during the version-0 graph and panel-backlog generation unless the evidence gap is already obvious and noncontroversial.
