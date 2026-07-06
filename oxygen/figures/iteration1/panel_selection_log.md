# Iteration 1 Panel Selection Log

This folder now contains the strict-alignment iteration-1 panel set. The alignment target is `/Users/4470246/Downloads/hypoxiaLTEE_captions.txt`; `figureCaptions.txt` preserves those intended panel ideas while marking unsupported panels as missing.

## Alignment Rule

- Keep panels that directly match the intended Downloads caption structure.
- Rename retained files so filename letters match intended panel letters.
- Remove weak proxy panels from `oxygen/figures/iteration1/`.
- In `figureCaptions.txt`, preserve intended missing panel text and mark it as `Fig. <N><letter>: missing`.
- Do not use a weakly related plot as evidence for an intended panel.

## Active Retained Panels

### Figure 1

- `fig1a_optimization_data_streams_overview.png`

Rationale:

The panel shows the matched SUM159 2N/4N in vitro and in vivo data streams, the time-frame contrast, and known in vitro versus model-implied in vivo O2/resource framing.

### Figure 2

No panels retained.

Missing intended panels:

- Fig. 2A: Central model schematic.
- Fig. 2B: Cost arm.
- Fig. 2C: Instability arm.
- Fig. 2D: Buffering arm.
- Fig. 2E: Integrated outcome.

Reason:

No existing report-generated or code-generated schematic adequately represents these mechanism panels.

### Figure 3

- `fig3a_invitro_growth_ploidy_burden_fit.png`
- `fig3c_invitro_viability_after_missegregation.png`

Rationale:

Fig. 3A directly supports the in vitro fit-summary panel. Fig. 3C directly supports the survival-buffering panel by showing model-implied viability after one-copy missegregation increasing with chromosome number.

Missing intended panels:

- Fig. 3B: Missegregation source, high-ploidy cells missegregate more often.
- Fig. 3D: Downward reshaping through viable chromosome-loss descendants.
- Fig. 3E: Negative control / rejected interpretation that the in vitro story is not simply 4N cells being killed.

### Figure 4

- `fig4a_invivo_effective_oxygen_trajectories.png`
- `fig4b_low_o2_main_feature_auc.png`
- `fig4c_intermediate_o2_main_feature_auc.png`
- `fig4d_high_o2_main_feature_auc.png`
- `fig4e_landscape_cluster_driver_violins.png`
- `fig4f_landscape_tsne_clusters.png`

Rationale:

Fig. 4A is a partial but accurate model-implied in vivo O2/resource panel. Fig. 4B-D are fixed-reference-O2 feature-AUC panels that map onto low, intermediate, and high O2 resource-regime interpretations. Fig. 4E directly supports the parameter-landscape idea by showing `p_mis_base` and `n_O` differences across in vivo clusters. Fig. 4F restores the t-SNE embedding of the pooled initial-sample cloud and best-fit in vivo/in vitro solutions, with landscape subclusters overlaid, so the reader can see where the driver-defined clusters sit in the fitted landscape.

Rename actions:

- `fig4c_low_o2_main_feature_auc.png` -> `fig4b_low_o2_main_feature_auc.png`
- `fig4d_intermediate_o2_main_feature_auc.png` -> `fig4c_intermediate_o2_main_feature_auc.png`
- `fig4e_high_o2_main_feature_auc.png` -> `fig4d_high_o2_main_feature_auc.png`
- `fig4g_landscape_cluster_driver_violins.png` -> `fig4e_landscape_cluster_driver_violins.png`

### Figure 5

- `fig5b_joint_o2_vs_proliferation_rate.png`
- `fig5c_joint_o2_vs_missegregation_rate.png`
- `fig5d_joint_post_ms_survival_in_vivo_vs_invitro.png`
- `fig5e_joint_parameter_ratio_all_soft.png`

Rationale:

Fig. 5B-D are the closest existing mechanism-function panels for the proliferation, CIN, and survival-filter axes. Fig. 5E is the broadest available integrated joint context contrast.

Missing intended panel:

- Fig. 5A: Joint-fit setup.

Caveat:

Fig. 5B-D come from `fit_joint_sigma1e6_fit_report_seed38.html`, whereas Fig. 5E comes from the all-soft ratio figure. These should be regenerated from a single canonical joint run before final manuscript use.

Rename actions:

- `fig5d_joint_o2_vs_proliferation_rate.png` -> `fig5b_joint_o2_vs_proliferation_rate.png`
- `fig5b_joint_post_ms_survival_in_vivo_vs_invitro.png` -> `fig5d_joint_post_ms_survival_in_vivo_vs_invitro.png`
- `fig5a_joint_parameter_ratio_all_soft.png` -> `fig5e_joint_parameter_ratio_all_soft.png`

### Figure 6

- `fig6a_pooled_fixed_o2_curve_class_examples.png`

Rationale:

This single panel shows fixed-O2 response classes, candidate behavior shapes, and the monotone-increasing class, so it supports Fig. 6A-C.

Missing intended panel:

- Fig. 6D: Selection overlay combining response class, joint-fit parameter landscape, fit quality, and biological plausibility.

## Removed Panels

The following weakly matched or extra panels were removed from `oxygen/figures/iteration1/` during the strict caption-alignment pass:

- `fig3b_invitro_predicted_ploidy_distribution.png`
- `fig3d_invitro_nonviable_fraction_vs_ms_rate.png`
- `fig3e_invivo_burden_fit.png`
- `fig3f_invivo_weighted_mean_ploidy_fit.png`
- `fig4b_mode_predictability_top3_auc_by_o2.png`
- `fig5e_joint_o2_vs_death_rate.png`
- `fig6b_objective_distribution_by_curve_class.png`

Reason:

Each removed panel was useful as supporting context or exploratory evidence, but did not directly match the current Downloads caption structure.

## Restored Historical Panels

On 2026-07-03, the historical first-pass panels requested by name were regenerated from the existing reports. The two panels that now align with the manuscript-facing in vitro figure were restored under their original filenames and included in the rebuilt assembled Fig. 3. The historical in vivo burden panel was evaluated and documented below, but is not retained in the current iteration folder because it does not match the intended Fig. 3E negative-control panel.

### `fig3b_invitro_predicted_ploidy_distribution.png`

Selected source:

- `/Users/4470246/Downloads/reports/in_vitro/fit_report_seed10.html`
- Embedded report figure: `Figure 2.7 Predicted Ploidy Distribution`

Candidates considered:

- `Figure 2.5 Aligned Growth, Chromosome Count, and Burden Fit`
- `Figure 2.6 Flow-Density Fit`
- `Figure 2.7 Predicted Ploidy Distribution`

Decision rationale:

`Figure 2.7` was selected because it most directly shows the model-predicted chromosome-state distribution across in vitro passages, including the deprived 2N high-ploidy/WGD expansion and later broadening/downward reshaping. The aligned fit and flow-density panels are useful validation views but are less direct for the distribution-trajectory question.

### `fig3d_invitro_nonviable_fraction_vs_ms_rate.png`

Selected source:

- `/Users/4470246/Downloads/reports/in_vitro/fit_report_seed10.html`
- Embedded report figure: `Figure 3.1 Nonviable Daughter Fraction vs MS Rate`

Candidates considered:

- `Figure 3.1 Nonviable Daughter Fraction vs MS Rate`
- `Figure 3.2 Death Rate vs Missegregation Rate`
- `Figure 3.3 Ploidy vs Viability After MS`

Decision rationale:

`Figure 3.1` was selected because it directly links missegregation rate to the fraction of nonviable daughter cells across reference ploidy states. It is mechanistically related to the intended missegregation/death counterforce panel. The death-rate panel is downstream of multiple mechanisms, and the viability panel was already retained separately as `fig3c`.

### Historical candidate: `fig3e_invivo_burden_fit.png`

Selected source:

- `/Users/4470246/Downloads/reports/in_vivo_top_10/rank01_fit_report_seed146.html`
- Embedded report figure: `Figure 1.1 Burden Trend Absolute (Real Scale)`

Candidates considered:

- `Figure 1.1 Burden Trend Absolute (Real Scale)`
- `Figure 1.2 Burden Live/Dead Decomposition`
- `Figure 1.3 Predicted Burden Live/Dead Decomposition (0-100/300/1000 day)`

Decision rationale:

`Figure 1.1` was selected because it is the clearest standalone in vivo burden-fit panel, showing observed tumor burden and fitted trajectories for the seed-146 separate in vivo fit. The live/dead decomposition panels contain useful mechanistic detail but are less direct as a compact burden-fit candidate.

Caveat:

This historical `fig3e` candidate is not the currently intended Figure 3E negative-control/rejected-interpretation panel described in `figureCaptions.txt`. It is therefore not retained in the current iteration folder and is not included in the rebuilt `assembled_fig3.png`.

## Generated Figure-Level Outputs

The assembly script `oxygen/figures/assemble_iteration_panels.py` creates figure-level PNGs from the retained panel files and writes them directly under `oxygen/figures/`:

- `oxygen/figures/assembled_fig1.png`
- `oxygen/figures/assembled_fig3.png` (current rebuild contains Fig. 3A-D; Fig. 3E remains an intended panel to be generated)
- `oxygen/figures/assembled_fig4.png`
- `oxygen/figures/assembled_fig5.png`
- `oxygen/figures/assembled_fig6.png`

No `assembled_fig2.png` is expected because all Figure 2 panels are currently missing.
