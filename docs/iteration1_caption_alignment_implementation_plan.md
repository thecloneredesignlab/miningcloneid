# Iteration 1 Caption Alignment Implementation Plan

## Purpose

Bring `oxygen/figures/iteration1/figureCaptions.txt` into closer alignment with the current intended figure structure in `/Users/4470246/Downloads/hypoxiaLTEE_captions.txt`, while preserving the factual accuracy of the current iteration-1 captions.

The key rule for the alignment pass is:

- If an intended panel has an adequate existing evidence panel, keep the panel and caption it accurately.
- If an intended panel does not yet have evidence, keep the intended caption idea from `hypoxiaLTEE_captions.txt`, but mark the panel as missing in the LaTeX caption using the panel ID, for example `Fig. 2A: missing`.
- If a saved iteration-1 panel is poorly matched to the intended Downloads caption structure, remove it from `oxygen/figures/iteration1/`.
- Do not let a weakly related existing panel substitute for a missing intended panel.

This is a cleanup/alignment plan, not a plan to generate new analyses.

## Source Documents

Structural source of truth:

- `/Users/4470246/Downloads/hypoxiaLTEE_captions.txt`

Accuracy guardrails:

- `oxygen/figures/iteration1/figureCaptions.txt`
- `oxygen/figures/iteration1/panel_selection_log.md`
- Existing panel files in `oxygen/figures/iteration1/`

## Deliverables

1. Revised `oxygen/figures/iteration1/figureCaptions.txt` using the Figure 1-6 structure from `hypoxiaLTEE_captions.txt`.
2. A pruned `oxygen/figures/iteration1/` containing only panels that are aligned with the intended caption structure.
3. Updated `oxygen/figures/iteration1/panel_selection_log.md` recording panels kept, panels removed, and missing intended panels.
4. Validation output confirming:
   - every non-missing panel named in `figureCaptions.txt` exists;
   - removed panels are no longer in `oxygen/figures/iteration1/`;
   - all figure environments and labels in `figureCaptions.txt` are balanced.

## Panel Alignment Inventory

### Figure 1: Matched SUM159 Lineages Reveal Environment-Dependent Ploidy Selection

Downloads structure:

- Figure title only; no detailed panel list.

Current candidate:

- `fig1a_optimization_data_streams_overview.png`

Decision:

- Keep.

Rationale:

The panel is well aligned with the Figure 1 title because it shows the matched SUM159 2N/4N system, in vitro and in vivo data streams, absolute time-frame contrast, and in vitro known-O2 versus in vivo latent-resource distinction.

Caption action:

- Keep as Figure 1A and use the current accurate caption language.

### Figure 2: Resource Limitation Creates Opposing Pressures on Ploidy Evolution

Downloads structure:

- Fig. 2A: Central model schematic.
- Fig. 2B: Cost arm.
- Fig. 2C: Instability arm.
- Fig. 2D: Buffering arm.
- Fig. 2E: Integrated outcome.

Current candidates:

- None.

Decision:

- Keep the intended panel ideas in the caption, but mark all panels missing.

Caption action:

- Add a Figure 2 LaTeX block with:
  - `Fig. 2A: missing. Central model schematic. Resource limitation reshapes ploidy evolution through competing forces.`
  - `Fig. 2B: missing. Cost arm. High ploidy is costly under resource stress through growth suppression and death.`
  - `Fig. 2C: missing. Instability arm. Stress can increase CIN/WGD pressure.`
  - `Fig. 2D: missing. Buffering arm. High ploidy buffers otherwise lethal chromosome missegregation.`
  - `Fig. 2E: missing. Integrated outcome. The outcome depends on growth suppression, death, CIN/WGD generation, and post-missegregation survival filtering.`

File action:

- No files to remove for Figure 2.

### Figure 3: In Vitro Ploidy Reshaping Occurs Without Strong High-Ploidy Elimination

Downloads structure:

- Fig. 3A: In vitro fit summary.
- Fig. 3B: Missegregation source: high-ploidy cells missegregate more often.
- Fig. 3C: Survival buffering.
- Fig. 3D: Downward reshaping through surviving chromosome-loss daughters.
- Fig. 3E: Negative control / rejected interpretation: the story is not "4N cells are killed off."

Current candidates:

- `fig3a_invitro_growth_ploidy_burden_fit.png`
- `fig3b_invitro_predicted_ploidy_distribution.png`
- `fig3c_invitro_viability_after_missegregation.png`
- `fig3d_invitro_nonviable_fraction_vs_ms_rate.png`
- `fig3e_invivo_burden_fit.png`
- `fig3f_invivo_weighted_mean_ploidy_fit.png`

Decisions:

- Keep `fig3a_invitro_growth_ploidy_burden_fit.png` for Fig. 3A.
- Mark Fig. 3B as missing.
- Keep `fig3c_invitro_viability_after_missegregation.png` for Fig. 3C.
- Mark Fig. 3D as missing unless a direct panel showing viable chromosome-loss descendants is found in existing reports.
- Mark Fig. 3E as missing.
- Remove `fig3b_invitro_predicted_ploidy_distribution.png` from iteration1 because it is only a partial proxy for downward reshaping and does not show surviving chromosome-loss daughters.
- Remove `fig3d_invitro_nonviable_fraction_vs_ms_rate.png` because it is mechanistically related but does not show the intended high-ploidy missegregation source panel.
- Remove `fig3e_invivo_burden_fit.png` and `fig3f_invivo_weighted_mean_ploidy_fit.png` because they are in vivo panels and the Downloads Figure 3 structure is explicitly in vitro.

Caption action:

- Rewrite Figure 3 around the intended A-E structure:
  - Fig. 3A: use accurate caption for `fig3a`.
  - Fig. 3B: missing; preserve the intended missegregation-source wording.
  - Fig. 3C: use accurate caption for `fig3c`.
  - Fig. 3D: missing; preserve the intended downward-reshaping wording.
  - Fig. 3E: missing; preserve the rejected-interpretation wording.

### Figure 4: In Vivo Fit-Derived Resource Regimes And Parameter Landscapes

Downloads structure:

- Fig. 4A: Oxygen framing.
- Fig. 4B: Low-O2 regime.
- Fig. 4C: Intermediate-O2 regime.
- Fig. 4D: High-O2 regime.
- Fig. 4E: Parameter landscape.

Current candidates:

- `fig4a_invivo_effective_oxygen_trajectories.png`
- `fig4b_mode_predictability_top3_auc_by_o2.png`
- `fig4c_low_o2_main_feature_auc.png`
- `fig4d_intermediate_o2_main_feature_auc.png`
- `fig4e_high_o2_main_feature_auc.png`
- `fig4f_landscape_tsne_clusters.png`
- `fig4g_landscape_cluster_driver_violins.png`

Decisions:

- Keep `fig4a_invivo_effective_oxygen_trajectories.png` for Fig. 4A, but caption it as partial evidence for oxygen/resource framing, not as the full conceptual oxygen-is-not-resource-stress panel.
- Keep `fig4c_low_o2_main_feature_auc.png` for Fig. 4B.
- Keep `fig4d_intermediate_o2_main_feature_auc.png` for Fig. 4C.
- Keep `fig4e_high_o2_main_feature_auc.png` for Fig. 4D.
- Keep `fig4g_landscape_cluster_driver_violins.png` for Fig. 4E.
- Remove `fig4b_mode_predictability_top3_auc_by_o2.png` from iteration1 because it is useful supporting context but is not part of the current Downloads Figure 4 panel structure.
- Remove `fig4f_landscape_tsne_clusters.png` from iteration1 unless a two-part Fig. 4E is explicitly desired. The Downloads panel idea specifically emphasizes biological cluster drivers such as `p_mis_base` and `n_O`, which are shown more directly by the violin plot.

Caption action:

- Keep the intended Figure 4A-E structure.
- Use the Downloads biological language first, then add accurate clauses describing what the actual retained panel shows.
- Do not overstate the feature-AUC panels as proof of mechanisms; describe them as classifier/feature-importance summaries of fitted parameter landscapes at fixed reference O2.

### Figure 5: Joint Fits Distinguish In Vivo And In Vitro Stress Responses

Downloads structure:

- Fig. 5A: Joint-fit setup.
- Fig. 5B: Proliferation axis.
- Fig. 5C: CIN axis.
- Fig. 5D: Survival-filter axis.
- Fig. 5E: Integrated contrast.

Current candidates:

- `fig5a_joint_parameter_ratio_all_soft.png`
- `fig5b_joint_post_ms_survival_in_vivo_vs_invitro.png`
- `fig5c_joint_o2_vs_missegregation_rate.png`
- `fig5d_joint_o2_vs_proliferation_rate.png`
- `fig5e_joint_o2_vs_death_rate.png`

Decisions:

- Mark Fig. 5A as missing. The current ratio plot is an outcome/contrast panel, not a setup panel.
- Keep `fig5d_joint_o2_vs_proliferation_rate.png` for Fig. 5B if the canonical-run caveat is acceptable; otherwise mark Fig. 5B as missing.
- Keep `fig5c_joint_o2_vs_missegregation_rate.png` for Fig. 5C if the canonical-run caveat is acceptable; otherwise mark Fig. 5C as missing.
- Keep `fig5b_joint_post_ms_survival_in_vivo_vs_invitro.png` for Fig. 5D if the canonical-run caveat is acceptable; otherwise mark Fig. 5D as missing.
- Keep `fig5a_joint_parameter_ratio_all_soft.png` for Fig. 5E because it most closely functions as the integrated context contrast.
- Remove `fig5e_joint_o2_vs_death_rate.png` because the Downloads Figure 5 structure does not include a separate death-rate panel.

Caption action:

- Explicitly state that panels B-D are provisional if retained from `fit_joint_sigma1e6_fit_report_seed38.html`.
- Avoid implying that the ratio plot and mechanism-function panels come from the same canonical joint run unless that is confirmed.
- If run consistency is prioritized over retaining provisional evidence, mark Fig. 5B-D as missing and keep only Fig. 5E.

### Figure 6: Oxygen-Linked Joint-Fit Behaviors Guide Model Selection

Downloads structure:

- Fig. 6A: Fixed-O2 response classes.
- Fig. 6B: Candidate behavior shapes.
- Fig. 6C: Monotone-increasing priority class.
- Fig. 6D: Selection overlay.

Current candidates:

- `fig6a_pooled_fixed_o2_curve_class_examples.png`
- `fig6b_objective_distribution_by_curve_class.png`

Decisions:

- Keep `fig6a_pooled_fixed_o2_curve_class_examples.png` as the shared evidence panel for Fig. 6A-C, because it shows fixed-O2 response classes, candidate shapes, and includes the monotone-increasing class.
- Mark Fig. 6D as missing.
- Remove `fig6b_objective_distribution_by_curve_class.png` from iteration1 if strict structural alignment is prioritized, because it is not the requested selection overlay. If retained, it should be labeled as supporting context rather than Fig. 6D.

Recommended implementation choice:

- For strict alignment, remove `fig6b_objective_distribution_by_curve_class.png` and mark Fig. 6D missing.
- For evidence usefulness, keep it outside the main Figure 6 caption as a candidate supplement. Because the user requested removing poorly matched panels from `iteration1`, use the strict-alignment choice in this pass.

Caption action:

- Figure 6A-C should accurately describe the one retained curve-class panel.
- Figure 6D should preserve the intended overlay language and be marked missing.

## File Removal Plan

Remove these poorly matched or extra panels from `oxygen/figures/iteration1/` during the implementation pass:

- `fig3b_invitro_predicted_ploidy_distribution.png`
- `fig3d_invitro_nonviable_fraction_vs_ms_rate.png`
- `fig3e_invivo_burden_fit.png`
- `fig3f_invivo_weighted_mean_ploidy_fit.png`
- `fig4b_mode_predictability_top3_auc_by_o2.png`
- `fig4f_landscape_tsne_clusters.png`
- `fig5e_joint_o2_vs_death_rate.png`
- `fig6b_objective_distribution_by_curve_class.png`

Keep these aligned panels:

- `fig1a_optimization_data_streams_overview.png`
- `fig3a_invitro_growth_ploidy_burden_fit.png`
- `fig3c_invitro_viability_after_missegregation.png`
- `fig4a_invivo_effective_oxygen_trajectories.png`
- `fig4c_low_o2_main_feature_auc.png`
- `fig4d_intermediate_o2_main_feature_auc.png`
- `fig4e_high_o2_main_feature_auc.png`
- `fig4g_landscape_cluster_driver_violins.png`
- `fig5a_joint_parameter_ratio_all_soft.png`
- `fig5b_joint_post_ms_survival_in_vivo_vs_invitro.png`
- `fig5c_joint_o2_vs_missegregation_rate.png`
- `fig5d_joint_o2_vs_proliferation_rate.png`
- `fig6a_pooled_fixed_o2_curve_class_examples.png`

## Caption Rewrite Plan

Rewrite `oxygen/figures/iteration1/figureCaptions.txt` into six LaTeX figure blocks:

1. Figure 1: one retained panel.
2. Figure 2: all intended panels marked missing.
3. Figure 3: retained A/C; B/D/E marked missing.
4. Figure 4: retained A-E, using the selected evidence panels.
5. Figure 5: A marked missing; B-D retained with provenance caveat; E retained as integrated contrast.
6. Figure 6: A-C retained as one curve-class evidence panel; D marked missing.

Use this missing-panel convention:

```latex
\textbf{(B) Fig. 3B: missing. Missegregation source.}
High-ploidy cells are expected to missegregate more often, but no iteration-1 panel directly visualizes this quantity.
```

Use this retained-panel convention:

```latex
\textbf{(C) Survival buffering.}
Panel file: \texttt{fig3c_invitro_viability_after_missegregation.png}. The existing report panel shows model-implied viability after a one-copy missegregation event as a function of chromosome number...
```

## Validation Plan

After the implementation pass:

1. List files:

```bash
find oxygen/figures/iteration1 -maxdepth 1 -type f | sort
```

2. Validate retained PNGs:

```bash
file oxygen/figures/iteration1/*.png
```

3. Validate captions:

```bash
python3 - <<'PY'
from pathlib import Path
text = Path('oxygen/figures/iteration1/figureCaptions.txt').read_text()
assert text.count('\\\\begin{figure}') == text.count('\\\\end{figure}')
assert text.count('\\\\label{fig:') == text.count('\\\\begin{figure}')
print('figure environment validation passed')
PY
```

4. Validate all non-missing panel filenames mentioned in caption comments exist. Missing panels should be explicit text entries, not missing file references.

5. Review `git status --short` and confirm only intended files changed or were removed.

## Open Decision Before Implementation

The main implementation choice is how strict to be with partially matched panels.

Recommended default:

- Be strict for `iteration1`: remove weak proxies and mark intended panels missing.

Reason:

The point of this pass is alignment to the current manuscript figure structure. Weak proxies make the folder look more complete than it is and can blur the difference between available evidence and intended figure content.
