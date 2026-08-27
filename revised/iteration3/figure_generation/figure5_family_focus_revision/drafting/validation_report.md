# Figure 5/6 legend and narrative validation

## Scope

- Starting manuscript: `revised/iteration2/manuscript/ltee_hypoxia_model_results_redraft_iteration3.tex`.
- Iteration-3 manuscript: `revised/iteration3/manuscript/ltee_hypoxia_model_results_redraft_iteration3.tex`.
- Requested figure changes: retain the exact approved Supplementary Figure 3 display as Figure 5D and revise the Figure 5 legend without regenerating the other figures.
- Results scope: preserve all three families, focus biological interpretation on C01 with an explicit operational rationale, distinguish cross-family conclusions from C01-specific parameter equality, and state missegregation quantities at the level represented by the model.
- Ecological interpretation scope: distinguish imposed in-vitro O2 from latent fitted in-vivo O2, treat common O2 values as standardized model-coordinate comparisons, and explain the high-O2 Figure 6 reversal without claiming that directional chromosome-state flux was quantified.
- A pre-edit manuscript checkpoint was committed as `6c30436` before the ecological and fixed-input interpretation revisions.

## Passed

- The iteration-3 Figure 5 PNG is byte-identical to the reviewed iteration-2 composite (`SHA-256: 1dbd0db6e9fd090cce6bb67c49906cf78fa95a7d748d04465e854698569ac785`).
- The standalone Supplementary Figure 3 PNG is also byte-identical across iterations (`SHA-256: 23ff0bcc38a4046eaab3a98f2d7abdc26055df7802d6f30b55461e80210a58ec`).
- Original-resolution inspection confirmed that Figure 5D contains the C01--C03 endpoint-class matrix and aligned endpoint-level context-ratio distributions, with no clipping or overlap.
- The Figure 5 Results subsection follows the results-text structure: question and fit-family scope, shared cross-family conclusions, C01-specific interpretation, function-level evidence, and a consolidated limitation and synthesis.
- C01 is described using the precise modeled quantities: approximately shared baseline per-chromosome missegregation and response steepness, but higher complete effective per-chromosome missegregation in vivo at standardized low O2 through the death-hazard-linked inducible branch.
- The Figure 5 Results and Discussion distinguish shared cellular machinery from different resource ecologies and do not treat equal numerical O2 inputs as biologically equivalent culture and tumor environments.
- The Figure 6 Results identify recurrent production and retention of chromosome-loss descendants as an encoded, model-consistent route to the high-O2 reversal while explicitly stating that directional transition fluxes were not decomposed.
- Quantitative C01 function values were checked against `c01_function_value_checks.tsv`; cross-family ranges were checked against `tao_joint_function_contrasts.tex`; family objectives and search accounting were checked against `supp_joint_search_summary.tex`.
- The canonical and explicitly named iteration-3 TeX sources are byte-identical.
- `scripts/agentRrunner.sh --check` passed with R 4.5.1.
- R parse checks passed for `draw_Figure5.R`, `draw_Supp_Figure5_1.R`, `data_Supp_Figure5_1.R`, and `run_all_figures.R`.
- `latexmk -pdf -interaction=nonstopmode -halt-on-error ltee_hypoxia_model_results_redraft_iteration3.tex` completed successfully and produced a 67-page PDF with no unresolved citations or references. Existing typesetting warnings remain in inherited tables and long Methods material, but the new 24-point Discussion overflow was removed.
- Rendered inspection of Results pages 7--10 confirmed that the complete Figure 5 and Figure 6 subsections are readable, the C01 focus precedes the function-level interpretation, and the high-to-low replenishment interpretation is bounded by the flux-analysis limitation.
- Rendered inspection of Discussion pages 17--18 confirmed that the ecological interpretation and survival/WGD comparison read continuously without clipping or overlap.

## Regeneration blocker

A fresh data-driven vector rebuild cannot be completed from the local archive. The following generated Supplementary Figure 3 tables are absent:

1. `within_pair_parameter_stability.tsv`
2. `between_pair_parameter_stability.tsv`
3. `analysis_config.tsv`
4. `selected_primary_family_pairs.tsv`
5. `soft_coupling_master_long.tsv`

The builder also requires the 500 per-seed `joint_soft_coupling.tsv` files for each displayed primary family. The available consume-only bundle does not contain those 1,500 context-specific endpoint tables. The Figure 5 worker fails explicitly on these missing inputs; that failure was not counted as a successful test.

## Layout validation

The inherited 0.90-text-height placement clipped the expanded legend. Figure 5 is now placed at 0.62 text height with a locally small caption font. The full figure and legend render together on page 15, and the prior `Float too large for page` warning is absent. Figure 6 and its legend render on page 16. This changes only manuscript placement; the full-resolution assets remain unchanged.
