# Drafting feedback intake

## Inputs and scope

| input | role |
|---|---|
| `ltee_hypoxia_model.tex` | Restored manuscript and historical caption context; read-only in Gate B |
| `agent-dev/manuscript_work_packages/ltee_hypoxia_model/ideation/ideation_report.html` | Approved Gate A review surface |
| `agent-dev/manuscript_work_packages/ltee_hypoxia_model/ideation/ideas.md` | Approved five-figure drafting contract |
| `agent-dev/manuscript_work_packages/ltee_hypoxia_model/ideation/feedback_decisions.csv` | Stable directive ledger |
| `agent-dev/manuscript_work_packages/ltee_hypoxia_model/ideation/source_inventory.csv` | Source and generator inventory |
| `agent-dev/manuscript_work_packages/ltee_hypoxia_model/ideation/tex_figure_requirements.csv` | Panel roles and claim guardrails |
| `agent-dev/manuscript_work_packages/ltee_hypoxia_model/ideation/fit_quality_evidence.csv` | Fit-quality sources and limitations |
| `agent-dev/manuscript_work_packages/ltee_hypoxia_model/ideation/feedback_manager_context.md` | Canonical feedback-manager handoff |

Package decision: `draft` Figures 1–5. Figure 6 is
`defer_to_future_user_decision`. Production manuscript edits are
`defer_to_manuscript_assembly`.

## Directive coverage plan

| directive_id | drafting instruction | evidence/output | final status |
|---|---|---|---|
| FD01 | Add chromosome trajectories to the matched-design figure. | `final_figures/recommended/figure1.png` | addressed |
| FD02 | Retain one integrated model overview. | `final_figures/recommended/figure2.png` | addressed |
| FD03 | Limit Figure 3A to growth rate over lineage passage. | `final_figures/recommended/figure3.png` | addressed |
| FD04 | Retain the historical control-versus-deprived comparison. | `final_figures/recommended/figure3.png` | addressed |
| FD05 | Do not draft Figure 3E or run a restricted-model comparison. | `not_drafted.md`; no Figure 3E artifact | addressed |
| FD06 | Lead Figures 3–5 with direct fit-quality evidence. | Figures 3A, 4A, and 5A in `final_figures/recommended/` | addressed |
| FD07 | Treat historical Figure 4A–E as in-vivo analyses. | `subagent_notes/figure4.md`; `subagent_notes/figure4_legend.md` | addressed |
| FD08 | Keep the fixed-O2 triptych and solution landscape in the main figure. | `final_figures/recommended/figure4.png` | addressed |
| FD09 | Use exact low/high-O2 parameter definitions and noncausal AUC language. | Figure 4C labels and `subagent_notes/figure4_legend.md` | addressed |
| FD10 | Use the six July joint-pair winners throughout Figure 5. | `final_figures/recommended/figure5_joint_fit_adequacy_and_context_functions.png`; `source_manifest.csv` | addressed |
| FD11 | Order Figure 5 as fit quality, ratios, proliferation, missegregation, survival. | Figure 5A–E | addressed |
| FD12 | Retain a full all-six parameter-ratio view. | Figure 5B | addressed |
| FD13 | Leave Figure 6 out. | `not_drafted.md`; no Figure 6 artifact | addressed |
| FD14 | Defer journal-specific constraints until polishing. | `not_drafted.md`; repository style contract applied | deferred |
| FD15 | Retain `vi_C01`, `vi_C02`, and `vi_C03`. | Figure 4D–E | addressed |
| FD16 | Preserve the pooled in-vivo/in-vitro embedding and saved geometry. | Figure 4D; pinned raster recorded in `source_manifest.csv` | addressed |
| FD17 | Keep main fit blocks compact; place optimizer diagnostics supplementary. | Figures 3A, 4A, 5A; `supplementary/fit_quality_optimizer_diagnostics.png` | addressed |
| FD18 | Draft only within the approved five-figure/exclusion boundary. | Five recommended figures plus declared supplementary diagnostic | addressed |
