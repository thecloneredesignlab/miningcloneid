# Round-3 mandatory analyses

This directory is an isolated, reproducible addendum for the round-3 manuscript revision. It does not overwrite the frozen figures, fitted objects, scripts, or manuscript in `revised/iteration1`.

## Scope

- `review_triage.md`: adjudication of the two review packages and the deep-research report.
- `scripts/`: deterministic local analyses and manuscript assembly.
- `results/`: machine-readable outputs grouped by review issue.
- `figures/` and `tables/`: manuscript-ready additions.
- `provenance/`: file hashes, software/session records, and input manifests.
- `hpc/`: sensitivity/refit task matrices only. No HPC job is submitted without separate authorization.
- `logs/`: captured local run logs.

## Evidence boundary

All new numerical analyses use existing data or fitted outputs already present in this repository or in the adjacent `soft_coupling` checkout. Analyses that require refitting are explicitly labeled pending and are not represented as completed evidence.

## Reproduction

From the repository root:

```bash
Rscript revised/add/round3_mandatory/scripts/00_freeze_input_manifest.R
Rscript revised/add/round3_mandatory/scripts/01_o1_o2_endpoint_adequacy.R
Rscript revised/add/round3_mandatory/scripts/02_full_sample_predictive_adequacy.R
Rscript revised/add/round3_mandatory/scripts/03_anchor_and_joint_function_audit.R
Rscript revised/add/round3_mandatory/scripts/04_fixed_o2_finite_time_validation.R
Rscript revised/add/round3_mandatory/scripts/05_boundary_necrosis_wgd_audit.R
Rscript revised/add/round3_mandatory/scripts/05b_expanded_grid_boundary_sensitivity.R
Rscript revised/add/round3_mandatory/scripts/06_build_integrated_manuscript.R
```

The integrated manuscript is generated as:

`revised/iteration1/manuscript/ltee_hypoxia_model_round3_integrated.tex`

For a temporary compile check from `revised/iteration1/manuscript`:

```bash
qa_dir=$(mktemp -d /tmp/ltee_round3_pdf.XXXXXX)
latexmk -pdf -interaction=nonstopmode -halt-on-error \
  -outdir="$qa_dir" ltee_hypoxia_model_round3_integrated.tex
```

## Decision records

- `review_triage.md`: review-by-review adjudication.
- `final_recommendation.md`: current submission-readiness decision and remaining required work.
- `hpc/`: unsubmitted refit and held-out-validation task specifications.
