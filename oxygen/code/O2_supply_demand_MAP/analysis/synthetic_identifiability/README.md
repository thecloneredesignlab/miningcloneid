# Synthetic recovery and experimental precision

| Workflow | Source directories | Role |
| --- | --- | --- |
| Production contracts and D0 | `sif_contract_validation/`, `culture_d0_recovery/`, `culture_conditional_d0_800/` | Frozen production replay, parameter mappings and conditional deterministic diagnostics. |
| Conditional growth and karyotype design | `culture_growth_karyotype_conditional/`, `culture_growth_karyotype_recovery/` | Generate and validate the production-compatible synthetic observations. |
| D19 deterministic recovery | `culture_D19_feasibility/`, `culture_D19_recovery_restart/`, `culture_D19_independent/` | Feasible blind starts, recovery and independent replay. |
| N20 noisy fitting | `culture_N20_independent/` | Existing engineering, cost checks, frozen release recovery, scoring and collection. The full-pilot collector also uses `culture_D19_independent/collect_closeout.py`. |
| Experimental OLS precision | `experimental/ols/` | Fits, likelihood profiles and recovery with the 114 OLS slopes. |
| Direct and correlated counts | `experimental/counts/`, `experimental/counts_ar1/` | Fits, profiles, recovery and compensation with the 221 post-seeding counts. The AR(1) workflow uses within-passage correlation. |
| Karyotype and death readouts | `experimental/karyotype_and_death/` | Assay alignment, A17 joint fitting, and comparison of matched death fractions with A7K. |
| A7K | `experimental/a7k/` | Inheritance perturbations and branch, ploidy-cohort, and individual 4N swaps. |

For experimental study-source mappings and local archival replay, see [the structural analysis index](../structural_identifiability/structural_first/README.md). The N20 working source, frozen release and checkpoints keep their established paths.
