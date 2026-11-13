# Tumor volume-to-cell conversion: identifiability review

Reviewed on 2026-11-09. Recommendation only: add a dedicated `rho_2N` assessment to the later tumor phase. The current culture pilot remains limited to `p_misseg` and `buffer_beta`. No tumor simulation, optimization or job submission was performed for this review.

## Parameter and observation equation

The parameter is **`rho_2N`**, in cells/mm^3, optimized as `log10_rho_2N`. Its reciprocal, `c_vol_2N_eff_mm3 = 1/rho_2N`, is an equivalent effective volume scale, not a second independent parameter. In the inspected production source, the same reciprocal density is assigned to every chromosome state. Despite the historical parameter name, the current conversion does not vary with ploidy.

Write live abundance as x_N, retained hypoxia-death material as d_H,N, and the retained buffer-death stock under production bookkeeping as d_E,N. The buffer-death label does not itself distinguish biological failed daughters from numerical losses booked into that stock. Under this observation convention:

\[
V(t)=\frac{\sum_N[x_N(t)+d_{H,N}(t)+d_{E,N}(t)]}{\rho_{2N}},
\qquad
f_D(t)=\frac{\sum_N[d_{H,N}(t)+d_{E,N}(t)]}{\sum_N[x_N(t)+d_{H,N}(t)+d_{E,N}(t)]}.
\]

Thus `rho_2N * V` is total model cell-equivalent material, including retained dead stocks. It is not the viable-cell count. At a matched time, the model identity for live abundance is `(1 - f_D) * rho_2N * V`. This identity is conditional on the production volume and dead-stock conventions; it does not turn histology into a direct viable-cell count. The effective density can also absorb discrepancies between this simplified observation model and real tissue composition.

The old `beta_size` entry must not be profiled as an active parameter of this production model. The volume helper returns a constant `1/rho_2N`; the tumor backend supplies `beta_size=0`, and the generator explicitly ignores that argument. A new ploidy-dependent volume law would be a separate model alternative requiring its own data justification and validation.

## Verified production evidence

The authority is the iteration5 Figure 4 seed25 package, not generic defaults in the current worktree. The archived source and runtime config MD5 values match the Figure 4 provenance record. The core R/C++ MD5 values also match that record and the frozen production-v2 source already retained for culture. Supporting helpers differ between the current and frozen source, so this review inspected the frozen conversion functions and tumor backend directly.

| Item | Archived seed25 value |
| --- | --- |
| Fitted `rho_2N` | 141109.43458562 cells/mm^3 |
| Search bounds | 50000–500000 cells/mm^3 |
| Base simulated initial abundance | 1000000 |
| Harvest-specific initial multiplier enabled | FALSE |
| Oxygen-demand reference `o2_Nref` | 1000000 |
| Oxygen feedback | TRUE |
| Crowding | FALSE |
| Soft prior enabled; weight | TRUE; 0.03 |
| Prior center for log10 density | 4.626668957 |
| Prior SD for log10 density | 0.15 |
| Prior center on natural scale | 42332.0165261 cells/mm^3 |

The prior center is below the production lower bound. This is a provenance/configuration finding, not an automatic claim of an error: retain the exact setting for replay and evaluate its influence separately. At the archived endpoint, the density's weighted prior term is approximately 0.182274; the total prior contribution across parameters is 0.579469. The total objective is 14.0508154 and the data contribution is 13.4713464. A fitted endpoint alone cannot establish a narrow data-supported density range.

The read-only audit exercised the two pure conversion functions from frozen source over all 133 chromosome states, checked constant state volumes and inverse-density scaling, and checked the live-count accounting identity on a toy stock vector. These checks are not a tumor trajectory replay or a structural-identifiability proof.

## Why a dedicated assessment is needed

1. **Absolute cellular output depends on the scale.** Similar volume predictions could imply different absolute live counts, dead material and daughter-production totals if `rho_2N` changes and the dynamical parameters compensate. Assess this before interpreting tumor predictions in absolute cell units. Stable per-division quantities may coexist with uncertain absolute counts.
2. **The baseline assumes an absolute initial abundance.** The production configuration fixes initial abundance at one million and disables harvest-specific multipliers. Whether that represents a known effective founding population, rather than nominal inoculation, needs an explicit observation/conditioning audit. A profile conditional on that assumption should be labeled accordingly.
3. **Resource and death parameters can substitute for scale.** Refit active growth/death, oxygen supply/demand, clearance, initial-state and observation-error parameters when profiling density. Do not freeze them using low old Sobol sensitivity. For fixed dynamics, common density factors cancel from viable chromosome fractions and necrotic fractions; these streams can constrain density indirectly through refitted dynamics, not through a direct fractional scale measurement.
4. **The old Sobol study does not answer this question.** Its 14-factor fixed-oxygen operator study excluded density and used normalized mean ploidy and a per-day eigenvalue. It provides no completed density profile, posterior or confidence interval.

## Proposed tumor-stage sequence

After the culture pilot and a separately approved tumor compute plan:

1. Freeze and replay the actual tumor production endpoint without optimization. Check volume, terminal chromosome and necrosis observations, times, initial stocks, clearance, density bounds and each penalty contribution. Do not reuse the culture unrestricted baseline as a tumor baseline.
2. Fit separate unrestricted baselines for the **data-only score** and the **exact production penalty policy**. Profile `rho_2N` over its declared positive production domain, reoptimizing all active nuisance parameters at each point with multiple starts. Report boundary and search limits, full predictions and separate score/penalty terms. The log10 coordinate is appropriate here; no exact-zero density null is meaningful.
3. Audit the fixed initial abundance and fixed demand-reference assumptions. First distinguish externally measured quantities from modeling conventions. If an initial-size sensitivity is warranted, introduce it as an explicit model/configuration extension with its own unrestricted baseline, not as an unannounced change to production. Compare inferred viable abundance, retained dead material, oxygen histories and integrated cell output across good fits.
4. If the density profile indicates important compensation, consider a limited joint density/initial-size or density/resource diagnostic before broad tumor profiles. Integrate the finding with the planned `mu_hp`/`k_clear` study. Stream removal, expanded density bounds and additional targets remain separately budgeted follow-ups. If an absolute scale remains unresolved, prioritize a matched independent count/density calibration with explicit viable, dead and tissue-volume definitions.

A useful **candidate scaling audit**, rather than a conclusion about the current fit, is simultaneous scaling of all live/dead stocks and initial abundance by a, `rho_2N` by a, and the oxygen-demand reference by a. In the ideal homogeneous count equations with crowding disabled and absolute numerical thresholds inactive, this preserves volume, normalized chromosome outputs, necrotic fractions and the live-demand/reference ratio. The production model fixes initial abundance and `o2_Nref`, so this is not a symmetry obtained by changing density alone in the present parameter space. Verify any claimed symmetry against the full numerical implementation, its floors/clipping, and the actual retained outputs before making a structural claim. Bounds and priors can also break an apparent symmetry without supplying experimental information.

For the later five-class COMBOS model, explicitly compare a known conversion scale against an unknown shared conversion scale. A symbolic model that fixes all volume weights cannot answer whether density itself is identifiable. Preserve retained-dead-material accounting and declare initial-abundance/resource assumptions. Conclusions remain conditional on that reduced model.

No formal likelihood-ratio thresholds, p values or confidence intervals follow from these proposed profile-objective diagnostics without calibration.

## Source and output locations

Read-only paper reference: `/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration5/data/Figures/Figure4/`, especially `source_fit_config/fit_config.rds`, `runtime_fit/seed25/{best_params.tsv,fit_summary.tsv,fit_config.rds}`, `invivo_parameter_table_seed25.csv` and `figure4a_runtime_config_provenance.tsv`.

Production tumor fit root recorded there: `/share/lab_crd/taoli/Project/soft_couping_org/oxygen/results/fit_invivo_unified_np256_500seed_all_xxlarge_r442_exact_20260828_145253`. The compact iteration5 `invivo_best_parameters_500seeds.tsv` omits `rho_2N`; it cannot supply the complete density/scale start vectors. Before a tumor multistart study, inspect the full seed packages and collect missing full endpoints from that source if required. Their current completeness has not been established by this review.

Current audit: `oxygen/results/structural_identifiability/20261109_tumor_volume_scale_review/volume_scale_audit.json` and `source_manifest.json`. Reproduce the pure-function/settings audit with `scripts/agentRrunner.sh oxygen/code/O2_supply_demand_MAP/analysis/structural_identifiability/audit_tumor_volume_scale.R <project> <iteration5> <output-directory>`. No production files or current culture fitting inputs were modified.
