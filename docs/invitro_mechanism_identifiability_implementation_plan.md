# Structural Identifiability and Mechanism Tests: Culture First, Then Tumors

**Status:** Draft for review; no implementation, fitting, or symbolic analysis has been run under this plan.
**Date:** 2026-09-23.
**Scope:** Two complementary reductions, starting with culture and extending to tumors and cross-context comparisons in explicitly gated stages. Preserve the current manuscript, figures, production fits, and data. eFAST remains a separate analysis.

**Review update:** The five-class recurrent symbolic model, retained dead biomass, extended model family, tumor observation design, and measurement-stream removal tests below are required parts of the staged plan. The one-generation ledger is an accounting check, not a substitute for the five-class model. The existing filename is retained for continuity.

## 1. Purpose and Decisions

The manuscript proposes that resource stress changes chromosome-variation generation, while ploidy changes which altered daughters survive and expand. The guiding question is:

> Can the observations distinguish selection against chromosome-rich cells from the generation and survival of chromosome-loss descendants, and can they distinguish how those processes operate in culture versus tumors?

The immediate aim is not to estimate every parameter independently. It is to determine which mechanisms and effective processes the observations constrain, and whether the biological conclusions survive simpler alternatives.

Pursue two complementary tracks:

1. **Full chromosome-state model family:** Retain the full grid and observation models and refit targeted alternatives. Begin with the reference culture model, no inducible missegregation, and ploidy-independent survival. Then test direct oxygen-linked missegregation, selection on pre-existing states, and shared versus context-specific response functions.
2. **Small symbolic identifiability study:** Use a five-live-class branching model with retained dead biomass and simplified rational response functions to expose parameter confounding with COMBOS, or DAISY if needed. Start with controlled culture oxygen and introduce latent tumor resources afterward. Examine the resulting ambiguities under the actual sampling schedule in the full model.

One small model is not expected to accomplish both jobs. Structural conclusions from the symbolic reduction and empirical support from refitted full-grid variants must remain separate.

The central estimand is the **production and subsequent expansion of viable chromosome-loss descendants**, compared with chromosome-gain and WGD-derived descendants. Generation probability and survival probability may be poorly separated while their combined output is constrained. This is a hypothesis to test, not an assumption to encode into the analysis.

The final report must answer:

- Can each restricted model reproduce growth and the observed chromosome distributions, including both starting-ploidy backgrounds and controls?
- Which individual parameters, function values, or combinations remain ambiguous?
- Are viable-daughter flux predictions more stable than the underlying parameters?
- Can different resource histories explain culture-tumor contrasts with shared cellular response functions, or does releasing particular response blocks improve the explanation?
- Which measurement streams distinguish death production from clearance, or variation generation from survival filtering?
- Which additional measurement would distinguish observationally similar explanations?

## 2. Existing Implementation and Reuse

These are inspected local reference paths, not a claim that their defaults reproduce the latest production analysis. First reconcile them with the version used for the manuscript.

| Component | Existing location | Reuse or change |
| --- | --- | --- |
| Input loading and configuration | `oxygen/code/in-vitro-utils/io.R` | Reuse fit-object and flow loading; require explicit paths and audit all sample matches. |
| Lineage construction | `oxygen/code/in-vitro-utils/lineage_adapter.R` and `objective.R` | Preserve parent-child links, assigned oxygen, passage inputs, and physical observation identifiers. |
| Culture propagation and reseeding | `oxygen/code/in-vitro-utils/runner.R` | Reuse the production passage rule; do not substitute the long-horizon Figure 6 passage protocol. |
| Observation models and optimizer mapping | `oxygen/code/in-vitro-utils/objective.R` | Preserve scoring; generalize the parameter specification to distinguish active, fixed, and inactive parameters. |
| Calibration backend | `oxygen/code/O2_supply_demand_MAP/util/o2_supply_demand_map_fit_invitro_backend.R` | Reuse optimizer execution, logging, and provenance rather than creating a second fitter. |
| Death, division, and daughter kernel | `oxygen/code/O2_supply_demand_MAP/model/model_O2_supply_demand_MAP.cpp` | Reuse unchanged where possible; the first two restrictions are already representable, while an oxygen-only missegregation mode needs an explicit tested dispatch. |
| Parameter normalization | `oxygen/code/O2_supply_demand_MAP/util/o2_supply_demand_map_common_semantics.R` | Check exact zeros survive normalization and worker serialization. |
| Kernel and event-accounting tests | `oxygen/tests/testthat/test-buffer-missegregation.R`, `test-event-accounting.R`, `test-invitro-defaults.R` | Extend rather than duplicate the existing oracles. |
| Existing profile runner | `oxygen/code/O2_supply_demand_MAP/optimizer/profile_likelihood_O2_supply_demand_MAP.R` | Reference only: currently invokes in vivo fitting and must not be used unchanged for culture profiles. |
| Tumor calibration and observation assembly | `oxygen/code/O2_supply_demand_MAP/util/o2_supply_demand_map_fit_invivo_backend.R` | Reuse the frozen production burden, terminal chromosome-distribution, and necrosis operators in the tumor stage. |
| Cross-context parameter mapping | `oxygen/code/O2_supply_demand_MAP/util/o2_supply_demand_map_fit_joint_backend.R` | Reuse context mappings and joint data assembly; add exact sharing and explicit release blocks rather than treating a soft penalty as equality. |
| Necrosis and joint-contract tests | `oxygen/tests/testthat/test-necrosis-loss.R` and `test-joint-soft-coupling.R` | Extend harvest-only observation, retained-biomass, and exact sharing tests. |

The local culture path consumes `fit_data.Rds`, `jobs_2N.Rds`, `jobs_4N.Rds`, and an attached G0/G1 density table. Candidate locations are `oxygen/ploidyOxygen/data/fit_objects/` and `oxygen/data/g0g1_ploidy_density_grid.csv`. Treat these as discovery locations, not hard-coded future inputs. Regenerate observation counts from the approved production snapshot instead of copying numbers from older reports.

## 3. Freeze the Experimental and Computational Contract

### WP0: Production Replay and Observation Audit

**Inputs:** Approved production source/container, parameter table, fitted endpoints, culture fit objects, flow densities, lineage metadata, and fitting command/configuration. Before the tumor stage, also require burden records, harvest-linked chromosome data, histology mapping, tumor fit configurations, and joint parameter-sharing/penalty definitions.

1. Record source commit and any local patch, container digest, R/package/compiler versions, input hashes, parameter transforms/bounds, objective weights, numerical settings, and seed identifiers. Obtain missing production inputs before calibration; do not silently substitute an older top fit.
2. Export an observation manifest with physical sample ID, modality, cohort, lineage/parent segment, passage, oxygen, interval/time, observed quantity, sample size where available, and source path/hash. Retain shared ancestry between nominal replicate branches.
3. Separate measured quantities from conditioning inputs and fitted quantities. Passage growth is interval information, not a continuously observed live-cell trajectory. Flow density grid points are not independent cells. Initial distributions and observation-error parameters may be fitted rather than known exactly.
4. Verify the actual passage-selection rule. The inspected runner selects a simulated day closest to a target final live-cell count and passes its composition forward. This is not necessarily propagation to the recorded passage duration. Record selected days and target-count use; preserve candidate days, tie-breaking, fallbacks, state propagation, and selected-duration growth scoring for all three refits. Any change to that rule requires a separate, labeled sensitivity analysis.
5. Replay at least one accepted production endpoint without optimization. Require agreement in component scores, scored observation IDs, predicted distributions, selected passage days, and growth predictions within numerical tolerances established before the comparison. For tumors, verify actual harvest/assay dates against any last-burden-time fallback in preprocessing; do not silently equate them. A date correction requires a separately documented baseline update.
6. Confirm which observations are actually scored. The inspected objective uses growth, chromosome-count samples, and flow densities; do not add a hypothetical death/CIN time series merely because a report displays a modeled one.

**Outputs:** `input_manifest.tsv`, `observation_manifest.tsv`, `production_contract.json`, `baseline_replay.tsv`, and `baseline_replay_report.md`.

**Stop condition:** Missing provenance, a nonreproduced baseline, inconsistent observation mapping, or an unresolved passage convention. Resolve before fitting; the symbolic feasibility study may proceed with explicitly provisional assumptions.

### Observation Contract Across Contexts

| Stream | Prediction and actual sampling | Requirement |
| --- | --- | --- |
| Culture growth | Growth over recorded passage intervals, using the verified production selection/scoring rule | Preserve assigned oxygen histories, lineage ancestry, and reseeding; do not treat daily simulated counts as observations. |
| Culture chromosome counts and flow | Full distributions at their measured passages | Preserve separate modality transformations and sample identities; no replacement by mean ploidy. |
| Tumor burden | Total burden at each measured time | Retain the production mapping from live and retained dead material to volume, inclusion filters, and scale parameters. |
| Tumor chromosomes | Viable chromosome-state distribution at the matching harvest only | Keep the harvest-specific endpoint and cohort; no invented longitudinal tumor karyotypes. |
| Tumor necrosis | Retained dead volume divided by total volume at the matching harvest only | Retain clearance and volume weighting; necrotic fraction is neither instantaneous death rate nor cumulative death count. |

Verify these contracts against the frozen code, not just report descriptions. Apply the same schedule to numerical checks of the five-class reduction, using declared projection/reconstruction operators for distributional observations and recording their approximation error. Do not refit a surrogate to a mean-only version of the dataset and describe it as a test of the full observation design.

## 4. Question-Driven Full-Grid Model Family

### First Batch: Three Culture Models

For chromosome state \(N\) and assigned oxygen \(O\), the inspected core uses

\[
p_N(O)=p_0+a\frac{\mu_N(O)}{\mu_N(O)+k},
\qquad
s_N=s_{\max}\exp\!\left[-\beta\left(\frac{44}{N}\right)^b\right],
\qquad S_N(m)=s_N^{|m|}.
\]

Here \(p_0\) is baseline per-chromosome missegregation probability; \(a\) is the maximal death-hazard-linked increment; \(k\) is a death-hazard half-saturation scale; and \(S_N(m)\) is survival after an error involving \(|m|\) chromosome copies. The value 44 is `2 * N_UNIT` in the current configuration. Retain production probability clipping and numerical boundary behavior in refits and record when they are active.

| Model | Exact restriction | Parameters retained or removed | Biological question |
| --- | --- | --- | --- |
| `reference` | None | Refit all production culture parameters. | What can the proposed architecture explain? |
| `no_inducible_ms` | `p_misseg = 0`, hence \(p_N=p_0\) | Retain fitted `p_mis_base`, WGD, growth/death functions, and buffering. Remove `p_misseg` from the active vector and mark `k_o_mis` inactive. | Do the observations require stress to increase per-chromosome error probability, beyond baseline variation and selection? |
| `flat_survival` | `buffer_beta = 0`, hence \(S_N(m)=s_{\mathrm{flat}}^{|m|}\) | Refit `buffer_smax = s_flat`; remove `buffer_beta` and mark `buffer_n_exp` inactive. Retain error-size dependence, inducible CIN, WGD, and direct ploidy-dependent growth/death. | Is a parental-ploidy buffering advantage required, beyond survival costs depending on error size? |

Important distinctions:

- `no_inducible_ms` is not a no-CIN model: baseline missegregation remains, and more chromosome copies still create more opportunities for errors.
- `flat_survival` is not an all-daughters-survive model. Setting survival to one would test a different hypothesis. Setting the exponent to zero while fitting both the ceiling and penalty would introduce a redundant parameter combination; use the restriction above instead.
- WGD remains a constant-probability division branch in all three models. Preserve one tracked doubled lineage for WGD versus two potential daughters for an ordinary division, and remove each dividing mother once.
- Do not change direct high-ploidy growth/death costs, initial-state flexibility, or observation models between alternatives. Their ability to compensate is part of the scientific test.

The inspected optimizer fits 20 quantities: 14 mechanistic parameters, four initial-distribution parameters, and two observation-error scales. Under that specification each restricted model has 18 active quantities. Recount from the frozen production configuration and export the active/fixed/inactive table; do not assume these counts apply to every archived version.

### WP1: Exact Restrictions and Tests

The current optimizer hard-codes the 20-parameter vector and positive log bounds. Implement a named model-variant specification with active parameters, exact fixed values, transforms, and inactive reasons. Reconstruct the full runtime parameter list before invoking the existing model. A tiny positive bound is not an exact zero null. Keep production defaults and parameter CSVs unchanged.

Tests must establish:

- The default reference round-trips between natural and optimizer scales and reproduces the replay.
- Zero inducible amplitude gives the same \(p_0\) at every oxygen/state and is invariant to the inactive \(k\).
- Zero buffering penalty gives identical survival across parental states for the same \(|m|\), while still penalizing larger errors and preserving survival one at \(m=0\).
- Fixed zero values survive normalization, saving/loading, and parallel worker execution; inactive parameters are absent from optimization and uncertainty summaries.
- Ordinary division, WGD, maternal removal, biological nonviability, and out-of-grid loss retain correct accounting.
- All models score the identical predefined observation set. Failed/nonfinite predictions incur a declared failure penalty, not silent removal from modality averages.

Exact nulls may lie outside the reference fit's original positive parameter bounds. Keep the production-bound reference as the primary replay/comparator and explicitly report this fact. If using a nesting argument, additionally evaluate the boundary-extended reference at the null solution and reoptimize it with exact zeros admitted; report that as a separate comparison. Do not claim numerical nesting merely because the equations have a boundary limit.

**Outputs:** Variant specifications, parameter-status table, regression tests, and a validation report. No full fitting batch until these pass.

### Second Batch: Two Additional Mechanism Alternatives

Both are required planned comparisons, with separate pilot/compute approval after the first batch. Retain the full chromosome grid, assigned culture exposures, passage rules, and observation models; refit all remaining parameters.

| Model | Definition | Question and safeguards |
| --- | --- | --- |
| `oxygen_linked_ms` | Replace the death-linked branch with \(p_N(O)=p_0+aK_{\mathrm{MS}}^{n_O}/(K_{\mathrm{MS}}^{n_O}+O^{n_O})\), independent of parental state at a given \(O\). | Can the observations distinguish experienced-stress feedback from a direct environmental response? Refit the baseline, amplitude, and oxygen half-response scale; retain the production oxygen-response exponent, buffering, WGD, and direct growth/death. Mark `k_o_mis` inactive. The new scale has oxygen units. This is a nonnested alternative, not a parameter-zero test. |
| `preexisting_selection` | Set baseline missegregation, inducible missegregation, and WGD probability to zero; the ordinary daughter kernel is \(2I\). | Can selection among initially present states explain remodeling without generating new states? Retain state-dependent division/death and refit initial minority populations under observational constraints. Remove inactive CIN and survival parameters. |

For `oxygen_linked_ms`, constrain probabilities and prespecify comparable flexibility and bounds. At constant oxygen, population-average per-chromosome missegregation must be invariant to composition in this alternative; the number of errors per division or per day can still change with chromosome content and division frequency. Compare at the actual exposure history and inspect constant-oxygen intervals for informative composition changes. Do not infer feedback merely from changes in an error-count rate.

For `preexisting_selection`:

- Use a small prespecified mixture on the full grid, for example \((1-f)g_{\mathrm{major}}+fg_{\mathrm{minor}}\), with minority lower-ploidy or WGD-compatible states as appropriate to the starting cohort. Add an intermediate component only if justified; do not freely estimate one initial weight per chromosome state. Bound weights and component shapes using initial karyotype/flow observations, sample sizes, and assay detection/error assumptions. Absence in a sparse initial sample is not proof of zero prevalence.
- Do not let later observations define the admissible initial support or its bounds after seeing which variant fits better. Refit mixture weights within those bounds, including an explicit sensitivity to uncertain detection limits.
- Give comparison models the same initial-mixture flexibility in a matched secondary comparison, so failure or success is not driven simply by changing the initial-state model. Keep the production-initialization replay separately labeled.
- Check physical founding plausibility if the explanation depends on fewer than one expected minority cell in an inoculum or passage. Evaluate founder sampling probabilities rather than treating arbitrarily tiny deterministic tails as guaranteed founders.
- Test that an exactly absent state remains absent in latent propagation; observation smoothing must not seed biological states. Composition-preserving reseeding cannot create new states, and nonviability parameters must not remain falsely active when errors are disabled. Ensure fitted mixtures reach both interactive and cached/batch model initializers.

Export `extended_variant_definitions.tsv`, `initial_mixture_constraints.tsv`, and the same fit/prediction tables as for the first batch.

### Tumor and Cross-Context Stage: Shared Responses Versus Different Resource Histories

First replay the approved separate tumor and joint fits under WP0, including retained dead biomass and harvest-only measurements. Test the selected mechanism alternatives in tumors on the full grid; culture oxygen stays an assigned input, whereas tumor oxygen is latent and inferred through the approved burden/resource dynamics.

Then fit the following joint models, with an explicit shared/free map:

| Joint model | Cellular parameters allowed to differ between culture and tumor |
| --- | --- |
| `shared_response` | None of the cellular response parameters. Allow context-specific resource histories and appropriate observation/initial-state nuisance parameters. |
| `release_growth_death` | `lam_max`, `alpha_o2`, `gamma_growth`, `mu_hp`, and `gamma_mu`, subject to the frozen production definitions. |
| `release_cin` | `p_mis_base`, `p_misseg`, and `k_o_mis`. Growth/death and survival remain shared. |
| `release_survival` | `buffer_smax`, `buffer_beta`, and `buffer_n_exp`. Growth/death and CIN remain shared. |
| `release_combined` | Combinations of the preceding blocks, followed by all three released as a benchmark. Start with single-block releases to expose substitution between explanations. |

Keep `O2_crit` and `n_O` shared initially because they affect multiple responses; release them only in a separately labeled response-shape sensitivity. Keep the per-division WGD probability shared and constant in this planned block comparison; a context-specific constant WGD probability is another declared follow-up if warranted. Do not introduce oxygen-dependent WGD generation. Intersect release blocks with each variant's active parameters. Export every other production parameter with an explicit assignment rather than silently leaving it context-specific.

Require common response-function semantics before declaring parameters shared. For example, reconcile `ploidy_O2_death` modes across contexts: identical numbers under different death-function switches are not identical cellular responses. Preserve distinct experimental protocols, initial states, geometry, clearance, resource histories, and measurement errors where appropriate; they are not all cellular-response parameters that must be equated.

Exact sharing requires a single common parameter or a delta fixed exactly to zero. Strong soft coupling is not the shared-response null; a released block requires genuinely independent unpenalized context copies in the primary release test, not merely a larger coupling scale. Hold the scoring convention and comparable physical bounds fixed, and report data-fit contributions separately from priors/regularization. Where the original contexts use incompatible bounds, obtain an approved common comparison domain before fitting. Compare exact-sharing/release models using the same non-coupling penalty policy; retain the production soft-coupled replay as a separate reference, not as an automatically nested model.

Audit score labels against their implementation: the current joint `objective_unpenalized` still incorporates the in vivo component, which can contain priors. Export data terms, priors, coupling, and phenotype constraints separately. Preserve production constraints for replay, but obtain approval to remove any conclusion-enforcing phenotype gate from the actual mechanism test; otherwise the comparison could assume the outcome it is meant to test.

Reoptimize latent resource parameters and all permitted nuisance parameters for each comparison. Track compensation between resource calibration, `O2_crit`/`n_O`, death, and CIN. Do not interpret identical values of the model's oxygen coordinate as experimentally equivalent culture and tumor environments. If a block release improves fit, ask whether the improvement remains under reasonable latent-resource constraints; otherwise the ecological-versus-cellular distinction may itself be unresolved. Initialize from multiple available fit families, not only the manuscript's preferred family.

**Outputs:** `context_parameter_map.tsv`, `context_release_comparison.tsv`, per-stream predictions, resource-response profiles, and a conclusion distinguishing resource-history sufficiency from evidence for context-specific response functions.

Before either extended stage, add tests for oxygen-only probability invariance, selection-only support preservation, minority-initial-state propagation, exact sharing/release round trips, and stream masks that leave cohort/timing unchanged. New mechanism modes must reach R diagnostics, C++ single-scenario and batch paths, cache signatures, and parallel workers consistently.

## 5. Empirical Model Comparison

### WP2: Fair Refits on the Full Chromosome Grid

**Dependency:** WP0 and WP1 approved and validated.

1. Fit both starting-ploidy backgrounds and all approved control/deprived branches jointly, as in the production culture calibration. Retain the full chromosome distributions and original observation transformations; do not reduce them to mean ploidy.
2. Reoptimize every active mechanism, initial-distribution, and error parameter. Parameter deletion at the reference optimum is only a diagnostic, not the null-model test.
3. Use projected production endpoints plus independent starts within the approved bounds. Benchmark a small pilot, for example five starts per variant, then agree a full search budget based on convergence. Keep search effort comparable by reporting evaluations, wall time, dimension, and attained score, not just equal seed counts. Do not assume five starts establish the best fit.
4. Preserve the production objective first. The inspected implementation averages Gaussian growth scores, within-sample chromosome log scores, and flow-density cross-entropies within modalities before weighting them. It is not automatically an ordinary independent-observation likelihood.
5. Export every component and physical-sample residual alongside the aggregate score. Compare interval growth, complete chromosome distributions, high-chromosome expansion and subsequent reduction in the 2N lineage, downward remodeling in 4N, and control stability. Define chromosome windows from the production convention or prespecify them before looking at model rankings; a high-chromosome window is not direct WGD lineage evidence.
6. Check numerical convergence, boundary occupancy, error-scale inflation, and compensation through direct growth/death costs. Expand numerical grids only as a separately declared check of truncation; do not silently change the biological model.

**Evidence standard:** Report optimized score differences and distribution/growth discrepancies first. Do not attach naive likelihood-ratio p-values, AIC values, or standard profile-likelihood confidence thresholds to the current averaged objective. The inactive \(k\) under the zero-amplitude null and boundary restrictions also invalidate a routine Wilks-theorem argument. Formal inference would require an approved generative observation model and appropriate calibration, not relabeling the existing score. [Composite-likelihood review](https://www3.stat.sinica.edu.tw/sstest/oldpdf/A21n11.pdf)

For validation, prefer withholding an entire physical sampling event or an independent terminal branch where feasible, keeping linked flow/karyotype measurements together. Document common ancestors and whether counts used to calculate a held-out growth observation also condition state selection. Such observations cannot simultaneously serve as supposedly unseen validation data without redesigning the conditioning rule. If no leakage-free validation split is possible, call the comparison in-sample and state that limitation.

Do not bootstrap flow grid points, optimizer endpoints, or individual karyotyped cells as independent biological lineages. A later simulation-based recovery study may generate data from each variant at the real schedule and refit all variants, but must first specify realistic sample-level noise and dependence. It is a gated extension, not part of the initial fitting batch.

**Outputs:** `fit_index.tsv`, `fit_components.tsv`, `observation_predictions.tsv`, `variant_comparison.tsv`, and `optimizer_convergence.md`, plus complete configurations and logs.

## 6. What to Estimate Besides Parameters

Use state-resolved trajectories to connect fit equivalence to the paper's mechanism. For live count \(x_N(t)\), division rate \(\lambda_N(t)\), WGD probability \(w\), and ordinary-division viable-daughter kernel \(B_{N',N}(t)\), export

\[
J_{N'\leftarrow N}^{\mathrm{ordinary}}(t)
=x_N(t)\lambda_N(t)(1-w)B_{N',N}(t).
\]

The existing kernel already counts expected daughters; do not multiply it by two again. Export the WGD flux separately with its one-lineage weight. Sum \(J\) over \(N'<N\) for viable chromosome-loss output and over \(N'>N\) for gain output. Report cells/day and integrated descendants produced over each modeled passage, with an additional division-normalized yield where useful. These are production fluxes, not counts of unique eventual descendants.

Export for every competitive fit:

- \(p_N(O)\), division/death rates, and survival for selected error sizes at matched states, including \(N=44,88\).
- Generated versus surviving loss/gain daughters, WGD output, and direct death, separating biological survival loss from numerical/out-of-grid removal.
- Cumulative viable-loss flux and the subsequent abundance/growth of the receiving chromosome states. Do not equate production with expansion or claim ancestry from state counts alone; explicit source tracking is needed for ancestry claims.
- Population-average per-chromosome missegregation and division-weighted error summaries, with denominators and units stated. They are different quantities.
- The ratio of survival near tetraploidy to survival near diploidy for a fixed error size, and function values over the observed oxygen/state domain.

Because the ordinary kernel generates symmetric \(N-M\) and \(N+M\) daughters and applies symmetric survival at fixed parent/error size, loss flux alone does not establish directional ploidy reduction. Check gain flux, differential expansion/death, WGD, and boundary effects as well. Confirm that chromosome-number displacement balances before selection and truncation in the symmetric ordinary-division test.

**Outputs:** `state_function_values.tsv`, `daughter_fluxes.tsv`, `passage_flux_summaries.tsv`, and `flux_accounting_validation.md`. Reuse existing kernel tests and any already validated flux helper; do not depend on completion of a separate Figure 6 analysis plan.

## 7. Small Symbolic Model for COMBOS/DAISY

### WP3: Five Live Classes with Simplified Response Functions

The full model includes a large chromosome grid, nonrational parameterizations, finite sampling, clipping, and passage rules. Do not submit it unchanged and assume the tool solves that complete problem. COMBOS is designed to identify parameters/combinations of rational ODE input-output models; inputs, measured outputs, and known versus unknown initial conditions are part of the specification. [COMBOS paper](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0110261)

#### State Classes and Branching

Use five recurrent live-cell classes: **below near-diploid, near-diploid, intermediate, near-tetraploid, and above near-tetraploid**. A concrete pilot uses chromosome intervals \(1\!:\!32,33\!:\!54,55\!:\!76,77\!:\!98,\ge99\), with fixed representative counts \(n_j=(22,44,66,88,110)\). These boundaries and representatives are approximation choices, not biological thresholds or inferred parameters. Freeze them before analysis and test alternative boundaries/lifting weights.

Let \(P\) aggregate chromosome states into classes and \(R\) lift class counts to a fixed within-class distribution, with \(PR=I_5\). Start with mass at the representative count; a fixed weighted lifting is a prespecified refinement. Neither bin edges nor within-bin weights are free fit parameters in the symbolic pilot.

The live equations are

\[
\dot x_i=\sum_j b_j\mathcal B_{ij}x_j-(b_i+d_i)x_i,
\qquad \mathcal B=(1-w)B^{\mathrm{ordinary}}+wW.
\]

Here \(b_j,d_j\) are division and direct-death rates, \(B^{\mathrm{ordinary}}\) counts expected surviving ordinary daughters, and \(W\) counts one retained doubled lineage. A dividing mother is removed once. Keep unchanged division and altered-daughter contributions separately labeled even when both stay in the same coarse class.

For each representative parent \(N=n_j\), build daughters on the chromosome scale before aggregation: \(M\sim\mathrm{Binomial}(N,p_j)\) gives the correlated pair \(N-M,N+M\); apply the same maternal-state survival \(S_N(M)\) to both. For \(M=0\), count two unchanged daughters. The WGD branch maps \(N\) to one \(2N\) lineage with probability \(w\), constant per division. Do not normalize viable daughter-kernel columns to sum to one or multiply already counted daughters by two.

All positive destinations, including upper-tail WGD products, need an explicit class destination. A separately declared biological viability rule may exclude states such as zero chromosomes. Computational truncation, tail pruning, or an omitted bin must not be interpreted as biological death. Track these separately, expand the calculation domain or provide validated tail routing, and fail the reduction check if material overflow remains unexplained. The full-model production replay retains its original boundary convention; this audit does not silently alter it.

#### Simplify Shapes Before Removing Dependencies

For the pilot, define \(z_j=n_j/44\), assigned oxygen \(O\), and

\[
h(O)=\frac{K_O}{K_O+O},\qquad
b_j(O)=\frac{b_{\max}}{1+\alpha h(O)z_j},\qquad
d_j(O)=d_Hh(O)z_j,\qquad
p_j(O)=p_0+a\frac{d_j(O)}{k+d_j(O)}.
\]

This fixes the resource Hill and ploidy exponents at one initially. It preserves resource-dependent growth/death, a state-dependent experienced death hazard, and death-linked inducible missegregation. Changing composition can change population-average \(p_j\) at unchanged oxygen. Require positive denominators and \(p_0\ge0,a\ge0,p_0+a\le1\), so parameter-dependent probability clipping is unnecessary in the symbolic model. \(K_O\) has oxygen units, \(b_{\max},d_H,k\) have day\(^{-1}\) units, and \(\alpha,p_0,a\) are dimensionless. These simplified functions are not replacements for production functions in full-grid refits.

Parameterize buffering by two per-copy survival anchors, \(s_2\) at 44 chromosomes and \(s_4\) at 88. Use

\[
s_N=(1-\theta_N)s_2+\theta_Ns_4,\qquad S_N(m)=s_N^{|m|},
\]

where the fixed interpolation weight is zero at/below 44, \((N-44)/44\) between 44 and 88, and one at/above 88. Thus \(S_N(0)=1\). The flat tails are explicit approximation assumptions, not evidence that biological survival saturates at those counts. Report both the near-diploid survival anchor and the increment \(s_4-s_2\), targeting the paper's survival-floor versus ploidy-dependent-increase distinction. Declare whether \(s_4\ge s_2\) is imposed by the buffering hypothesis; if imposed, do not present the sign as independently learned. Retain event-size dependence and maternal-state symmetry.

The symbolic reference retains all mechanisms; the initial nulls are \(a=0\) with \(k\) inactive, and \(s_2=s_4\). Keep WGD constant. Fixed integer chromosome counts, fixed interpolation weights, and finite binomial sums give rational/polynomial expressions, although high polynomial degree can still make elimination expensive. Do not freely fit Hill exponents, replace the death link with a direct oxygen link, or silently drop large-error events just to make a symbolic result easier to obtain.

#### Retained Dead Biomass and Tumor Resources

Add retained direct-death and failed-daughter biomass stocks, initially \(Q_H,Q_E\):

\[
\dot Q_H=\sum_j v_j d_jx_j-c_HQ_H,\qquad
\dot Q_E=(1-w)\sum_j e_jb_jx_j-c_EQ_E.
\]

The fixed conversion weights \(v_j\) map live cells to biomass/volume; \(e_j\) is the expected retained biomass of biologically nonviable ordinary daughters per division, derived from the same error/survival kernel. Reconcile these weights and retained fractions with the frozen production observation convention rather than equating a missing daughter with a full parent-cell volume. Numerical overflow contributes to neither stock. If explicit biological WGD failure is part of the approved model, account for its failed biomass separately; do not invent such a branch from an upper numerical boundary.

For a declared volume convention, total burden is \(V=\sum_jv_jx_j+Q_H+Q_E\) and necrotic fraction is \((Q_H+Q_E)/V\). Aggregate the two stocks only if their clearance laws, observation weights, and resource-feedback roles permit it; common linear clearance and observation solely through their sum are a sufficient case. Otherwise retain separate stocks, or class-resolved dead states if clearance depends on dead-state composition. Equal aggregate necrosis does not reveal its separate causes. Fixed representative volumes need their own approximation check because coarse projection need not preserve biomass or chromosome moments.

Begin with prescribed culture oxygen; dead stocks remain unobserved unless a measured culture assay actually observes them, and culture passage resets follow the verified protocol. Introduce a latent tumor oxygen/resource state only in the next stage. One rational burden-linked surrogate is

\[
\tau_O\dot O=O_{\min}+\frac{O_S-O_{\min}}{1+\eta V/V_*}-O.
\]

Here \(O_S\ge O_{\min}\ge0\), \(\eta\ge0\), \(\tau_O>0\), and \(V_*\) is a fixed reference volume. Declare initial oxygen known or unknown; hold one resource-scale convention fixed to avoid an artificial scaling symmetry. Use the frozen production burden definition in \(V\); including retained dead volume in a burden-based supply proxy does not mean it consumes oxygen. This rational surrogate changes the production resource law and requires validation. Its identifiability result is not a proof for the original latent-oxygen dynamics, which remain unchanged in the full-grid context comparisons.

#### Accounting Benchmark and Reduction Checks

Retain the one-generation ledger below as a unit-level algebraic check, not as the primary reduced model. For one parent and fixed error size \(0<m<N\), let \(q\) be an ordinary-division error-event probability, \(s\) symmetric conditional daughter survival, and \(w\) WGD probability:

| Outcome | Expected count per dividing mother |
| --- | --- |
| Mother removed | \(1\) |
| Unchanged \(N\) daughters | \(u=2(1-w)(1-q)\) |
| Retained \(N-m\) daughters | \(\ell=(1-w)qs\) |
| Retained \(N+m\) daughters | \(\ell=(1-w)qs\) |
| Doubled \(2N\) lineage | \(w\) |
| Biologically nonviable daughters | \(v=2(1-w)q(1-s)\) |

Check \(u+2\ell+w+v=2-w\). Every destination is included; there is no numerical or boundary loss in this diagnostic ledger. Here \(q\) is an event probability per ordinary division, **not** the full model's per-chromosome \(p_N\); \(s\) is conditional event survival, not an original buffering parameter.

First audit this algebraic coefficient map. For example, if \(u,w\) were independently known and \(w<1\), then \(q=1-u/[2(1-w)]\). Also \(2\ell+v=2(1-w)q\). These identities show why unchanged daughters or nonviability can provide information beyond \(qs\); they do not establish that existing culture observations recover those quantities.

For the five-class model, construct the frozen approximation \(A_c=PAR\) using chromosome-scale daughter accounting. With representative lifting this gives the branching equation above. If a weighted lifting uses distinct fine-state division rates within one bin, form rate-weighted daughter fluxes before projection; \(PBR\) combined with an independently averaged division rate is generally not the same operator. Demonstrate exact closure before using that term; it is not expected here.

Test no-error propagation, symmetric loss/gain accounting before biological filtering, WGD multiplicity, within-bin errors, intermediate-state occupancy, and lower/upper-tail routing. Compare projected full-grid trajectories with reduced trajectories for both starting cohorts and multiple resource conditions, including initial high-state expansion and subsequent reduction. Perturb bin edges/representatives. Fixed lifting may erase gradual within-bin chromosome loss or repeated upper-tail WGD; if this changes the mechanism or coexistence pattern, refine the closure or add needed states/moments before using it. Do not manufacture downward transitions to improve the match, and do not replace the five classes with a two-state 2N/4N model.

### Output and Initial-Condition Ladder

| Case | Outputs/assumptions | Purpose |
| --- | --- | --- |
| S0 | Ideal continuously observed total live abundance in culture, or total burden in the tumor extension | Establish what growth/burden alone can distinguish; do not conflate these two outputs. |
| S1 | S0 plus all five viable class fractions, with four independent fractions in the calculation | Test the additional information from distributions; also test fractions without total abundance. |
| S2 | S1 plus total necrotic fraction as a hypothetical continuous output, then direct division, death, error-event, or daughter-survival measurements one at a time | Test death/clearance and generation/survival ambiguity. A continuous necrosis output is a richer hypothetical experiment than terminal histology. |
| S3 | Repeat with known versus unknown initial mixtures/abundance, clearance, volume scales, and (in tumors) latent-resource parameters; compare one versus both starting ploidies | Expose conclusions dependent on external knowledge rather than measured trajectories. |

For every case state whether the input is one constant oxygen level, several separately constant experiments sharing parameters, or a switching schedule. A proof that assumes an arbitrarily varying input must not be presented as a result for the observed exposure history. Test purely parental versus unknown mixed initial states explicitly; generic nonzero initial-state results may not apply to the proposed experiment. Any symbolic passage extension uses externally specified times and declared resets, not target-count matching. Its conclusions apply to that observation scheme, not automatically to the implemented culture objective.

Retained dead stocks are model states, not automatically separate measured outputs. Label continuous outputs above as idealized, noise-free experiments. The real culture interval/distribution and tumor measured-time/harvest-only operators are evaluated separately in WP4, including clearance and unknown initial states. Numerical sensitivity/rank checks at those finite schedules address local distinguishability, not global structural identifiability. [Observation assumptions and identifiability methods](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0027755)

### Tool and Reduction Validation

1. Check available COMBOS/DAISY installation or service and reproduce a published identifiable and nonidentifiable example. Archive tool/version, exact symbolic input, raw output, assumptions, and runtime. Only symbolic equations need be supplied to a service, not experimental records. [DAISY paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC2888537/)
2. Classify results as globally identifiable, locally identifiable, nonidentifiable with explicit combinations, or unresolved/tool-limited. A timeout is not nonidentifiability.
3. Verify any proposed equivalence by substitution into the input-output coefficients; use numerical trajectory agreement only as an additional check. Record positivity/probability restrictions and special boundary cases separately from generic results.
4. Run the five-class rational model first with controlled oxygen, then its latent-resource extension. Use the ledger or additional fixed-parameter cases to diagnose solver complexity, not to replace the main question. Archive unresolved cases, bin/closure checks, and all approximations if the five-class formulation is intractable.
5. Do not transfer a structural proof from the reduced model to the original simulator. Use its combinations to nominate targets for the full-model checks below. Known fixed shape constants in a reduced model must be listed as additional assumptions, not claimed as identifiable original parameters.

**Outputs:** `symbolic_models/`, `five_class_definition.tsv`, `reduction_validation.tsv`, `dead_biomass_contract.md`, `symbolic_case_manifest.tsv`, `identifiability_results.tsv`, `identifiable_combinations.md`, and `reduction_scope.md`.

## 8. Bridge to the Actual Sparse Observation Schedule

### WP4: Local Distinguishability and Focused Profiles

**Dependencies:** WP0, WP1, and symbolic candidates from WP3; use multiple converged full-model fits from WP2 when available.

1. Define the finite observation map by composing the production passage propagators, reset/selection rules, and growth/karyotype/flow operators at the actual observation events. In the tumor stage, include burden only at measured times and viable chromosome distributions and necrotic fractions only at their matched harvests. Include fitted initial conditions and observation-error parameters as nuisance quantities. Do not replace sparse karyotypes or terminal necrosis with continuous observation.
2. Evaluate scaled sensitivities and singular values of this observation map near several competing fits, with perturbation-size and numerical-tolerance checks. Remove normalization redundancies. Include observation-scale effects through the declared score/information calculation rather than pretending growth noise scales are parameters of the noiseless mean trajectory.
3. Identify weak directions involving inducible amplitude/half-saturation, survival shape, direct growth/death, WGD, and initial composition; then add dead-biomass clearance and latent resource parameters in tumors. Compare with symbolic combinations without assuming exact correspondence between reduced function anchors and original parameters.
4. Audit changes in selected passage day under perturbation. The nearest-target-count selection can create nonsmooth objectives; a smooth Jacobian within one selection region does not describe all alternatives. Use profiles and explicit state-selection logs, and stop any derivative-based claim that is unstable to this behavior.
5. Run a small number of targeted **profile-objective** scans, reoptimizing all other active parameters with multiple starts. Prioritize inducible amplitude and ploidy-dependent survival, then one representative function/flux combination if technically feasible. Direct flux profiling requires a constrained optimization formulation; do not substitute an arbitrary penalty without validation.
6. Show whether similar-quality fits support different generation/survival parameters yet similar viable-loss output and receiving-state expansion. State any objective-tolerance cutoff and test sensitivity to it; it is not a confidence region without calibration.

If a normalized generative likelihood is subsequently approved, extend to likelihood profiles and simulation-calibrated intervals. Until then, use profiles to demonstrate compensation and practical resolution, not formal parameter confidence limits. Local numerical rank is not a global structural-identifiability proof. [Profile-likelihood framework](https://pubmed.ncbi.nlm.nih.gov/19505944/)

### WP4b: Remove One Measurement Stream at a Time

This is a required analysis of information contribution, distinct from removing a biological mechanism. Refit after each stream removal and repeat local-rank/profile and flux-resolution checks. Use the reference mechanism first; extend to ambiguous alternatives only where the contrast is informative.

| Context | Measurement sets | Main question |
| --- | --- | --- |
| Culture | All three streams; omit growth, karyotype, or flow in turn; also growth alone | Do chromosome distributions separate generation from survival beyond what growth tells us? Are karyotype and flow complementary or redundant? |
| Tumor | Burden plus harvest chromosomes plus harvest necrosis; omit each stream in turn | Does necrosis constrain death production versus clearance? Do terminal chromosomes resolve mechanisms that burden alone cannot? |
| Joint | Remove one context-specific stream while holding the other context's measurements fixed | Which measurements support a claimed context difference rather than allowing one context to inherit its constraints from the other? |

Freeze the included lineages/tumors and all retained stream weights before applying score masks. Do not let `paired_only` or missing-stream filters silently change the cohort in each comparison. Mark nuisance parameters that become inactive when their observation stream is removed. Keep removed outcomes out of warm-start selection and state conditioning for any prediction-validation claim; otherwise label the result as a conditional information audit, not a fully withheld-data test.

Compare changes in weak directions, profile shapes, fitted function/flux ranges, and predictions for the omitted stream. Total objectives computed from different measurement sets are not directly comparable. If the production clearance rate is fixed, repeat a prespecified free-clearance sensitivity before claiming that necrosis separates production and removal. Retain clearance as potentially active when necrosis is omitted, because burden still includes retained dead biomass. A terminal necrotic fraction may still leave production and removal confounded; report the result rather than presuming resolution.

**Outputs:** `finite_observation_map.md`, `local_sensitivity_diagnostics.tsv`, `profiles/`, `parameter_vs_flux_resolution.tsv`, `stream_ablation_manifest.tsv`, `stream_information_comparison.tsv`, and `practical_resolution_report.md`.

## 9. Translate Ambiguity into a Measurement Recommendation

### WP5: Small Experiment-Discrimination Study

Use parameter sets or model variants that fit the existing data similarly but imply different mechanisms. Compare their predictions for candidate additional measurements at matched oxygen and starting ploidy. Select observation times where predictions separate within an experimentally accessible interval, not merely at an asymptotic endpoint.

| Unresolved distinction | Candidate additional measurement | What it separates |
| --- | --- | --- |
| Frequent errors with low survival versus rare errors with high survival | Division-resolved imaging that measures segregation errors before daughter loss and follows both daughter fates | Error generation from conditional retention; align event-level and per-chromosome units. |
| Error probability versus division frequency | Error observations with a denominator of scored divisions, plus division timing | Per-division errors from event flux per day. |
| Direct stress death versus loss of altered daughters | Fate tracking that distinguishes death without a scored error from post-error arrest/death, with a defined follow-up window | Direct hazard from survival filtering; avoid treating all missing daughters as deaths. |
| WGD generation versus preferential expansion of pre-existing high-ploidy cells | Direct WGD/cell-division event observations with initial-state and subsequent lineage tracking | Generation from selection and starting heterogeneity. |
| Viable-loss output versus later expansion | Paired short-interval chromosome-state/fate measurements and longitudinal growth | Production of altered descendants from their subsequent fitness. |

Test candidate outputs symbolically where possible, then numerically with plausible measurement noise on the full model. Additional chromosome snapshots may improve trajectory resolution without directly measuring error generation; report that distinction. Do not claim that endpoint transcriptomic stress signatures or multinucleation directly resolve per-chromosome missegregation or daughter survival.

**Outputs:** `candidate_measurements.tsv` and a short `recommended_next_measurement.md` identifying the leading ambiguity, discriminating observable, required timing/denominator, and expected limitations. This plan designs comparisons; it does not authorize new experiments.

## 10. Implementation Layout and Work Order

Proposed new analysis entrypoints/helpers:

```text
oxygen/code/O2_supply_demand_MAP/analysis/structural_identifiability/
  README.md
  analysis_config.example.yml
  audit_inputs.R
  model_variants.R
  run_variant_fits.R
  summarize_variant_fits.R
  export_daughter_fluxes.R
  build_symbolic_cases.R
  validate_coarse_model.R
  symbolic/                       # portable equations and tool input templates
  run_observation_map_checks.R
  run_targeted_profiles.R
  run_stream_ablation.R
  run_context_comparisons.R
  compare_candidate_measurements.R
```

Keep shared parameter conversion/fit changes in the existing culture, tumor, and joint helpers. Entry scripts should call those helpers, not copy the simulator. Add focused tests under `oxygen/tests/testthat/`. No dependency on an uncommitted analysis document or another repository is required.

Each approved run gets an immutable output directory:

```text
oxygen/results/structural_identifiability/<run_id>/
  manifests/
  audit/
  fits/<context>/<mechanism_variant>/<sharing_model>/<measurement_set>/
  symbolic/
  profiles/
  fluxes/
  stream_ablation/
  context_comparison/
  measurement_design/
  report/
```

Configuration must include production source/container, input paths/hashes, context and mechanism variant, exact shared/released parameter map, observation-stream mask, objective/bounds/penalties, passage/harvest convention, observation manifest, numerical tolerances, seeds/search budgets, output root, coarse classes/lifting/biomass rules, and symbolic output/initial-condition cases. Enumerate approved combinations in a run manifest rather than launching the full Cartesian product of variants, releases, and stream masks. Refuse an existing completed run ID unless explicitly resuming it after manifest compatibility checks. Keep bulk trajectories/results out of commits; commit the plan, code, tests, compact manifests, and approved summaries.

Recommended order:

1. Approve this scientific scope and locate the exact production culture snapshot.
2. Complete WP0 replay and WP1 restriction tests; in parallel, define and validate the five-class culture model in WP3 and benchmark the symbolic tool. The ledger is only a unit check.
3. Review a pilot compute estimate and approve the initial three-model culture refits. Perform focused WP4 checks and culture stream-removal comparisons.
4. Approve the second culture batch for oxygen-linked missegregation and pre-existing selection, using matched initial-state assumptions.
5. Freeze/replay the tumor and joint baselines, then approve tumor variants, exact-sharing/block-release comparisons, and tumor/joint stream-removal runs. Extend the symbolic model with retained dead biomass observations and latent resources without inventing longitudinal harvest measurements.
6. Combine empirical and symbolic results in focused profiles and WP5 measurement recommendations. Obtain a separate budget for broad profiles or synthetic refit studies, and produce a compact evidence report before proposing manuscript changes.

Run R through `scripts/agentRrunner.sh` in the pinned environment. Smoke-test one small unit before any batch. Use the approved cluster/container for fitting or work expected to exceed approximately five minutes; check current scheduler limits before submission, cap concurrency, and archive commands/job IDs/resource use. No fitting, large symbolic elimination, or manuscript regeneration is authorized by drafting this plan.

## 11. Acceptance and Interpretation

Deliver a short report with: model/parameter-status table; observed-versus-predicted growth, burden, distributions, and necrosis by applicable variant; symbolic output-to-combination table; parameter-versus-flux resolution summary; stream-information and context-release comparisons; and the recommended additional measurement. Any plots are new diagnostics, not replacements for manuscript panels without approval.

| Possible result | Defensible conclusion |
| --- | --- |
| A restricted model fits equally well after adequate optimization | The observations in the evaluated context do not require the removed feature within the tested family. |
| A restricted model consistently fails specific observed trajectories across well-converged refits | Evidence favors that feature relative to this alternative, subject to observation-model and validation limitations; this is not proof against every alternative mechanism. |
| Generation and survival vary widely but viable-daughter output is stable | Emphasize the constrained effective production/expansion process; do not overinterpret individual parameters. |
| Even viable-loss output differs across similar-quality fits | The route to ploidy reduction remains unresolved by these data; prioritize the discriminating measurement. |
| A reduced model is structurally identifiable but full-model profiles are broad | Ideal distinguishability does not imply practical resolution in the actual experiment. |
| Shared cellular responses fit both contexts with different resource histories | The tested data do not require intrinsic response differences; this is sufficiency within the model family, not proof of identical biology. |
| Releasing a response block improves the cross-context explanation | Evidence favors that block relative to shared-response alternatives only after checking resource, initial-state, penalty, and observation-model compensation. |
| Removing necrosis or a distribution stream widens a relevant ambiguity | That stream contributes information under the tested design; quantify the remaining ambiguity rather than declaring complete identification. |
| Symbolic analysis is intractable or depends on unrealistic observation assumptions | Report a feasibility boundary, not a full-model identifiability verdict. |

Completion requires baseline replay, exact-null and sharing tests, validated coarse-state/biomass accounting, consistent observations across variants, documented optimizer adequacy, explicit symbolic assumptions, stream-information tests, a separation of numerical from biological loss, and conclusions limited to the evidence actually obtained. Report completion separately for culture, tumor, and joint stages; finishing the first culture batch does not complete the full plan. A successful reference fit alone does not establish that its architecture was necessary.

## 12. Coverage of the Reviewed Brief

| Requested element | Plan location |
| --- | --- |
| Two complementary reductions with the selection-versus-generation/survival question | Section 1. |
| Five live classes, retained dead biomass, one-lineage WGD, and honest coarse-graining | WP3: state/branching, retained biomass, and reduction checks. |
| Rational growth/death, death-linked inducible errors, survival anchors plus event-size dependence, constant WGD | WP3: response-function definitions. |
| No inducible errors; flat survival; oxygen-only errors; pre-existing selection | Section 4: first and second culture batches. |
| Defensible minority initial populations and refitting all remaining parameters | Section 4: mixture constraints; WP2. |
| Shared responses across contexts followed by selective releases | Section 4: tumor and cross-context stage. |
| Real culture passage/distribution schedules and measured-time/harvest-only tumor observations | Section 3 observation contract; WP0 and WP4. |
| Ideal continuous SI separated from actual finite-schedule numerical checks | WP3 output ladder; WP4. |
| Leave-one-stream-out tests, especially necrosis/clearance and distributions/generation/survival | WP4b. |
| Measurable combinations and selection of the next discriminating experiment | Sections 6 and 9. |
