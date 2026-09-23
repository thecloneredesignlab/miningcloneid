# Profile Objectives and Calibrated Likelihood Inference

**Status:** Implementation plan only. No fitting, simulation, profile calculation, or statistical calibration has been performed under this plan.
**Date:** 2026-09-23.
**Scope:** Practical identifiability in the existing culture, tumor, and joint models, beginning with a small culture pilot. Preserve production outputs and manuscript assets.
**Companion:** [Structural identifiability and mechanism tests](invitro_mechanism_identifiability_implementation_plan.md), especially WP0, WP1, WP4, and WP4b. This plan specifies the profile workflow and inference gates in more detail; it does not replace the symbolic COMBOS/DAISY study or require that study to finish first. eFAST remains separate.

## 1. Questions and Deliverables

The immediate question is whether the observed growth and chromosome distributions constrain particular mechanisms after other parameters are allowed to compensate. Start with:

1. Does the culture dataset constrain the maximal inducible per-chromosome missegregation increment, including whether zero inducible missegregation can explain the observations?
2. Does it constrain the ploidy dependence of post-missegregation survival, including a survival function that is independent of parental ploidy but still depends on error size?
3. Which measurement streams provide those constraints, and which apparent constraints instead come from bounds, priors, or coupling penalties?
4. Are effective missegregation, survival contrasts, or viable chromosome-loss flux better constrained than individual parameters?

Deliver verified **profile-objective curves** first. A curve describes the best fit attainable at each imposed target value, not the effect of changing one parameter while keeping the others fixed. Profile likelihood is an established approach to practical identifiability, but likelihood-based interval interpretation requires the appropriate statistical objective and assumptions. [Raue et al., 2009](https://pubmed.ncbi.nlm.nih.gov/19505944/).

Formal confidence intervals are a later, separately approved deliverable. The present weighted scores must not acquire likelihood-ratio confidence labels simply because the existing script is called `profile_likelihood`.

## 2. Definition of a Profile Study

Let a study specify a model variant, admissible parameter domain, observation mask \(S\), retained modality weights, and penalty policy \(r\):

\[
Q_{S,r}(\theta)=\sum_{s\in S}w_s Q_s(\theta)+R_r(\theta),
\qquad
\widehat Q_{S,r}=\min_{\theta\in\Theta_{S,r}}Q_{S,r}(\theta).
\]

For target parameter \(\psi\) and all remaining active unknowns \(\eta\), compute

\[
Q^{\mathrm{prof}}_{S,r}(\psi)
=\min_{\eta:(\psi,\eta)\in\Theta_{S,r}}Q_{S,r}(\psi,\eta),
\qquad
\Delta Q^{\mathrm{prof}}_{S,r}(\psi)
=Q^{\mathrm{prof}}_{S,r}(\psi)-\widehat Q_{S,r}.
\]

These are ideal minima. Numerically, report the lowest verified feasible objective attained and document how thoroughly it was searched. Multiple starts do not certify global optimality.

Required rules:

- Every change in data mask, weights, model variant, bounds, or penalty policy defines a new study with its own unrestricted refit. Archived endpoints are starting candidates, not automatically the correct baseline.
- Use the minimum across validated starts as the profile value. Means, ranges, and the number of optimization endpoints describe search behavior, not inferential uncertainty or posterior mass.
- If a constrained fit improves on its unrestricted baseline beyond the replay tolerance, restart the unrestricted search from that solution and update every delta. Never clip negative deltas to zero to hide a stale baseline.
- A data-only profile requires new constrained and unrestricted optimization with the penalties removed. Subtracting a penalty from penalized optima is not equivalent.
- Modality scores evaluated along the combined-objective optimum are **component traces**, not modality-specific profiles. A modality-only profile requires its own refits.
- Bounds, fixed initial-state assumptions, and hard-shared response parameters still constrain a data-only study. State them explicitly; absence of an added penalty is not absence of assumptions.

## 3. Existing Code: Reuse and Required Changes

The following observations concern the inspected local code. Reconcile it with the approved production version before implementation; do not silently mix older local defaults and newer container-generated fits.

| Component | Inspected behavior | Planned change |
| --- | --- | --- |
| `oxygen/code/O2_supply_demand_MAP/optimizer/profile_likelihood_O2_supply_demand_MAP.R` | Per-seed launches, score collection, adaptive stepping, best-fit selection, and plotting already exist. | Reuse execution/logging concepts and tested helpers; add an explicit objective-diagnostic mode with context adapters. Do not copy the entire runner into another independent fitter. |
| `run_single_seed_fit()` in that runner | Hard-codes `--fit_invivo`. | Dispatch explicitly to culture, tumor, or joint fitting; verify the effective model and observation contract in each result. |
| `make_fixed_parameter_table()` and in vivo configuration synchronization | Fixing a coordinate changes `init/lower/upper`; these fields can also determine prior defaults and, for `o2_S0`, a model oxygen cap. | Separate immutable model/prior configuration from target constraints and optimizer bounds. Fixing a coordinate must not silently change the objective or an ancillary model setting. |
| Culture `ivt_optimizer_spec()` and parameter conversion in `oxygen/code/in-vitro-utils/objective.R` | The optimizer uses a hard-coded 20-parameter positive-log representation. | Implement the companion plan's single active/fixed/inactive parameter map; share it between mechanism restrictions and profiles. Do not rely on an `estimate=FALSE` CSV field that the culture vector does not honor. |
| Profile continuation and the culture/in vivo fit backends | A saved best map does not ensure that the neighbor's full nuisance vector reaches the next optimizer; culture DE initialization also needs explicit wiring. | Pass, validate, and save the actual optimizer initial population/vector. Test that baseline and neighboring solutions are really used. |
| `read_single_seed_result()` | A finite objective and reported target can be treated as completion. | Require run fingerprints, target agreement, feasibility, exit/termination status, expected observations, and independent objective reevaluation. Do not accept stale output from an interrupted rerun. |
| Existing profile plots/reducers | Some primary curves use means across seeds; raw modality components are displayed along combined-objective fits. | Plot verified minima for inference; keep search dispersion and component traces separately labeled. No synthetic completed seed reconstructed from a summary fallback. |
| Existing retry/grid logic | Reducing the step and repeating the seed list is not an additional search at the failed coordinate. | Separate fixed-coordinate restarts from adding/moving grid points; give every restart its own reproducible ID and seed. |
| Existing thresholds | Defaults include `1.92`, `3.84`, CI interpolation, and CI-related refinement/stopping. | Disable all of these in objective-diagnostic mode. Use declared grid, numerical, and compute criteria instead. Preserve legacy artifacts without relabeling them. |

Supporting backends are `oxygen/code/O2_supply_demand_MAP/util/o2_supply_demand_map_fit_{invitro,invivo,joint}_backend.R`. Continue to use their simulation and observation machinery. The culture kernel already represents the first two exact restrictions; the main initial engineering work is parameter mapping and reliable optimization orchestration, not a new simulator.

## 4. WP0: Freeze Inputs, Observation Design, and Objective

**Inputs:** Approved source/container, parameter table, fitted endpoints from distinct solution regions, fit objects, flow-density inputs, lineage metadata, and production configuration. Reuse the companion plan's provenance exports when available.

1. Record source commit plus any local patch, image digest, R/package/compiler versions, all input hashes, numerical settings, transforms, bounds, and score definitions. An unavailable current production input is a blocker, not permission to substitute an older report.
2. Enumerate physical observation IDs, lineage/parent links, cohort, passage or time, modality, sample size where known, and assigned oxygen. Retain distributions and actual sampling times. Daily simulated values are not additional observations.
3. Export a machine-readable objective ledger: raw component score, normalization denominator, weight, prior, soft-coupling term, phenotype penalty, and failure penalty. Verify that it sums to the reported objective. Do not assume a field named `objective_unpenalized` excludes a prior nested inside a context score.
4. Replay at least one accepted production endpoint without optimization, comparing predictions, selected passage days, scored observation IDs, component scores, and total score. Establish numerical tolerances before profiling.
5. Freeze the scored observations per study. A nonfinite prediction is a recorded failure under the declared policy, not permission to remove that sample from a modality average.
6. Freeze conditioning inputs separately from observations used in the score. In the inspected culture runner, a final live-cell count selects a simulated passage-end day and therefore the state propagated to descendants. Preserve this rule for the primary study and log every selected day, tie-break, and fallback. Its discrete changes can create genuine objective kinks.
7. For later tumor studies, verify burden times, actual harvest-linked chromosome and necrosis samples, retained-dead-biomass observation rules, latent-resource parameters, and fitted scale/initial-condition parameters. Do not turn harvest-only observations into continuous measurements.

The current culture score averages growth contributions and sample-level chromosome and flow scores before combining modalities. A flow density grid is not a set of independent cells. Keep these conventions for the first diagnostic study, while explicitly describing the result as a profile of that score.

**Outputs:** `production_contract.json`, `input_manifest.tsv`, `observation_manifest.tsv`, `objective_components.tsv`, `baseline_replay.tsv`, and `baseline_replay_report.md`.

**Gate:** No biological profile jobs until the production replay, observation mapping, and immutable objective contract pass review.

## 5. WP1: Targets, Exact Nulls, and Nuisance Parameters

### Initial Culture Targets

Profile the full culture observation model using both starting-ploidy cohorts and their control/deprivation histories. In the inspected notation,

\[
p_N(O)=p_0+a\frac{\mu_N(O)}{k+\mu_N(O)},
\qquad
S_N(m)=\left[s_{\max}\exp\{-\beta(44/N)^b\}\right]^{|m|}.
\]

Retain production clipping, the configured chromosome unit, and the constant-probability, one-cell WGD branch. These equations identify the targets, not a replacement model.

| Target | Biological interpretation | Exact null and nuisance handling |
| --- | --- | --- |
| `p_misseg` (\(a\)) | Maximal death-hazard-linked increment in per-chromosome missegregation probability; dimensionless. | At zero, keep baseline `p_mis_base`, WGD, growth/death, and survival free. `k_o_mis` is data-inactive in this culture model. This is not absence of all CIN. |
| `buffer_beta` (\(\beta\)) | Dimensionless strength of parental-ploidy dependence in survival. | At zero, fit `buffer_smax`; retain error-size dependence. `buffer_n_exp` is data-inactive. This is not survival equal to one. |

Use a **separate, explicitly boundary-extended study** admitting exact zero for both targets in its unrestricted and constrained fits. Preserve the production-bound replay/comparison unchanged. Otherwise a null outside the original positive bounds cannot be interpreted as a nested profile point. Register the extension before any new fits.

For these two coordinates use a zero-capable optimizer mapping, such as a fixed-scale `log1p(theta / scale)`, with scales recorded in configuration. Fixed zeros bypass logarithmic conversion entirely. Other coordinates retain the approved production transforms. A small positive number is not an exact zero test.

Remove an inactive coordinate only when it is inactive in the **complete objective**. In penalized or joint studies, a parameter absent from one context's predictions may still enter a prior, another context, or a coupling term. Retain it or explicitly minimize its residual penalty contribution. Do not drop that contribution or silently change the prior. A prior that excludes zero must be reported as excluding the null, not bypassed with epsilon.

### What Must Remain Adjustable

- All other active mechanistic parameters, including baseline missegregation, the inducible half-saturation scale when active, WGD, growth/death responses, and survival parameters.
- The fitted 2N and 4N initial-distribution parameters and relevant observation-error scales. For the inspected culture version these include `sigma_growth` and `sigma_kary`; verify the effective production specification rather than assuming a universal count.
- In later tumor studies, fitted initial burdens/distributions, scale factors, latent oxygen parameters, clearance, and observation nuisance quantities whenever they affect retained data.
- In joint studies, both context vectors subject to the declared sharing map and penalty policy. Optimizer family labels are starting candidates, not hard constraints for the global profile. Family-conditional diagnostics may be reported separately.

Measured oxygen histories, physical sample identities, fixed assay transformations, and the approved passage rule are not arbitrary nuisance parameters. Changing them defines a separate model/observation sensitivity study.

**Outputs:** `parameter_status.tsv` for every target/null, `profile_domain.json`, and tested conversion between free coordinates and full runtime parameters.

## 6. WP2: Bounded Pilot and Reliable Search

### First Batch

Run one full-data, data-only culture study first, with production modality weights and the declared boundary-extended domain. Replay any production penalties separately; do not silently redefine the original fit. If the production culture objective already has no penalties, avoid a duplicate policy.

- Refit its unrestricted baseline with at least five genuinely distinct starts, including good production endpoints from different solution regions and independent starts. Continue if results reveal unresolved basins.
- For each of the two targets, propose approximately ten values: exact zero, the refitted target value, and positive values spanning a declared small positive probe to the approved upper bound. Log-space the positive values, deduplicate, and fill large remaining gaps. Save the actual values before submission.
- The smallest positive probe is grid resolution, not a new positive lower bound. Choose it from the approved parameter scale; refine the gap to zero if it matters. Do not connect an unresolved large gap as if it had been evaluated.
- At each point, perform five nuisance refits: the unrestricted best vector, a neighboring optimum, a distinct alternative endpoint, and two independent starts. When a neighbor or distinct endpoint is unavailable, replace it with a unique independent start. Save actual starts after projecting onto the fixed-target constraint.
- Search outward from the optimum in both directions. Reserve additional reverse-continuation and independent restarts for nulls, apparent sharp rises, local minima, or basin switches. A point where all runs fail stays unresolved; add restarts at that coordinate before merely moving the grid.

The initial budget is about **100 constrained optimization runs**, plus unrestricted searches, verification restarts, and any approved grid refinement. An optimization run may itself contain global and local stages. The legacy maximum of 20 steps in each direction with 20 starts can require roughly 800 runs per parameter; do not inherit it as a default.

### Acceptance and Refinement

For each candidate require the correct source/data/configuration fingerprint, fixed-target agreement on the natural scale, exact zero when requested, nuisance feasibility, a valid prediction for every retained observation, and independent score reevaluation. Record launch status and optimizer termination separately. A valid checkpoint from a budget-limited run may be retained as a provisional candidate, but must not masquerade as a completed search.

Predeclare `target_tolerance`, `score_replay_tolerance`, and `optimization_agreement_tolerance` from replay and toy-test results. Compare independently initialized or reverse-continuation solutions at biologically important points. Agreement is a numerical check, not a proof of a global minimum.

Refine grids for unresolved curvature, basin changes, selected-passage-day changes, and gaps near zero, within an approved additional-fit budget. Do not impose monotonicity, smooth away kinks, or stop at an inherited confidence cutoff. Boundary-limited and compute-limited curves must be labeled as such.

Every retry needs a new attempt identifier and, when intended as an independent search, a new deterministic seed. Resume only if the complete run fingerprint matches; otherwise allocate a new study/run directory. Never overwrite production runs.

### Execution

Use `scripts/agentRrunner.sh` for R entry points and the frozen environment for workers. Benchmark one representative constrained fit before approving the batch; record walltime, peak memory, and optimizer evaluations. Compute the batch estimate from measured cost and the actual number of fits, not just the number of plotted points.

Use approved cluster arrays for batch or long-running work; inspect current QOS, walltime, and memory limits first. Do not inherit the legacy 62-core setting on a local machine. Parallelize independent starts/points within an explicit total-core limit and avoid nested worker oversubscription. Record scheduler IDs, failures, and checkpoints. Approval to write this plan is not approval to launch fitting jobs.

**Outputs:** Verified unrestricted baselines, per-start results, lower-envelope profiles, nuisance trajectories, observation predictions, search diagnostics, and a runtime report.

**Gate:** Review this two-target pilot before expanding streams, contexts, or targets.

## 7. WP3: Data-Stream and Penalty Comparisons

After the pilot, prioritize the smallest comparisons that address its ambiguities. The complete planned culture set is:

| Study mask | Question |
| --- | --- |
| Growth + karyotype + flow | What does the full score constrain? |
| Karyotype + flow; no growth score | What additional constraint comes from scored proliferation? |
| Growth + flow; no karyotype score | What do direct chromosome counts add beyond growth and flow? |
| Growth + karyotype; no flow score | What do flow-distribution observations add? |
| Growth only | Can distributions separate variation generation and survival beyond growth alone? |

Each mask needs a new unrestricted refit and new nuisance refits at each target. Keep retained modality weights and within-stream denominators unchanged; do not renormalize the remaining weights to sum to one. Share a core natural-scale target grid for comparison and add mask-specific optima/refinements as needed. Audit active nuisances again when a stream is removed.

Dropping a stream changes the objective. Compare where fit quality deteriorates, how compensation changes, and which predictions become unresolved; do not read cross-study score heights as a common significance scale. A removed stream can move or expose a minimum rather than simply broadening a curve.

These are initially **score-deletion information audits**, not independent predictive validation. In particular, deleting the growth score may leave its source counts involved in passage-end selection and downstream propagation. Export a table of every removed observation still used in conditioning. A genuinely held-out prediction requires a separately approved leakage-free observation/propagation protocol. Save omitted-stream predictions as diagnostics without claiming they are independent holdouts.

For tumor studies, compare the full score with removal of burden, terminal viable chromosome distributions, or harvest necrosis, one stream at a time. Retain measured times and scenario identities; dropping a score must not accidentally drop the entire tumor simulation. Clearance may remain active through total burden even when necrosis is omitted.

For each context, first compare data-only against the exact production penalty policy when different. Only then consider a small, prespecified penalty-strength sensitivity. In joint fits decompose context priors, soft coupling, and phenotype constraints separately. If a phenotype constraint enforces the biological direction being tested, a fit retaining it cannot independently establish that direction. Any release of a hard constraint is an explicitly labeled new study.

Do not launch the full target-by-mask-by-penalty Cartesian product without a second compute review. At the pilot settings, each additional two-target configuration adds roughly 100 constrained runs before baseline and validation work.

## 8. WP4: Tumor, Joint, and Function-Level Extensions

Proceed only after culture diagnostics are reliable:

1. **Tumor death and clearance:** Profile `mu_hp` and `k_clear` first, with terminal necrosis present and omitted, while retaining other active latent-resource and initial-condition parameters. If clearance was fixed in production, releasing it is a separately declared sensitivity with its own unrestricted baseline, not a profile of the original fixed-clearance fit. Ask whether the observations separate dead-material production from persistence.
2. **Joint comparisons:** Use the full admissible parameter space across starting families for the primary lower envelope. Compare independently refitted data-only and production-penalized results. Report family-conditional profiles only as additional diagnostics. A narrow curve within one family is not global practical identifiability.
3. **Derived targets:** Profile effective per-chromosome missegregation at declared chromosome/oxygen inputs, the survival increase between 44 and 88 chromosomes, and integrated viable chromosome-loss production over a specified culture interval. Use the configured chromosome unit if different from 22. Report function values and units explicitly.

For a derived target \(c(\theta)\), the required calculation is

\[
Q^{\mathrm{prof}}(z)=\min_{\theta:c(\theta)=z}Q(\theta).
\]

Implement this through a validated reparameterization or an equality-constrained optimizer with reported constraint residuals. An arbitrary finite penalty or a scatterplot of saved endpoints is not a derived profile. Do not profile a function contrast merely by freezing its contributing parameters.

For viable-loss production, define the event ledger before fitting: maternal division rate times the ordinary-division branch probability times the expected number of surviving daughters with fewer chromosomes than the mother conditional on ordinary division, weighted by live maternal counts and integrated over the specified interval. State whether the target is a total expected count or normalized per cumulative division; exclude passage dilution, the one-cell WGD branch, and numerical out-of-grid loss from biological chromosome-loss production. Distinguish this flux from subsequent expansion or endpoint abundance of those descendants.

A well-constrained derived quantity can coexist with weak individual parameters. For example, where \(\mu\ll k\), inducible missegregation is approximately \((a/k)\mu\), so amplitude and half-saturation scale may compensate. Likewise, daughter generation and survival can compensate while their viable output changes little. Test these possibilities rather than assuming them.

Constraining one anchor contrast or one flux is not necessarily equivalent to removing an entire mechanism. Preserve exact mechanism tests from the companion plan and establish any claimed equivalence mathematically. At common model oxygen coordinates, context contrasts remain standardized model evaluations, not matched measurements of intracellular resource stress.

## 9. Reporting and Interpretation

Produce three coordinated diagnostic views for each target:

1. **Profile:** Best verified \(\Delta Q\) versus the target, with exact zero clearly shown, evaluated points visible, and full-data/selected stream-deletion curves identified. Use a zero-capable display or a clearly separated zero point; do not hide zero on a logarithmic axis. No confidence band or probability label in diagnostic mode.
2. **Compensation:** Active nuisance estimates along the lower envelope, highlighting boundary hits, alternative minima, changing initial/error parameters, and selected-passage-day changes. Separate any family-conditional curves.
3. **Biological consequences:** Growth and full distribution predictions, effective missegregation, survival functions, and viable-daughter flux at selected profile points. These are predictions conditional on constrained refits, not posterior predictive intervals.

Show objective components and optimizer dispersion in accompanying QC tables/plots, explicitly as component traces and search diagnostics. All profile minima must link to the underlying parameter vector and predictions.

Classify results descriptively as `two_sided_rise`, `flat_over_tested_range`, `one_sided_unresolved`, `bound_limited`, `penalty_dependent`, or `optimizer_unresolved`, allowing more than one flag. Declare any descriptive score-rise threshold in advance and state that it is not a confidence threshold. A flat numerical profile does not prove structural nonidentifiability; a rising profile does not establish global structural identifiability or mechanism uniqueness.

At exact nulls, report the best refitted score difference, changes in each stream, and whether the observed remodeling remains reproducible. Until statistical calibration, use language such as "the best verified no-induction fit fits these observations less well under the declared score," not "inducible CIN is statistically required."

## 10. WP5: Later Statistical Calibration

This phase requires an explicit statistical design and compute approval. It is not a plotting option on the pilot curves.

### Choose the Inferential Target

Either specify a coherent generative likelihood, or justify component likelihoods and their dependence-aware calibration. Arbitrarily averaged mixed losses do not become a composite likelihood by renaming them. Even a valid composite likelihood generally needs adjusted inference rather than ordinary likelihood-ratio cutoffs; sensitivity and score-variability information need not coincide. [Varin, Reid, and Firth, 2011](https://www3.stat.sinica.edu.tw/sstest/oldpdf/A21n11.pdf).

For a generative or simulation-calibrated approach, specify observation errors and any process variation for the actual design: repeated passage growth within related lineages; physical chromosome samples; flow acquisition/density estimation rather than independent density-grid points; repeated tumor burdens; and harvest necrosis with its mouse/section sampling structure. Include only supported variance components and state which cannot be learned from the available replication. Do not invent missing assay counts.

Resampling units must follow the biological design and shared ancestry. Cells within one sample, technical acquisitions, or optimizer seeds must not be substituted for independent lineages or tumors. Small numbers of independent units may prevent reliable empirical composite-likelihood adjustments; assess that limitation before selecting the method.

### Calibration and Coverage Protocol

1. Define the statistic and its scale explicitly, for example \(T(z)=Q^{\mathrm{prof}}(z)-\widehat Q\), or twice that difference for a declared negative log likelihood. Do not mix conventions across datasets or penalty policies.
2. At a proposed target value, fit the restricted nuisance parameters under the selected inferential model. Generate replicate datasets at the real schedules and sample sizes, including declared dependence and the complete observation-processing pipeline. If counts influence passage selection, reproduce that dependence in the generator and fitting procedure.
3. For every synthetic dataset, redo unrestricted and target-constrained optimization with the same bounds, parameter activity rules, and validated search budgets. Calibrate the statistic's cutoff from these simulations, not from optimization-endpoint dispersion.
4. Treat exact-zero nulls separately. Parameters such as the inducible half-saturation scale become inactive there, so do not assume an ordinary chi-square cutoff or a universal boundary-mixture correction.
5. Check calibration across nuisance choices and competing solution regions, not just one plug-in optimum. Include exact-null, near-null, interior, weak-information, and boundary-limited scenarios. State clearly where any calibration is conditional on a particular generator or nuisance choice.
6. Validate coverage on an independent synthetic batch. Report coverage, false rejection at zero, interval/set width, disconnected or unbounded sets, bound hits, and optimizer failures, with Monte Carlo uncertainty. Count and investigate failures rather than dropping them from the denominator.
7. Begin with a small debugging batch without publishing confidence claims. Set final replicate counts from the required precision of quantiles and coverage, with a new cost estimate. For scale, 500 independent coverage trials give roughly two percentage points of 95% Monte Carlo half-width near 95% coverage; reliable tail-cutoff estimation may need more.
8. Invert the validated tests to form intervals or sets without forcing a single connected interval. If retaining the original score, call the output a simulation-calibrated compatibility set under the declared generator, not automatically a profile-likelihood confidence interval. Record the calibration version alongside every reported limit.

A change in observation model or score requires new fitted baselines and profiles. Existing objective curves cannot simply be relabeled as likelihood intervals. It is acceptable to finish with useful diagnostic profiles and a documented reason why calibrated intervals remain unsupported.

## 11. Scoped Implementation and Tests

### Proposed File Changes

These are future edits, not files created by this plan:

| Location | Responsibility |
| --- | --- |
| Existing `optimizer/profile_likelihood_O2_supply_demand_MAP.R` | Add explicit `objective_diagnostic` mode and context/config dispatch, with no inherited CI behavior in that mode. Preserve explicit legacy behavior and archived outputs. |
| New `util/o2_supply_demand_map_profile.R` under `oxygen/code/O2_supply_demand_MAP/` | Shared study validation, grid/start generation, profile orchestration, candidate acceptance, minimum reduction, checkpointing, and diagnostic outputs. Extract reusable helpers without unrelated optimizer refactoring. |
| Culture `objective.R` and relevant fit backend adapters | Shared active/fixed/inactive parameter mapping, exact zeros, actual warm-start injection, and frozen observation masks. Reuse the companion plan's mapping implementation if already available. |
| In vivo/joint backend adapters, in the later stage | Immutable prior/model settings, complete objective-component export, context-specific nuisance mapping, and exact target constraints. |
| New `oxygen/tests/testthat/test-profile-objective.R` | Fast analytic/mock tests; extend existing kernel tests for biological nulls. |
| Study configuration under the analysis output root | Parameter grids, source/input references, score policies, tolerances, search budgets, and execution manifest. No edits to production parameter CSVs. |

Keep the entry point thin and use existing fitting backends rather than implementing a second optimization stack. Defer a dedicated likelihood-calibration module until WP5's observation model and method are approved.

### Required Fast Test Oracles

| Test | Expected result |
| --- | --- |
| Nuisance refitting: \(Q(t,u)=(t+u-3)^2+(u-1)^2\) | \(u^*=2-t/2\); profile \((t-2)^2/2\). Holding \(u\) fixed must not pass. Use bounds containing the analytic optima. |
| Generation-survival ridge: \(Q(q,s)=(qs-c)^2\), \(q,s\in[0,1]\), \(0<c<1\) | Profile \([\max(c-q,0)]^2\); flat for \(q\ge c\), with a bound-induced rise below. |
| Penalty removal: add \(\lambda(s-s_0)^2\) to the ridge | Where interior, \(s^*=(qc+\lambda s_0)/(q^2+\lambda)\). Data-only refitting recovers the ridge; subtracting the penalty from penalized solutions generally does not. |
| Stream deletion: \(Q_1=(t-1)^2\), \(Q_2=(t+1)^2\) | Full-data baseline at zero; deletion baselines at +1 or -1. Every study's delta has its own zero minimum. |
| Competing mocked basins | Results \((t-1)^2\) and \((t+1)^2\) produce their pointwise minimum, invariant to seed order; not their mean. Failed/infeasible candidates cannot win. |
| Derived constraint: \(Q(a,b)=a^2+b^2\), \(a+b=z\) | Profile \(z^2/2\), attained at \(a=b=z/2\), with the constraint satisfied numerically. |
| Objective rescaling | Multiplying a score by a positive constant preserves optima but scales profile heights. Diagnostic mode must not silently reuse a confidence cutoff. |

Additional integration tests must verify exact-zero kernel behavior; full parameter reconstruction; data-inactive versus penalty-active handling; unchanged priors/model caps when only the target is fixed; real start-vector injection; fixed sample masks under numerical failures; wrong-target/stale/interrupted-result rejection; baseline-improvement recovery; fingerprint-checked resume; and absence of CI fields/lines/stopping rules in diagnostic outputs. Production replay is required again after integration.

## 12. Configuration, Outputs, and Completion Gates

Use configurable absolute input paths and an isolated output root, for example:

```text
agent-dev/major_analyses/<date>_profile_identifiability/<study_id>/
  config/       study.json, profile_grid.tsv, parameter_status.tsv
  manifests/    production_contract.json, inputs.tsv, observations.tsv
  baselines/    unrestricted fits and replay checks
  fits/         <target>/<point_id>/<start_id>/<attempt_id>/
  tables/       profile_points.tsv, profile_starts.tsv, nuisance_paths.tsv
                objective_components.tsv, predictions.tsv, derived_outputs.tsv
  figures/      profile, compensation, prediction, and optimizer-QC views
  qc/           validation_report.md, conditioning_audit.tsv, failures.tsv
  reports/      findings.md, runtime.md, rebuild_commands.txt
```

Required configuration fields: `study_id`, context, model variant, source/container fingerprints, input paths, observation mask, immutable model settings, parameter domain/transforms, exact targets, penalty definitions, modality weights/normalizations, start-bank paths, master seed, grid specification, target/score/search tolerances, optimizer budgets, restart/refinement limits, cores/memory/walltime, output root, and resume policy. Diagnostic mode must explicitly set inference calibration to absent. Reject incomplete scientific settings rather than guessing them at launch.

Minimum table contracts:

- `profile_starts.tsv`: study/target/point/start/attempt IDs, requested and reported target, actual initial-vector path/hash, parameter-result path, seed, fingerprint, all objective terms, reevaluation difference, scored-sample counts, target residual, optimizer/launch status, evaluations, runtime, scheduler ID, and acceptance reason.
- `profile_points.tsv`: natural target value, best verified score, matched unrestricted baseline ID/value, delta, winning run, attempted/accepted start counts, restart-agreement diagnostics, boundary/search flags, and component-trace links. Preserve unresolved points as unresolved.
- `nuisance_paths.tsv`: winning full parameter vectors with units, transforms, activity status, and bound indicators.
- `predictions.tsv`: physical sample IDs, modality, observed/predicted values or distributions, scored/omitted/conditioning status, and selected simulated passage day where applicable.
- `derived_outputs.tsv`: precisely defined function values, survival contrasts, and fluxes with evaluation inputs, units, and interval/normalization definitions.

Store hashes and rebuild commands for all outputs. Keep bulky fit directories outside commits unless separately approved; commit compact plans/configuration/provenance through scoped changes. Do not overwrite historical reports, replace manuscript figures, or edit Results as part of this implementation.

### Ordered Approval Gates

1. **Contract:** Approved production snapshot, objective ledger, exact-null domain, and successful replay.
2. **Engineering:** Fast tests and one constrained-fit smoke test pass; measured resource estimate approved.
3. **Culture pilot:** Two full-data targets, about ten values and five starts per value, plus baselines and targeted validation. Review curves and optimization adequacy before expansion.
4. **Information audit:** Refit selected stream-deletion and penalty studies, then complete the planned masks as warranted. Report score deletion separately from predictive validation.
5. **Context/function extension:** Approve tumor, joint, and derived targets based on the ambiguities revealed, with a fresh compute estimate.
6. **Statistical inference:** Approve a defensible observation model/calibration design and validate coverage before reporting confidence limits.

The first useful endpoint is a short scientific report stating which mechanisms or combinations are constrained under the declared score, what compensation remains possible, which measurements supply the constraint, and what new measurement would resolve the ambiguity. It need not wait for a full-model structural proof or confidence intervals for every parameter.
