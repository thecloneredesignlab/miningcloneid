# Meeting-aligned identifiability master plan: structural analysis before synthetic recovery

Date: 2026-11-12, America/New_York.

Scope reviewed: 2026-11-13, America/New_York. This plan covers the culture and tumor ploidy-evolution model and its identifiability analysis only.

Sequencing revised: 2026-11-13, America/New_York, following Noemi's clarification that structural-identifiability results should inform the design of synthetic practical-identifiability experiments.

Status: **Revised documentation for structural-first planning.** The earlier synthetic-first sequencing is superseded by this clarification. This document update authorizes no implementation, numerical experiments, HPC submissions, controller changes, or new compute budget.

## Purpose and governing priorities

The first scientific task is structural identifiability: identify theoretical parameter ambiguities and identifiable combinations under explicit model, input, observation, and initial-condition assumptions. These results must inform which quantities the subsequent synthetic practical-identifiability study attempts to recover.

Noemi's later clarification, supplied by the user on 2026-11-13, governs the sequence wherever it differs from the earlier AI meeting summary or the preceding message. The two original plans supply the mechanism questions, observation contracts, and reduced-model starting point. The revised order is:

1. Audit culture sources, all fourteen mechanistic parameters, nuisance roles, and the actual experiment.
2. Formulate and validate the structural-analysis model and observation/input cases; establish a working symbolic or analytical route.
3. Identify structural ambiguities, identifiable parameters/combinations, and unresolved cases, with their exact assumptions.
4. Map those findings back to the full simulator and review the resulting synthetic estimands, parameterization, observation sets, and initial-state assumptions. This is a release gate for designing and running the practical study.
5. Only then register generating truths, implement the reviewed finite-sampling design, and test deterministic and realistic recovery.
6. Expand to repeated datasets, further regimes, paired stream ablation, and experimental-design comparisons under finite budgets.
7. Use focused profiles and function/flux recovery to explain remaining practical ambiguity; undertake formal interval calibration later.
8. Extend to tumors with their own structural/input-output review before designing tumor synthetic recovery.

Source and assay audits may proceed together, and basic artifact inspection can continue during tool discovery. These supporting tasks do not bypass the structural review or finalize the synthetic recovery design. If the structural route is unresolved, stop that decision at review rather than automatically proceeding under the former synthetic-first policy.

The existing real-data profiles are supporting evidence. They do not substitute for synthetic recovery, and the current study must not become another open-ended search for a global optimum on the real data. A reproducible feasible fit can supply a generating parameter set even when independent recovery of the real-data minimum remains unresolved.

Successful simulation and refitting under the same model establishes conditional self-consistency, not the biological correctness or necessity of that model. The original oxygen-linked-missegregation, pre-existing-selection, and shared-versus-released-response alternatives remain later mechanism comparisons guided by the recovery findings. They are not silently discarded, and they are not added to the first synthetic batch. Any generator/estimator misspecification study is an explicitly registered later extension rather than an unreported change to the primary generator.

## Evidence, precedence, and current limitations

### Primary requirements

- Noemi's later clarification, supplied by the user on 2026-11-13: start with structural identifiability so that theoretical non-identifiabilities and parameter combinations inform the design of the synthetic practical-identifiability study. This clarification supersedes contradictory sequencing in the earlier summary and this document's original version.
- Meeting summary supplied by the user: /Users/4482173/.codex/attachments/74d1a2bd-650d-41ff-aa6b-ce908738dbfd/Pasted text.txt. Retain its finite-sampling, recovery, measurement-contribution, and design requirements, interpreted with the structural-first clarification.
- Earlier Noemi message: begin with a small culture profile pilot for inducible missegregation and ploidy-dependent survival; structural work can follow or proceed alongside it; defer formal confidence intervals until statistical calibration. Preserve this as historical context, not the current sequencing instruction.
- [Original profile plan](profile_likelihood_implementation_plan.md).
- [Original mechanism and structural-identifiability plan](invitro_mechanism_identifiability_implementation_plan.md).

All fourteen mechanistic parameters remain in the inventory. Structural analysis determines which should be targeted individually, through combinations, or only as explicitly unresolved diagnostics. The later practical study evaluates finite sampling, synthetic recovery, and measurement contribution for those reviewed targets; it must not require unique recovery of parameters already shown to be structurally confounded.

### Reviewed recent documentation

The review window is 2026-11-11 through 2026-11-12 in America/New_York, using creation or modification time rather than filenames alone. The twelve existing documents below qualified. Times describe the local files inspected for this plan, not live HPC state.

| Document in docs/ | Last local modification, EST | Use in the new plan |
| --- | --- | --- |
| culture_profile_optimization_review_plan.md | Nov 11, 01:33 | Search limitations, compensation questions, and bounded review history |
| culture_baseline_recovery_nojson_plan.md | Nov 11, 10:57 | Verified culture mapping, search orchestration, and provenance practices |
| culture_baseline_recovery_execution_log.md | Nov 11, 23:20 | Completed recovery evidence and boundary limitations |
| culture_seven_point_profile_nojson_plan.md | Nov 12, 00:02 | Exact-null semantics and constraints of the seven-point review |
| culture_seven_point_profile_execution_log.md | Nov 12, 18:36 | Final scientific outcome, diagnostic checks, accounting, and collection |
| culture_resolution_aware_stopping_plan_and_handoff_20261112.md | Nov 12, 12:39 | Separation of numerical integrity and scheduling resolution |
| culture_resolution_aware_revision_review_and_execution_plan_20261112.md | Nov 12, 13:00 | Vector-specific evidence, actual-origin recovery, and migration safeguards |
| invivo_identifiability_parallel_handoff.md | Nov 11, 15:12 | Tumor observations, scale assumptions, and ownership isolation |
| invivo_baseline_repair_plan_20261111.md | Nov 12, 03:33 | Known feasible points, optimizer failures, and scale compensation |
| invivo_improvement_driven_continuation_plan_20261112.md | Nov 12, 09:31 | Checkpoint and scoring contracts; its unbounded dispatch policy was subsequently revised |
| invivo_resolution_aware_stopping_plan_and_handoff_20261112.md | Nov 12, 14:22 | Finite 220/360-core-hour policy for the existing tumor study |
| invivo_resolution_aware_review_and_next_steps_20261112.md | Nov 12, 15:20 | Consolidated tumor policy implementation and accounting details |

The two original plans remain explicit reference documents even though they fall outside this filesystem-time inventory. Preserve all original documents and approved execution snapshots.

### Dated status relevant to planning

- The earlier culture baseline recovery obtained cross-family agreement under its own final reference. The later seven-point review found a lower reference and ended with independent reference recovery unresolved. The final B gap was approximately 0.00441833 against a 0.001 requirement. No neighbor continuation or final point repeat was released; zero of seven full point gates passed, although four passed their independent-point component. The journal records final verified delivery at 18:35:53 EST on Nov 12 and 1382.551944444 actual core-hours for that stage. These are distinct stages, not contradictory outcomes.
- Eight selected culture candidates passed diagnostic kernel, trajectory, and event-accounting checks. Their functions and predictions can guide the new study, with the original optimization flags retained.
- A local in vivo report updated at 23:27:57 UTC on Nov 12 still reported both baselines not accepted and profiles pending. Printed best scores were approximately 14.019672899630 and 13.409650415535, with independent-uniform gaps approximately 0.178941391630 and 0.124792185787. A report's closed-parent snapshot does not establish current task completion; later local state may contain queued work. Only the execution owner should reconcile live HPC state.
- This document neither restarts the closed culture study nor pauses, cancels, or expands the existing tumor study. Its prior 220/360 allowances and culture's prior allowances are not a budget for new synthetic experiments.

## Study-wide contracts

### Parameter roles

Verify this table against the frozen production-v2 source and actual runtime configuration before implementation. The inspected culture optimizer has twenty coordinates, not fourteen total unknowns.

| Role | Natural-scale parameters | Count |
| --- | --- | ---: |
| Division | lam_max | 1 |
| Chromosome variation | p_mis_base, p_misseg, k_o_mis, p_wgd | 4 |
| Daughter survival | buffer_smax, buffer_beta, buffer_n_exp | 3 |
| Oxygen-linked growth | alpha_o2, gamma_growth | 2 |
| Stress/death response | mu_hp, gamma_mu, O2_crit, n_O | 4 |
| Initial composition nuisance | init_mean_2N, init_sd_2N, init_mean_4N, init_sd_4N | 4 |
| Fitted observation nuisance | sigma_growth, sigma_kary | 2 |

The first five rows inventory the fourteen mechanistic parameters. Step 05 releases the actual individual or combination estimands after structural review. The primary realistic recovery configuration must retain nuisance estimation appropriate to the actual study. Known-nuisance configurations are diagnostic controls, not evidence that the parameters are recoverable when those nuisances are unknown. The flow kernel has separate configured/derived behavior; do not invent an additional fitted flow-error coordinate. Parameters that govern observation variance require an appropriate probabilistic observation analysis; a noiseless mean-output calculation alone cannot establish their identifiability.

At exact p_misseg=0, k_o_mis is inactive for the culture data-only model; at buffer_beta=0, buffer_n_exp is inactive. Do not demand recovery of an inactive parameter. Additional assay-error or dead-material parameters require explicit roles and assumptions rather than silently extending the fourteen-parameter claim.

### Experimental units and observation generation

Keep latent population dynamics, physical samples, technical observations, and optimizer starts separate. A large latent population does not imply a large observed sample. Density-grid points are not sampled cells. Shared lineage ancestors and reused cell counts can induce dependence.

Before synthesis, resolve whether passage timing is fixed, count-triggered, or a mixture. The inspected fitter selects a simulated day near an observed target live count and propagates that composition. An observed terminal count cannot be held fixed as a noiseless conditioning input while an independent noisy version is generated as though unrelated. Reproduce the actual joint stopping/count/time process or explicitly label a conditional approximation. Any discrepancy between the physical protocol and the production fitter must be exposed, not hidden by calling both schedules identical.

The primary generator follows the verified experiment. Preserve a production-compatible estimator for an explicit estimator-recovery study; if a faithful-time or other protocol alternative is needed, register it as a separate configuration with its own validation and interpretation.

### Scoring and inference

The inspected culture objective combines modality-weighted averages of growth, sample-level karyotype, and flow-density scores. It is not automatically a normalized joint likelihood. The first synthetic study may evaluate recovery of this existing estimator; it need not wait for a new inferential likelihood. It does require a specified generative observation model.

Distinguish four quantities throughout:

1. Numerical artifact validity and repeatability.
2. Adequacy of finite optimizer search.
3. Parameter and prediction recovery across independently generated datasets.
4. Formal confidence-set coverage, which belongs to later calibrated inference.

Retain all failures and native limits. Never select a fit because it is closest to the known truth. Do not interpret optimizer-start dispersion as a sampling distribution or confidence interval.

## Step 01 - Reconcile sources, ownership, and reusable evidence

### Inputs

Noemi's structural-first clarification, the relevant ploidy-model sections of the meeting summary and earlier message, the fourteen source documents above, completed culture packages, current tumor metadata, and their manifests.

### Execution steps

1. Build a requirement-to-analysis table showing how each applicable ploidy-model identifiability request is implemented or deferred, with structural-first precedence explicit.
2. Inventory the actual source, image, data, observation manifests, fits, numerical checks, failed gates, and collection evidence available for reuse.
3. Record ownership for the initial culture structural analysis, its full-model bridge, later synthetic recovery, and existing tumor execution. Structural analysis is part of the main scientific sequence. Do not launch another dispatcher against an owned run.
4. Preserve the closed culture outcome. Ask the existing tumor owner to reconcile its current work when an implementation handoff is separately authorized; this planning step itself changes no jobs.
5. Define a new synthetic namespace and ledger rather than modifying old scientific payloads.

### Outputs

source_inventory.tsv, requirement_traceability.tsv, ownership_and_scope.md, reusable_artifacts.tsv, and a dated status snapshot.

### Acceptance gate

Every primary requirement has an assigned step, every reused numerical object has provenance, and completed, unresolved, and live work are distinguished. Existing execution allowances are not represented as unused synthetic-study funds.

### Stop conditions

Hold dependent implementation for missing source identity or conflicting ownership. Do not repair these problems by overwriting old runs or assuming a local monitor reflects current HPC state. Independent documentation and data discovery may continue.

## Step 02 - Audit parameters, nuisance roles, and candidate estimands

### Inputs

The frozen optimizer specification, natural-scale parameter table, model switches, initial-state mapping, observation operators, and Step 01 inventory.

### Execution steps

1. Verify the fourteen mechanistic names, units, transformations, bounds, and actual runtime use.
2. Enumerate the four initial-distribution and two observation-error coordinates, plus every fixed setting.
3. Inventory known-versus-unknown initial-state and nuisance cases for structural analysis. Draft, but do not yet release, known-nuisance and realistically free-nuisance recovery configurations.
4. Record candidate parameter-specific error scales and summaries. Finalize individual versus combination estimands in Step 05 using structural findings, before generating recovery datasets. Relative error at zero is undefined.
5. Inventory important function/prediction targets, including induced error response, survival contrasts, and viable-daughter output. Audit active/inactive roles in every later mask or null.

### Outputs

parameter_roles.tsv, fixed_settings.tsv, candidate_estimands.tsv, structural_nuisance_cases.tsv, and parameter_map_validation.md.

### Acceptance gate

Counts reconcile to the actual executable mapping. Known and unknown quantities are explicit for structural cases, and candidate recovery targets remain provisional until Step 05. No fixed convention or inactive coordinate is presented as an identified parameter.

### Stop conditions

Hold dependent analysis if parameter semantics, units, active roles, or source versions disagree. Do not reduce the parameter inventory merely because a tool is difficult to run or recovery is difficult. A necessary new nuisance parameter must be declared before interpreting structural or practical results.

## Step 03 - Audit the actual experiment and observation contract

### Inputs

Culture fit objects, lineage metadata, passage records, karyotype cell lists, growth/count definitions, flow inputs and acquisition metadata, and any live/dead records.

### Execution steps

1. Enumerate physical samples and their cohort, ancestor, oxygen, passage, time, modality, and measured sample size. Record missing values explicitly.
2. Separate exogenous settings from stochastic observations and fitted quantities. Trace how observed counts affect selected passage days and downstream states.
3. Verify whether growth errors share underlying count measurements across passages; identify biological and technical replicate units.
4. Locate raw flow event counts, gating, and density processing where available. Distinguish acquisition count from retained G0/G1 events and from grid length.
5. Determine whether live/dead counts exist. Define assay timing, denominator, detached-cell handling, and any classification error; verify correspondence with model states.
6. Resolve fixed-time versus count-triggered passage design. Specify any conditional approximation or separate protocol sensitivity needed.

### Outputs

observation_manifest.tsv, experimental_schedule.tsv, conditioning_graph.md, modality_availability.tsv, dead_cell_observation_contract.md, and unresolved_inputs.tsv.

### Acceptance gate

Each available modality has a reproducible sampling/processing definition. The generator and fitter's use of counts, timing, and shared ancestry is explicit. A typical karyotype count of twenty is not substituted for sample-specific metadata without a labeled assumption.

### Stop conditions

Block claims of faithful real-design recovery for a modality whose sampling or protocol mapping is unresolved. Missing live/dead records do not block source/assay audits or structural cases using verified existing outputs; an additional dead-cell output remains explicitly hypothetical. They do not waive the structural-first release gate. Do not fabricate assay counts or treat tumor necrosis as culture live/dead data.

## Step 04 - Perform culture structural-identifiability analysis first

### Inputs

Steps 01-03 parameter/input/output inventories, the original five-live-class reduction with retained dead material, exact full-model equations and switches, and the actual citation for the proposed Mathematica package.

### Execution steps

1. Obtain the exact package/paper, establish access and license requirements, and verify supported model classes using authoritative documentation. Do not assume that the Mathematica package is COMBOS or that either accepts the complete simulator unchanged. Validate a known identifiable and a known nonidentifiable example before project-specific claims.
2. Formulate a tractable culture model, using the original five-live-class recurrent model as the starting point. Preserve the relevant error-generation, daughter-survival, WGD, growth/death, and dead-material accounting. Record rational response substitutions, fixed exponents, coarse states, and omitted mechanisms in a reduction map; the one-generation ledger is only an algebraic unit check.
3. Specify ideal outputs, inputs, and initial conditions for each case. Compare controlled constant-oxygen experiments sharing parameters, justified switching inputs, known versus unknown initial mixtures/abundance, and applicable passage/reset assumptions. A proof for an arbitrary varying input or continuously observed output must not be labeled a result for the actual sparse schedule.
4. Analyze an output ladder: growth/live abundance alone; abundance plus viable composition; and additional dead-cell or event/fate outputs only as declared observed or hypothetical cases. Remove redundant normalized fractions and distinguish live abundance from total retained material.
5. Classify globally identifiable, locally identifiable, nonidentifiable, and unresolved/tool-limited quantities under each case. Export explicit identifiable combinations or input-output-preserving transformations where obtained. Treat exact-zero boundaries and inactive parameters separately from generic interior conclusions.
6. Verify claimed equivalences by analytical substitution or other appropriate structural evidence. Numerical trajectory agreement and rank checks can corroborate or nominate candidates, but are not global structural proofs.
7. Use a finite, pre-registered case/time/memory allowance. If one package is unavailable, consider a validated alternative symbolic or analytical route within the same structural scope. Preserve raw inputs, versions, outputs, failures, and assumptions.

### Outputs

tool_access_and_compatibility.md, symbolic_models/, symbolic_case_manifest.tsv, reduction_validation.tsv, structural_identifiability_results.tsv, identifiable_combinations.md, structural_transformations.tsv, and structural_limitations.md.

### Acceptance gate

A structural results package has been reviewed before practical-study design. Each conclusion identifies its model, outputs, inputs, initial conditions, parameter domain, and proof/status. Theoretically identifiable parameters/combinations and unresolved quantities are explicit. A tool smoke test or a result for an unrelated toy model does not satisfy this project-specific gate.

### Stop conditions

A package timeout, inaccessible license, or intractable reduction is an unresolved structural case, not proof of nonidentifiability and not automatic permission to proceed with recovery. At the finite allowance, stop and review an alternative route, a justified model reduction, or a separately approved restricted diagnostic scope. Source and assay audits can continue; the former rule that practical recovery cannot be delayed is withdrawn. Do not silently drop mechanisms, fix unknowns, or change the pinned production model to obtain a convenient result.

## Step 05 - Translate structural results and release the practical-study design

### Inputs

The reviewed Step 04 results and assumptions, the full simulator and observation contract, the parameter/nuisance inventory, and the actual finite sampling schedule.

### Execution steps

1. Map each reduced-model parameter and identifiable combination to the original fourteen mechanisms and relevant nuisances. Label exact correspondences, approximations, fixed shape assumptions, and unmapped quantities. Neither identifiability nor nonidentifiability of a reduction transfers automatically to the full simulator.
2. Where structural equivalence is demonstrated for the relevant full input-output map, choose identifiable coordinates/combinations for primary recovery. If a representative parameter is fixed to break an equivalence, label that as a convention or external-information assumption, not identification by the data.
3. For findings supported only by the reduction, define targeted full-model checks and retain their limited status. Numerical checks at actual observation events address local distinguishability, not a global structural proof. Do not remove original parameters solely on the strength of an unvalidated reduced correspondence.
4. Build the final estimand table: individually estimable targets, identifiable combinations, inactive quantities, parameters retained only for diagnostic comparison, and unresolved cases. Keep all fourteen original mechanisms accounted for; do not demand unique recovery from a structurally confounded coordinate.
5. Use the structural output ladder to choose baseline modalities, informative additions, known/unknown nuisance cases, and initial generating regimes. Preserve realistic finite sampling and noise; theoretical identifiability under ideal outputs does not establish practical recovery.
6. Register the parameterization, observation/conditioning design, truth-selection rule, error metrics, finite search strategy, and staged compute proposal before implementing the synthetic generator or launching new recovery fits.
7. Record an explicit scientific review decision and the unresolved limitations. A narrowed practical diagnostic after incomplete structural work requires an explicit scope decision; package inconvenience is not a standing bypass.

### Outputs

structural_to_full_model_map.tsv, recovery_estimands.tsv, nuisance_configurations.tsv, estimand_decision_log.tsv, structural_informed_synthetic_design.md, full_model_bridge_checks.tsv, and structural_review_gate.tsv.

### Acceptance gate

The structural-to-practical review has accepted a coherent design for the specified targets and model. Its coordinate changes, external assumptions, applicable findings, and remaining limitations are explicit. Only after this gate may Steps 06-14 implement and run that practical design, subject also to the separate implementation and compute authorization.

### Stop conditions

Hold practical-design release for an unexplained parameter mapping, a reduction that changes the mechanisms of interest, or an unsupported claim of structural identification. If the relevant case remains unresolved, stop for review rather than returning automatically to an all-parameter recovery study. Do not select favorable truths or noise levels to conceal a structural ambiguity.

## Step 06 - Register the first generating truth and simulate latent dynamics

### Inputs

Verified existing culture fits, Steps 02-05 contracts and accepted structural review, the pinned simulator, and full initial/error settings.

### Execution steps

1. After Step 05 release, choose one reproducible existing full parameter vector consistent with the reviewed estimands and representative conventions. Use provenance and scientific representativeness, not favorable synthetic recovery. Disclose unresolved global-search status, structural-equivalence assumptions, and boundaries.
2. Freeze all generating coordinates, fixed settings, initial states, source hashes, and any conditioning assumptions. Do not choose only fourteen values while omitting the rest of the generator.
3. Simulate the complete lineage and passage history under assigned oxygen and the audited protocol. Preserve latent live states and required dead/event outputs.
4. Verify state propagation, reseeding, ordinary division, one-lineage WGD, and event accounting. Separate numerical out-of-grid loss from biological nonviability.
5. Define whether this first experiment uses deterministic expected dynamics with observation noise only. Any later biological/process variability is a separate declared model, not hidden in independent measurement draws.

### Outputs

truth_registry.tsv, truth_parameters.rds, truth_configuration.tsv, latent_trajectories.rds, latent_event_audit.tsv, and truth_selection.md.

### Acceptance gate

The complete generating object replays with the frozen model and yields feasible trajectories under the declared protocol. A globally optimal real-data fit is not required. No synthetic outcome has influenced truth selection.

### Stop conditions

Hold a truth with irreproducible provenance, invalid dynamics, or unresolved observation-state correspondence. Do not silently replace it with a favorable case. Valid boundary truths remain usable with explicit boundary interpretation.

## Step 07 - Implement realistic finite-sample observation generation

### Inputs

The Step 05 accepted synthetic design, registered latent trajectories, physical-sample manifest, nuisance truth values, assay processing, and independently registered data-generation RNG streams.

### Execution steps

1. Generate karyotype cells/counts using each sample's audited size and the declared chromosome-observation kernel. Preserve the distinction between true state and measurement error.
2. Generate passage growth using the actual definition and time/count process. Prefer raw-count generation when required to reproduce shared errors; otherwise label and justify a passage-level error approximation.
3. Generate live/dead counts through the verified assay operator, using finite counts rather than exact latent fractions. A conditional binomial model may be an initial assumption only when appropriate; document dependence or misclassification sensitivity separately.
4. For flow, sample events and reproduce gating/normalization/density estimation where metadata supports it. Distinguish physical measurement spread from smoothing bandwidth to avoid double-counting noise. If only a processed density is available, label an approximation and vary unsupported assumptions later.
5. Keep assigned culture oxygen as an input. Do not create an additional noisy oxygen observation or claim observed intracellular stress.
6. Export artificial data in the fitting pipeline's expected schema, retaining true latent values separately. Verify finite counts, probabilities, IDs, missingness, and empirical generator behavior over cheap observation-only draws.

### Outputs

observation_generator_contract.md, synthetic_dataset.rds, synthetic_observation_manifest.tsv, sampling_audit.tsv, generator_checks.tsv, and RNG_manifest.tsv.

### Acceptance gate

Every generated observation maps to a declared physical sample and sampling process. Finite sample counts are respected; exact latent distributions are not substituted for realistic data. Proposed dead-cell or flow scenarios remain distinguishable from verified existing measurements.

### Stop conditions

Stop a modality generator on inconsistent denominators, impossible counts, duplicate sample IDs, unsupported noise, or hidden reuse of observed outcomes as fixed inputs. Do not run the large refitting batch before the generator passes. Observation-only debugging may continue without allocated model fitting when later implementation is authorized.

## Step 08 - Test deterministic/noise-free recovery

### Inputs

The frozen truth, deterministic expected observations, the fitted observation mapping, declared numerical tolerances, the known-nuisance diagnostic configuration, and the upstream structural expectations for individual parameters or combinations.

### Execution steps

1. Construct expected observations without finite-sample randomness. If the existing interface only accepts integer cells, use a clearly labeled expected-observation engineering adapter rather than calling a large random sample noise-free.
2. Fix or separately control observation-error scales in this diagnostic so shrinking a fitted scale cannot manufacture an apparent recovery result. Do not set a positive-log nuisance to an invalid zero.
3. Replay the generating vector and independently recompute its predictions and score. An expected score need not be numerically zero.
4. Run a finite blind multi-start recovery test. Keep any truth-initialized positive control separate from eligible independent starts.
5. Compare observation recovery, objective attainability, and recovery of the structurally reviewed parameters/combinations. Check alternative equal-quality vectors against predicted equivalences and reduction limitations. Investigate differences without tuning toward the known answer.

### Outputs

noisefree_case.tsv, truth_replay.tsv, noisefree_fit_results.tsv, prediction_identity.tsv, alternative_solutions.tsv, and noisefree_diagnosis.md.

### Acceptance gate

The generator/fitter bridge and attainable prediction/score are reproduced to predeclared numerical tolerances. Failure of unique parameter recovery does not by itself fail this gate when different vectors reproduce the same observations; it becomes a documented ambiguity for subsequent study.

### Stop conditions

Do not proceed to expensive noisy replication if blind fits cannot approach a known feasible generating solution and the discrepancy is unexplained. Distinguish implementation, optimization, protocol, and possible model-equivalence causes. A finite numerical test cannot establish a global structural symmetry.

## Step 09 - Fit one realistic synthetic dataset and benchmark cost

### Inputs

One frozen realistic dataset, Steps 02-08 validation including the structurally informed design, the released nuisance/coordinate configuration, and a separately authorized pilot compute allowance.

### Execution steps

1. Pre-register a finite start bank, seeds, bounds, actual search controls, and one bounded rescue allowance. Four independent complete starts are a planning candidate, subject to measured cost and the approved implementation budget.
2. Estimate the Step 05 individual parameters or combinations and required nuisance coordinates, retaining all mechanisms in the forward model. Report the status of all fourteen original parameters. Any fit in redundant original coordinates is an explicitly labeled diagnostic, not a requirement to recover each coordinate uniquely. Known-nuisance runs remain separately labeled controls.
3. Keep truth-based scoring or initialization as separate diagnostics. On noisy data the fitted optimum can outperform the generating vector and differ from it; that is not automatically failure.
4. Choose the primary estimate by the lowest verified feasible score, never by parameter distance to truth. Preserve competing fits and original termination/evidence.
5. Measure whole-chain allocated cost, evaluations, wall time, memory, feasibility, and checks. Separate optimizer variation within this one dataset from later sampling variation.
6. Freeze dataset-level adequacy rules and finite stopping before confirmation replicates. Do not transplant real-data 0.001/0.002 tolerances without checking the new score and desired recovery resolution.

### Outputs

single_dataset_fit_index.tsv, recovered_parameters.tsv, optimizer_quality.tsv, prediction_recovery.tsv, measured_cost.tsv, and pilot_review.md.

### Acceptance gate

The realistic dataset is correctly consumed, numerical evidence is valid, and remaining optimizer uncertainty is measured and compatible with the intended recovery comparison, or explicitly classified unresolved. Unbiased recovery of every parameter is not required to advance a valid diagnostic study.

### Stop conditions

Stop at the registered pilot/rescue limits rather than repeatedly adding starts. Withhold parameter-identifiability conclusions dominated by unresolved optimization. Stop immediately on source/target/score identity failure. A budget-limited estimate can remain a provisional result with its limitations.

## Step 10 - Repeated independent synthetic datasets at one truth

### Inputs

The validated generator and estimator, frozen truth, dataset-level search protocol, measured pilot cost, and a registered replicate manifest.

### Execution steps

1. Start with a small development/debugging batch, provisionally about twenty independent datasets. Decide a larger batch, for example about one hundred or more, from Monte Carlo precision and measured cost rather than treating those counts as automatic authority.
2. Separate data-generation seeds from optimizer seeds. Give every dataset its own observation draw and complete refitting record.
3. Freeze the estimator after development corrections. New confirmation datasets use untouched seeds; development results are not relabeled as prospective confirmation.
4. Report bias, natural/log/scaled error, RMSE or robust alternatives, estimator quantiles, boundary occupancy, correlations, optimizer success, and prediction recovery for the released estimands. Original-coordinate values that depend on a structural representative convention are reported as such, not scored as independently recoverable truths.
5. Retain failed datasets in total denominators and show both unconditional failure rates and conditional successful-fit summaries. Predefine any bounded retries.
6. Distinguish Monte Carlo uncertainty on a reported rate from parameter confidence intervals. State that estimator quantiles describe recovery under a fixed generating truth.

### Outputs

synthetic_replicate_manifest.tsv, replicate_fit_index.tsv, parameter_recovery_summary.tsv, prediction_recovery_summary.tsv, failure_summary.tsv, and monte_carlo_precision.tsv.

### Acceptance gate

Independence units and seed separation are verified, all attempts are accounted for, and metric precision is adequate for the stated descriptive claims. A valid conclusion can be poor or biased recovery under the specified design.

### Stop conditions

End the batch at its predeclared dataset/compute limit or precision endpoint. Do not add datasets until a preferred recovery result appears. If numerical failure dominates, stop interpretation, diagnose on development fixtures, and register a fresh confirmation batch after any method change.

## Step 11 - Test additional generating parameter regimes

### Inputs

Step 10 results, the reviewed structural cases and estimands, existing fit/regime inventory, declared scientific contrasts, and an independently budgeted scenario registry.

### Execution steps

1. Select a small number of meaningfully distinct existing-fit regimes: differing inducible error and survival responses, weak-information regions, and relevant boundaries.
2. Freeze the additional truth-selection rule before their recovery outcomes are seen. Record which choices are informed by the exploratory first case.
3. Repeat replay, generator, and deterministic checks for each truth, then the agreed realistic replicate design.
4. Include exact-null cases only with their correct inactive roles. Retain baseline missegregation at the no-induction null and error-size survival at the ploidy-survival null.
5. Compare recoverability across regimes without pooling incompatible estimands or claiming a single favorable truth is representative of the whole domain.

### Outputs

regime_registry.tsv, regime_recovery_summary.tsv, boundary_and_activity_summary.tsv, and generalizability_report.md.

### Acceptance gate

Each regime is feasible and reproducible, has correct parameter activity and applicable structural status, and has a clear conditional interpretation. Revisit special boundary or input cases not covered by the initial structural result. Conclusions identify which regimes were actually tested.

### Stop conditions

Do not launch an unrestricted parameter-space sweep or full Cartesian product. Stop invalid truths or incompatible activity maps. Additional regimes require remaining registered allowance, not favorable outcomes or unused historical budgets.

## Step 12 - Paired data-stream ablation

### Inputs

The same registered truths and synthetic datasets, modality availability, active-parameter rules, and a focused mask registry.

### Execution steps

1. Start with growth alone, growth plus dead-cell data where defined, and all available modalities, as requested in the meeting.
2. Retain the actual available-data configuration, growth plus karyotype plus flow, when dead-cell measurements are unavailable. Treat dead-cell augmentation as a proposed experiment rather than an existing-data result.
3. Add targeted no-karyotype and no-flow contrasts where they answer remaining ambiguities. Avoid automatically multiplying all truths, masks, noise levels, and profiles.
4. Revisit the structural output case for each mask before refitting. A quantity identifiable with all ideal outputs may become structurally confounded after stream removal. Refit the appropriate active parameters/combinations and nuisance quantities, preserving sample identities and legitimate conditioning dependencies.
5. Use paired replicate comparisons to measure changes in parameter, function, and prediction recovery. Do not compare total score heights across different masks as a significance scale.
6. If removed outcomes still condition passage selection, label a conditional information audit. Do not call it leakage-free withheld prediction.

### Outputs

stream_mask_manifest.tsv, mask_parameter_roles.tsv, paired_recovery_comparison.tsv, stream_information_summary.tsv, and measurement_contribution.md.

### Acceptance gate

Every mask is applied to the intended same experimental units, has its structural output case reviewed, receives its own valid refits, and supports a paired information comparison. Missing streams and hypothetical augmentations remain explicit.

### Stop conditions

Hold masks with outcome leakage, changed cohorts, unsupported stream generation, or inactive coordinates left spuriously fitted. Stop at the finite mask/replicate allowance. Lack of improved recovery is a reportable result, not a reason to keep changing weights.

## Step 13 - Evaluate feasible experimental-design changes

### Inputs

The leading recovery limitations, paired ablation results, feasible sample-size/time ranges, and experimental cost constraints.

### Execution steps

1. Use the structural findings and practical recovery limitations to change a small number of candidate factors: karyotype count, live/dead count, measurement time, passage schedule, or starting-ploidy coverage. Revisit the relevant structural case when an input, output, or initial-condition assumption changes.
2. Use illustrative count levels such as 20, 50, and 100 only after checking feasible assay-specific values. Do not equate a flow event count with a karyotype count.
3. Compare targeted designs before combined designs, using paired truth/replicate generation where useful and recording any nested-sampling dependence.
4. Re-simulate dynamics when a schedule or reseeding design changes; do not merely resample an incompatible old trajectory.
5. Rank designs by improvement in the unresolved biological target, recovery reliability, and measurement burden. Include prediction accuracy when parameters remain weak.

### Outputs

design_registry.tsv, design_cost_assumptions.tsv, design_recovery_comparison.tsv, and recommended_measurement_design.md.

### Acceptance gate

Recommendations follow from simulated observations at feasible schedules and sample sizes and state their assumed generator and uncertainty. They identify the ambiguity improved, rather than simply requesting more data.

### Stop conditions

Stop infeasible designs and unregistered combinatorial expansion. A recommendation does not authorize new wet-lab experiments. If no tested design resolves the ambiguity, document that limitation and identify the next discriminating observable.

## Step 14 - Focused profiles, compensation, and biological predictions

### Inputs

Verified synthetic fits, the earlier real-data profile packages with their flags, recovery correlations, and the validated event/function diagnostics.

### Execution steps

1. Retain the biological questions represented by p_misseg and buffer_beta. Select individual-coordinate or combination profiles using the structural results and recovery findings; label redundant-coordinate profiles as diagnostics. Profiling all fourteen across every replicate is not the default.
2. At each imposed target, refit all other active quantities with registered finite search. Exact zero must use the tested zero-capable map and inactive semantics.
3. Link profile behavior to nuisance compensation, boundaries, selected-day changes, and synthetic parameter error. Keep real-data and synthetic-study curves separately identified.
4. Assess recovery of error probabilities, survival contrasts near the configured 2N/4N states, growth/death responses, and viable loss/gain daughter output.
5. Preserve ordinary-division daughter counting, the one-lineage WGD branch, and numerical-versus-biological loss distinctions. Production flux is not the number of unique eventual descendants or proof of ancestry.
6. Report whether similar-score vectors imply stable or divergent biological predictions. A derived profile requires a validated constraint formulation; candidate scatter alone is not one.

### Outputs

focused_profile_registry.tsv, profile_results.tsv, nuisance_paths.tsv, function_recovery.tsv, flux_recovery.tsv, and parameter_vs_prediction_resolution.md.

### Acceptance gate

Each curve has verified candidates and search limitations, every function/flux has a definition and accounting check, and conclusions distinguish individual parameters from combinations and predictions. Existing unresolved real-data profiles retain their status.

### Stop conditions

Stop at the registered target/grid/search allowance and retain unresolved points. Do not use 1.92/3.84 cutoffs or attach confidence labels to uncalibrated scores. Do not force monotonic curves, remove inconvenient boundaries, or expand domains to obtain a desired mechanism result.

## Step 15 - Calibrate formal inference after recovery validation

### Inputs

The accepted structurally informed estimands, a validated generator/estimator and focused profiles, an explicit inferential target, and a separate calibration/coverage simulation allowance.

### Execution steps

1. Specify a coherent generative likelihood or an explicitly simulation-calibrated score procedure for the reviewed individual parameters/combinations. Structural analysis is already an upstream input, not work deferred to this step.
2. Register statistic scaling, nuisance treatment, null-specific parameter activity, and the permitted shapes of confidence sets. Do not promise a unique finite interval for a structurally nonidentifiable original coordinate.
3. Generate independent calibration and coverage batches under the actual sampling design, including relevant nuisance regimes and special boundary cases.
4. Refit required unrestricted and constrained models under the same validated finite-search protocol. Retain optimizer failures and Monte Carlo uncertainty.
5. Validate coverage before reporting formal confidence sets. Keep estimator quantile ranges, Monte Carlo intervals for performance rates, and parameter confidence sets distinct.

### Outputs

inference_calibration_protocol.md, calibration_statistic_registry.tsv, coverage_results.tsv, confidence_set_results.tsv, and inference_limitations.md if this later phase is executed.

### Acceptance gate

Formal parameter or combination intervals require completed statistical calibration and independent coverage validation for the declared model, design, and estimand. Recovery quantiles alone do not satisfy this gate.

### Stop conditions

Defer formal interval claims if structural scope, generative assumptions, optimization adequacy, calibration budget, or coverage precision is inadequate. Stop at the registered finite allowance. Do not repair failed coverage by relabeling old profile scores or removing failed simulations.

## Step 16 - Review tumor structural identifiability, then extend recovery

### Inputs

The reviewed culture structural and practical findings, a working culture generator/refitter with understood limitations, tumor equations and observation manifests, fit packages, scale/initial-state assumptions, and the original in vivo owner's reconciled status.

### Execution steps

1. Enter this new synthetic extension after the culture workflow is validated and its main limitations understood. The gate does not require every culture parameter to be individually recoverable.
2. Preserve the existing tumor study and its outcomes. Reconcile any ongoing approved work through its owner; do not silently restart it or spend its 220/360 budget on synthetic cases.
3. Before designing tumor recovery, perform the tumor-specific structural/input-output review, extending the validated reduction only with explicit resource and retained-material assumptions. Compare known/unknown volume conversion, initial abundance, demand reference, clearance, and latent-resource quantities. Culture structural conclusions do not automatically transfer to tumors.
4. Translate those structural results into tumor estimands and sampling/measurement comparisons, with a tumor structural-to-practical release decision analogous to Step 05. Then define generation at measured burden times and matched harvests using verified volume-error, chromosome-sampling, and histology measurement processes.
5. Preserve animal/section-level sampling structure. Terminal necrotic fraction is not a continuously observed death rate or necessarily a binomial sample of twenty cells.
6. Treat tumor oxygen as latent model dynamics, unlike assigned culture oxygen. Audit rho_2N, fixed N0 and demand reference, retained dead material, initial states, and priors.
7. Begin with one tumor truth and the same deterministic/single-noisy/repeated sequence. Compare data-only and production-penalized estimators separately, refitting each.
8. Prioritize density/initial-size/resource and death/clearance compensation. Any newly free scale or changed volume law is a declared extension. Cross-context sharing/joint mechanism comparisons follow independently validated context results.

### Outputs

invivo_structural_results.tsv, invivo_structural_to_practical_gate.tsv, invivo_synthetic_contract.md, invivo_truth_registry.tsv, invivo_recovery_results.tsv, scale_and_clearance_recovery.tsv, and a scoped cross_context_followup_plan.md.

### Acceptance gate

Tumor-specific structural findings and their bridge have been reviewed before the synthetic design is released. The generator reproduces actual observation units and timing, its free/fixed maps and numerical evidence are valid, and claims distinguish absolute abundance, volume, normalized composition, retained material, and prior influence.

### Stop conditions

Do not launch tumor synthetic recovery with an unreviewed structural/input-output case, unresolved observation/scale contracts, or without a separate authorized synthetic budget. An unresolved structural case requires review; it is not an automatic bypass. Do not import culture point estimates as fixed tumor truth merely to simplify recovery. Failure of scale recovery is a result, not permission for unlimited fitting or silent model changes.

## Dependency and release map

| Work package | Dependencies | Meaning of progression |
| --- | --- | --- |
| Sources, parameter and observation audit: Steps 01-03 | Clarification, source identity, and experimental records | Establish structural inputs and candidate estimands; do not finalize recovery design |
| Culture structural analysis: Step 04 | Explicit model, outputs, inputs, initial conditions, validated tool/analytical route | Deliver reviewed parameter/combination and ambiguity results with exact scope |
| Structural-to-practical bridge and design: Step 05 | Step 04 results and validated full-model mapping | Mandatory release gate for the practical study; unresolved scope returns to review |
| First generator and controls: Steps 06-08 | Accepted Step 05 design and separate execution allowance | Establish numerical and observation consistency for structurally informed estimands |
| Realistic single-case pilot: Step 09 | Generator, deterministic bridge, and released targets | Establish estimator behavior and cost |
| Repeated recovery: Step 10 | Frozen method, budget, dataset manifest | Quantify recovery of the reviewed parameters/combinations |
| Regimes, ablation, design: Steps 11-13 | Valid first-case framework and structural cases for changed outputs/inputs | Avoid a full Cartesian expansion and unsupported transfer of structural status |
| Focused profiles/predictions: Step 14 | Structural estimands and valid candidate fits | Explain remaining practical ambiguity without replacing theoretical analysis |
| Formal calibration: Step 15 | Valid generator/estimator and separate inferential protocol | Later evidence for confidence claims |
| New tumor synthetic study: Step 16 | Culture framework plus tumor-specific structural analysis and bridge | Separate context, review gate, and budget |

## Compute, stopping, and version governance

### New budgets must be derived from measured synthetic cost

This document authorizes no allocations. Writing the plan is not permission to run a symbolic elimination or a numerical benchmark. After a separate implementation/compute instruction, define a finite structural/tool/reduction allowance first. Benchmark practical fitting only after the Step 05 design release and within a separately authorized package.

For a proposed package estimate:

    fit_cost = number_of_truths * number_of_datasets * number_of_masks
               * starts_per_dataset * measured_complete_chain_core_hours
    total_cost = fit_cost + replay/check/generation/report costs
                 + explicitly bounded rescue commitments

Use the actual selected combinations, not a default full product. Count data generation, feasible screening, failed jobs, dependent local polishing, and verification. For illustration only, eighty chains at fifty core-hours each already require 4000 fitting core-hours before checks; the prior culture chain costs must not be assumed to be synthetic costs or treated as available authorization.

Maintain separate ceilings for source/observation engineering, structural cases and reduction validation, the full-model bridge, the first recovery dataset, development replication, confirmation replication, and later regimes/masks/designs. A stalled structural case ends at review, not an unannounced switch to practical fitting. Development batch size around twenty and confirmation size around one hundred are proposals, not claims of sufficient statistical precision or automatic release.

### Finite optimizer policy

Pre-register actual start banks, maximum global/local work per dataset, at most the declared rescue allowance, and a stage cap. Every failed or interrupted attempt remains counted. Whole-dataset quality decisions must use a declared score scale and intended recovery resolution; old real-data score thresholds do not automatically transfer to a new mask or generator.

Keep a lowest verified candidate separate from a checked representative, without transferring evidence between vectors. Preserve real origin identities and original native endings. Numerical reproducibility and practical search agreement remain distinct. Known truth may support a separate diagnostic, but cannot select the reported estimate or count as an independent recovery.

Stopping a valid case with poor parameter recovery is permitted. Stopping an inadequately optimized case requires an optimizer-unresolved label. Do not turn either outcome into a structural proof or continue indefinitely until all parameters look recoverable.

### Source, runtime, and ownership

Reuse the pinned production-v2 numerical model, valid mapping/event checks, and current R/Python/shell orchestration where compatible. First verify exact versions and scientific contracts rather than copying entire legacy controllers. Keep the existing Java/JavaScript/JSON-free core-workflow conventions, English artifacts, typed TSV/CSV, native RDS, and static scientific figures. Symbolic-tool dependencies are documented as part of the initial structural stage; they do not silently change the pinned simulation runtime.

The previously pinned HPC image is:

    /share/lab_crd/taoli/Docker/o2_supply_demand_map_r442_hpc_exact.sif
    SHA256: dfbd7ebe1209814ca25358e56ae4d5553c6dc19852f5567ce369ab23bebe3372

At future execution, verify that image and source identity, scheduler limits, red/xxlarge allocation, no-requeue policy, and worker/thread isolation. Preserve Slurm-managed concurrency rather than inferring a running-core throttle from a total budget. None of these operations is performed by writing this document.

Proposed new code and result locations, to be created only during authorized implementation:

    oxygen/code/O2_supply_demand_MAP/analysis/structural_identifiability/structural_first/
    oxygen/results/structural_identifiability/YYYYMMDD_culture_symbolic_<study_id>/
    oxygen/code/O2_supply_demand_MAP/analysis/synthetic_identifiability/
    oxygen/results/synthetic_identifiability/YYYYMMDD_culture_recovery_<study_id>/
    oxygen/results/synthetic_identifiability/YYYYMMDD_invivo_recovery_<study_id>/

Use isolated immutable releases, input/source hashes, append-only candidate histories, atomic state/receipt handling, and a sole dispatcher per run. Do not overwrite old source snapshots, closed culture results, manuscript figures, or the other context's state. Store canonical large payloads on approved HPC storage with verified local collection appropriate to available space.

## Deliverable and interpretation standard

The main report must first state the structural findings and their exact assumptions, then answer which reviewed parameters/combinations are practically recoverable, which weak quantities leave predictions stable, and which uncertainties change biological conclusions. It must also identify informative measurements and tested feasible design improvements. A reduced-model result and a full-simulator conclusion must never be merged without a justified bridge.

Required report components:

- Source/design/parameter contract and actual tested scenario coverage.
- Structural equations, case-specific proof/status, identifiable combinations, reduction limitations, and the accepted structural-to-practical mapping and design gate.
- True-versus-estimated parameter plots and bias/error/quantile tables by truth and mask.
- Optimizer success, failure, boundaries, and per-start quality separate from sampling variability.
- Correlation/compensation summaries, with transformed and natural units stated.
- Recovery of growth, full distributions, relevant survival/error functions, and audited daughter-output predictions.
- Paired ablation and design comparisons, including hypothetical versus existing observations.
- Clearly labeled focused profiles, structural-tool findings, and any later calibration evidence.
- Allocated actual/committed costs, all finite stops, unavailable evidence, and reproducibility/collection status.

| Observed result | Permitted interpretation |
| --- | --- |
| Expected observations cannot be reproduced even at the frozen truth | Generator/adapter/model-contract failure; do not interpret recovery |
| Truth replay works but blind fits cannot approach attainable predictions | Optimizer adequacy remains unresolved under the tested effort |
| Different vectors reproduce deterministic observations | Candidate parameter ambiguity; investigate combinations and structural assumptions |
| Valid repeated noisy fits have broad or biased estimates | Limited estimator recovery under the specified generator/design; examine noise, nuisance, bounds, and loss effects |
| Parameters vary but relevant functions/predictions recover well | Emphasize prediction or combination resolution within tested regimes |
| Similar scores imply different biological predictions | The tested observations do not resolve that conclusion |
| Added observations improve paired recovery | Evidence for measurement value conditional on the tested generator and estimator |
| Reduced symbolic model has an identifiable combination | A result for that reduction and ideal observation assumptions, not automatically the full simulator |

Estimator percentile ranges, Monte Carlo intervals for performance rates, and confidence sets for parameters are different outputs. Label each correctly. No formal p value, confidence threshold, or global-identifiability claim follows merely from a low score or a finite recovery test.

## Appendix A - Unresolved inputs that must remain visible

| Input or decision | Required resolution | Work that need not wait |
| --- | --- | --- |
| Exact Mathematica paper/package | Obtain citation/name and investigate authoritative compatibility/licensing information; use a validated alternative structural route if needed | Source, parameter, and assay audits; practical design still requires structural review |
| Unresolved structural case or reduced/full-model mapping | Record the exact unresolved scope and review a tractable route or explicitly limited follow-up | Evidence preservation and independent input audits; no automatic synthetic-study release |
| Sample-specific karyotype sizes | Derive from physical cell lists and sample manifest | Source/assay audit and ideal-output structural cases with explicit assumptions |
| Live/dead availability and assay mapping | Confirm counts, timing, denominator, classification, and model-state correspondence | Existing-output structural cases; additional dead-cell outputs stay explicitly hypothetical |
| Flow acquisition/gating/event counts | Recover metadata or declare an approximation and sensitivity range | Non-flow structural cases and metadata audit |
| Passage timing versus target-count conditioning | Reconcile actual protocol and fitting selection rule | Source/parameter inventory; not full real-design recovery claims |
| Appropriate error/dependence model | Use available replication and assay definitions; label unsupported assumptions | Structural analysis of the declared ideal outputs and assay documentation |
| Structural and synthetic compute ceilings | Register bounded structural cases first; measure fitting cost only after design release | Documentation and existing-result inspection |

Questions for data owners or collaborators are future workflow items. This documentation request does not authorize sending messages or requesting package access on the user's behalf.

## Appendix B - Implementation handoff boundary

The user requested this documentation revision in response to Noemi's structural-first clarification. The current action is **documentation only**. No code, controller, numerical artifact, remote state, or job has been changed by preparing or revising it.

For a future implementation handoff, read this plan together with Noemi's later clarification, the original profile/mechanism plans, and current ownership logs. Complete Steps 01-03, then the Step 04 structural analysis and Step 05 structural-to-practical review before releasing the synthetic design. Register a separately authorized finite structural package first; practical compute is a later release. Do not resume real-data global-optimum searches as a prerequisite, reinterpret old unresolved profiles as completed identifiability evidence, or spend old budgets on the new program.

## Appendix C - Step-number migration and existing-work transition

This revision has sixteen steps. Existing files, receipts, run IDs, and approvals using the original fourteen-step numbering retain that historical meaning. Do not rename an old run or infer that an existing Steps 04-06 artifact is evidence for the newly numbered structural Steps 04-05.

| Original step | Revised step | Change |
| --- | --- | --- |
| 01 Sources/ownership | 01 | Add structural-first precedence and ownership |
| 02 Parameters/targets | 02 and 05 | Audit first; finalize estimands after structural findings |
| 03 Observation contract | 03 and 05 | Inventory actual design first; release recovery design after structural review |
| 13 Structural-tool component | 04 | Move to the initial scientific analysis |
| No separate bridge step | 05 | Add mandatory structural-to-full-model and practical-design review |
| 04 Generating truth | 06 | Move after the structural gate |
| 05 Finite observations | 07 | Use the reviewed synthetic design |
| 06 Deterministic recovery | 08 | Check the reviewed parameters/combinations |
| 07 Single noisy dataset | 09 | Fit released coordinates and nuisances |
| 08 Repeated datasets | 10 | Quantify practical recovery after structural review |
| 09 Additional regimes | 11 | Check structural scope for each regime |
| 10 Stream ablation | 12 | Revisit structural output cases |
| 11 Experimental design | 13 | Combine structural and practical findings |
| 12 Profiles/functions | 14 | Interpret remaining practical compensation |
| 13 Calibration component | 15 | Keep formal confidence inference late |
| 14 Tumor extension | 16 | Require tumor-specific structural review before recovery design |

Before another practical-stage release, the execution owner must reconcile any already authorized or completed work against this revision. Preserve valid deterministic, generator, and optimizer artifacts as engineering or historical diagnostic evidence. They do not retroactively satisfy the new structural gate and must not be relabeled as a structurally informed confirmatory design.

This documentation update does not cancel or pause live allocations, change an existing controller, or authorize a new one. Any required operational transition remains with the original execution owner under the applicable user instruction. The next proposed stage and its budget must cite the revised step numbers and the status of the structural review.

Originally prepared on 2026-11-12; scope and structural-first sequencing revised on 2026-11-13, America/New_York. No commit, push, HPC task, external message, or code implementation is included in this document revision.
