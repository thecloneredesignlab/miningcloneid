# Integrating the completed Sobol reference into the culture pilot

Reviewed on 2026-11-09. This is a source and diagnostic-plan supplement to the [profile plan](profile_likelihood_implementation_plan.md) and [mechanism-identifiability plan](invitro_mechanism_identifiability_implementation_plan.md). It reports no new culture-profile conclusions or symbolic identifiability results.

## Decision

Use the completed **global Sobol N131072** results to prioritize diagnostics, connect parameters to mechanisms, and formulate hypotheses for later measurement comparisons. Keep the approved production-v2 culture model, objective, extended parameter domain, initialization policy, optimization gates and HPC configuration unchanged. No new fitting jobs are required by this integration. Existing authorized culture work retains its own execution and acceptance gates.

The directory name `soft_couping_eFAST` is historical. The high-resolution reference is Sobol, with five independent scramble repeats. Early global eFAST results are not treated as a stable cross-method validation. No Sobol/eFAST recomputation or archived neighborhood continuation is included.

## Evidence and source authority

The reference lives under `/Users/4482173/Documents/GitHub/soft_couping_eFAST/oxygen/results/eFAST/sobol/global_convergence_20261021/comparison_N65536_N131072/`.

The bridge audit checked all 47,838 comparison-table cells: two outputs, 201 oxygen values, 14 S1, 14 ST and all 91 unique S2 pairs per output and oxygen. All N131072 cells pass the recorded operational criteria: five-repeat range at most 0.05 and adjacent-resolution mean difference at most 0.02. It independently reconstructed means and ranges for 20,100 N131072 cells from the supplied five-repeat source table and checked ST oxygen-band integral averages. That repeat table covers ST and 36 focal S2 pairs; it does not cover every S1 or S2. The original full-export validation is retained as historical evidence rather than represented as newly repeated here.

The compact `parameter_and_pair_summary.tsv` also contains only ST and those 36 S2 pairs, involving `p_mis_base`, `mu_hp` or `gamma_mu`. Other pairs below were extracted from `all_index_comparison.tsv.gz`.

An important provenance distinction: `inputs/source_audit.json` describes neighborhood-derived preparation and has `source_root` ending in `neighborhood10pct`. It is not sufficient evidence for the global bounds. For the completed global result, the frozen global plan, its global case and the case's hashed `parameter_bounds.tsv` are the authority. The latter exactly matches the supplied `inputs/global_bounds.tsv`; its SHA-256 is `8d2ea8ef048459f0570e56ef1fb2434900f93e40708c7f132c7e7788ef9f2b97`. No neighborhood result was resumed or used as a global result.

## Compatibility with production v2

| Aspect | Old global Sobol | Current culture pilot | Integration rule |
| --- | --- | --- | --- |
| Mechanistic core | Chromosome grid 22–154, N_UNIT 22, oxygen-dependent growth and ploidy-dependent death | Same core R/C++ file hashes and corresponding switches in frozen production v2 | Reuse audited parameter meanings and function connections |
| Protocol | One externally fixed oxygen value for each matrix calculation, 0–5% in 0.025% steps | Both initial-ploidy cohorts, assigned oxygen histories, finite passages and inherited states | Do not substitute fixed-oxygen eigenanalysis for lineage simulation |
| Output | Leading-eigenvalue net growth and leading-eigenvector mean ploidy | Growth and complete chromosome/flow distributions at selected passage endpoints | Evaluate compensation using actual culture predictions |
| Active quantities | 14 independently sampled mechanism parameters | 20 unrestricted fitted coordinates: the same 14 plus four initial-distribution and two observation-error parameters | Keep all culture nuisance coordinates available |
| Sampling or inference | Uniform/log-uniform variance decomposition over declared independent bounds | Reoptimization of a weighted average data score, with protocol-failure barriers | Sampling distributions do not become priors; sensitivity indices do not become profile scores |
| Exact mechanism nulls | Neither target includes zero | Exact zero admitted for both targets; nuisance dimensions adjusted | Use the new unrestricted extended-domain baseline and its constrained fits |
| Nonlinear dynamics | Fixed live-state linear operator M = G − diag(mu) | Finite-time simulator, configured crowding and passage selection/propagation | Core identity does not establish protocol or objective identity |

The current worktree, old worktree and frozen production-v2 core C++ SHA-256 all equal `868f4abd6f62e0fdcb3cae5bf4c32985fb216cfa83fc0d80f8db648f01435bf4`; the core R files also match. Supporting code does not uniformly match. The old `common_semantics` matches frozen v2 but differs from the current worktree; the old culture objective and backend differ from frozen v2. The pilot adapter therefore continues to use the explicitly frozen and replay-validated production backend. No old runner or sampler is copied into it.

Material parameter-domain differences include:

| Parameter | Old Sobol domain | Culture extended domain |
| --- | --- | --- |
| `lam_max` | 0.005–0.3 | 0.005–1.5 |
| `p_mis_base` | 0.000001–0.1 | 0.000001–0.005 |
| `p_misseg` | 0.005–0.5 | **0**–0.4 |
| `k_o_mis` | 0.001–0.02 | 0.00001–0.02 |
| `buffer_beta` | 0.01–10 | **0**–10 |
| `O2_crit` | 0.1–2.5 | 0.0001–2.5 |
| `n_O` | 0.5–5 | 0.0001–5 |

The bridge exports a comparison for all 20 culture coordinates, including transforms. Optimization transforms are not statements about parameter probability distributions. Bounds remain explicit assumptions; low old ST does not justify fixing any nuisance coordinate.

## Diagnostics to carry forward

1. **Inducible missegregation.** Track `mu_hp`, `gamma_mu`, `k_o_mis`, `p_mis_base`, `O2_crit` and `n_O` along the `p_misseg` profile. At matched oxygen and chromosome count, compare death hazard, saturation `mu/(mu+k_o_mis)`, the induced increment and the total per-chromosome error probability. In the induction formula, scaling `mu_hp` and `k_o_mis` together preserves that saturation at fixed oxygen-shape parameters; direct death still changes, so this is not an invariance of the full model. Verify any compensation against complete culture predictions and accounted daughter output.
2. **Ploidy-dependent survival.** Track `buffer_smax`, `buffer_beta` and `buffer_n_exp` jointly, and compare `s_N = buffer_smax * exp(-buffer_beta*(44/N)^buffer_n_exp)` and `S_N(m) = s_N^m` for positive error size, with `S_N(0)=1`. Ask whether survival functions or effective viable output remain stable while individual coordinates move. Small old pairwise S2 does not exclude compensation in culture.
3. **Null interpretation.** At `p_misseg=0`, retain basal errors, WGD, growth/death and survival; `k_o_mis` is data-inactive. At `buffer_beta=0`, refit `buffer_smax` and retain error-size dependence; `buffer_n_exp` is data-inactive. Canonical inactive values are not estimates. The production score has no parameter prior contribution that would reactivate these nuisance coordinates.
4. **Observed and effective output.** Retain complete growth, chromosome and flow predictions, selected passage days, all nuisance paths, boundaries and optimizer flags. Interpret probability, expected errors per division, viable daughter production, biological death and numerical losses with their own denominators and verified event ledger. Parameter paths and function traces are diagnostics of parameter-constrained fits, not separate profiles of derived quantities.
5. **Later measurements and COMBOS.** Use oxygen response regions as candidate conditions only. Select later measurements by finite-time prediction separation between good culture fits under realistic observation noise and protocols. Targeted stream removal requires a separate unrestricted baseline and nuisance refits for each mask. Reuse the mechanism map when preparing the five-live-class reduced symbolic model, then validate reduction error and assumptions; the map itself is not a structural proof.

Verified old mean-ploidy references are `ST(p_misseg)` = 0.510525 over 0–1% and 0.316903 over 3–5%, and `ST(buffer_beta)` = 0.107296 and 0.153382 over those bands. These are interval-integral averages, not averages over a measured exposure probability distribution. The largest mean-ploidy S2 for `p_misseg × mu_hp` is 0.0566468 at 0% oxygen; for `p_misseg × gamma_mu` it is 0.0228533 at 0.525%. These values motivate tracking, without establishing direction or compensation.

## Interpretation boundaries and artifacts

The old collection and comparison scripts provide useful engineering patterns: explicit file allowlists, staged transfer with source/destination SHA-256 verification, preservation of existing artifacts, complete output/oxygen/parameter-pair keys, and independent reconstruction of summary statistics. The current pilot already applies corresponding source, result, checkpoint and scoring checks. Its culture fitter and scheduler remain the current implementation. The old comparison script includes a Chinese report generator, so it is referenced as historical source rather than copied or executed to generate new outputs; this bridge and all current reports use English.

S1/ST measure output variance under the old domain and sampling assumptions. S2 measures nonadditive interaction without its direction. Neither establishes identifiability or compensation. ST includes overlapping interaction contributions and cannot be summed as disjoint variance shares. Operational convergence and repeat ranges are not significance tests or parameter confidence intervals. Old optimization endpoints are not posterior samples or culture baselines. The current weighted score still has no calibrated 1.92/3.84 thresholds, likelihood-ratio p values or formal confidence intervals.

Reproduce the bridge with `oxygen/code/O2_supply_demand_MAP/analysis/structural_identifiability/audit_sobol_culture_bridge.py`, passing the current project, old project and pilot study paths. It uses only completed tables and source files, with zero model calls and zero job submissions. Under the current pilot's `sobol_bridge/`, it writes English `audit.json` with source hashes, `code_comparison.tsv`, `parameter_domain_comparison.tsv`, `target_ST_reference.tsv`, `candidate_pair_reference.tsv` and `profile_diagnostic_hypotheses.tsv`. External historical source files, including the Chinese interpretation, remain unchanged.
