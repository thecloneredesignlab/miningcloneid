# Structural and practical identifiability workflows

Start with the production or experimental passage workflow for the reported in vitro structural results. Earlier model and certificate directories record the supporting analyses and their assumptions.

| Directory | Question and entries | Local evidence |
| --- | --- | --- |
| [01 Model and tools](in_vitro/01_model_and_tools/) | Whether the reduced model and available methods are suitable. See the initial reduction, reduction validation, Stage A review, and public tool checks. | The culture symbolic, reduction revision, Stage A, and public tools run directories. |
| [02 Structural certificates](in_vitro/02_structural_certificates/) | Finite-state sensitivity ranks and independent checks. See `finite_state/case_worker.py`, `boundary_review/`, `controls/`, `initial_observer/`, and the certificate reviews. | The Stage B, C04, four controls, and observer run directories. |
| [03 Full model bridge](in_vitro/03_full_model_bridge/) | Step 05 prepares a correspondence between structural findings and the full production model. Entry `prepare.py` takes the repository root and a new study directory. | `20261116_culture_step05_bridge_design_01` |
| [04 Production schedule](in_vitro/04_production_schedule/) | Production passage crosswalk, model-selected endpoint growth, the karyotype observer, and the global-boundary controls. Entries include `local_rank/certify.py`, `local_rank/final_checks.py`, and `independent_rank/verify_rank.py`. | `20261117_culture_actual_schedule_step01` through `step03`, and `20261117_culture_independent_rank_verification_01`. |
| [05 Experimental passages](in_vitro/05_experimental_passages/) | Reproduction of 114 slopes from 335 counts, unknown versus fixed seeding, donor inheritance, and the direct 221-count rank check. Entries include `ols_observation/reproduce_ols.R`, `unknown_seeding/verify_results.py`, `fixed_seeding/verify_fixed.py`, `donor_inheritance/rank_cases.py`, and `direct_count_rank/summarize_stage1.py`. | The `20261118_culture_experimental_*` run directories. |
| [In vivo Step 16](in_vivo/step16_metadata_audit/) | Tumor metadata and source audit. Entries `capture_authority.py`, `export_contract.R`, and `verify_audit.py`. | `20261113_invivo_structural_01`; authority comes from the frozen continuation archive. |

The [synthetic and practical workflows](../../synthetic_identifiability/README.md) contain D0, D19, N20 and the experimental precision analyses. The C1 rank replay reads the saved count fit and profile points from that pipeline.

## Experimental passages from 335 recorded counts to two observation maps

The 335 recorded `correctedCount` values comprise 114 day-0 seeding counts, 107 intermediate counts, and 114 final counts. The intermediate and final counts give the **221 post-seeding counts**. Each of the 114 observed passages contributes one **experimental OLS slope**, fitted with an intercept to the log counts from that passage, including day 0. The two analyses therefore use the same records in different ways. The C1 count model uses the 114 day-0 values to set the numbers of live cells at seeding and uses the 221 later values as individual observations.

Paths below are relative to the repository root. The dated result directories are local archives under `oxygen/results`.

### 1. Audit and align the 335 recorded counts

- **Script** [count_audit/audit_counts.py](in_vitro/05_experimental_passages/count_audit/audit_counts.py)
- **Reads** The CloneID snapshot `../soft_couping_org/data/InVitroData/cloneid_passaging_sum159_snapshot_20260731.tsv`, the prepared `count_growth_observations.tsv` and dead-cell annotations `sum159_dead_cell_db_enriched_20260731.tsv`, the recorded observations in `data/InVitroData_LTEE/invitro_passage_observations.tsv`, and `oxygen/data/metadata.csv`. It joins these to `growth_raw_count_manifest.tsv` and `experimental_schedule.tsv` in the archived `synthetic_identifiability/20261113_culture_contract_steps01_03`, the production `passage_crosswalk.tsv` in `structural_identifiability/20261117_culture_actual_schedule_step01`, and the lineage, counts, and fit metadata in `structural_identifiability/20261109_seed144_wp0/source/seed144`.
- **Writes** `count_event_mapping.tsv` and `audit_checks.tsv` in `oxygen/results/structural_identifiability/20261118_culture_experimental_OLS_step01/`. The mapping identifies each count by record ID, passage, calendar day, and model time step.

### 2. Reproduce the 114 slopes and build the model observation map

- **Scripts** [reproduce_ols.R](in_vitro/05_experimental_passages/ols_observation/reproduce_ols.R) reads the audit mapping, raw CloneID snapshot, and archived experimental slopes. For each passage it fits `lm(log(pmax(1, correctedCount)) ~ days)`, with days measured from seeding.
- [build_calendar_map.R](in_vitro/05_experimental_passages/ols_observation/build_calendar_map.R) reads the mapping, newly reproduced slopes, and frozen seed144 `fit_summary.tsv`. It uses [observation_operator.R](in_vitro/05_experimental_passages/ols_observation/observation_operator.R) to turn the recorded days into OLS weights and align them with model time steps. [build_passage_followup.R](in_vitro/05_experimental_passages/ols_observation/build_passage_followup.R) reads the production crosswalk and raw snapshot to trace candidate donor links.
- **Writes** `experimental_ols_reproduction.tsv`, `ols_calendar_map.tsv`, and `passage_gap_followup.tsv` in `oxygen/results/structural_identifiability/20261118_culture_experimental_OLS_step02/`. The 114 reproduced slopes agree with the recorded slopes to a maximum absolute difference of about `5.55e-15/day`.

### 3. Test independently initialized karyotypes with unknown and fixed seeding quantities

| Case | Scripts and inputs | Local results |
| --- | --- | --- |
| A. Unknown seeding quantities | [unknown_seeding/certify_experimental.py](in_vitro/05_experimental_passages/unknown_seeding/certify_experimental.py) uses the count mapping, OLS calendar map, passage follow-up, and the production rank archive and seed144 parameters. It allows independent starting karyotype compositions at 14 segments and lets the seeding quantities vary. [verify_results.py](in_vitro/05_experimental_passages/unknown_seeding/verify_results.py) checks the saved certificates. | `oxygen/results/structural_identifiability/20261118_culture_experimental_OLS_rank_01/`. In each visibility specification the OLS rank is 114 before and after adding the 14 mechanism columns. The all-count rank is 335 before and after. |
| B. Fixed seeding quantities | [fixed_seeding/certify_fixed.py](in_vitro/05_experimental_passages/fixed_seeding/certify_fixed.py) reuses the same design and the saved A witnesses and fixes the live-cell seeding quantities to the recorded values. The same 14 starting karyotype compositions can still vary independently. [verify_fixed.py](in_vitro/05_experimental_passages/fixed_seeding/verify_fixed.py) checks the saved certificates. | `oxygen/results/structural_identifiability/20261118_culture_experimental_fixed_inoculum_01/`. The OLS rank remains 114 before and after adding the mechanism columns. The count rank remains 221 because the 114 day-0 rows are constant. |

Both cases report **0/14 locally identifiable mechanism parameters**. In A, changes in seeding quantities and starting karyotypes can exactly compensate mechanism changes. In B, changes in the starting karyotypes alone suffice. In each archive, open `ols_shared/rank.tsv` for the nuisance-only and joint ranks, `parameter_verdicts.tsv` for the mechanism results, and `verification_summary.tsv` for the replay checks.

### 4. Restore 12 donor links and certify the OLS ranks in C1–C8

- **Scripts and inputs** [donor_inheritance/inherit.py](in_vitro/05_experimental_passages/donor_inheritance/inherit.py) reads the A/B passage design, B's `ols_shared/cut_edges.tsv`, the seed144-derived mechanism witness from the production rank archive, and the frozen `invitro_observed_kary.tsv`. It restores karyotype inheritance from the recorded donor at 12 transfers and uses independent karyotype proxies for the two earlier `dp` starts. Run it with `--root hist` and `--root gaussian` to generate the two sets of derivative tables. [rank_cases.py](in_vitro/05_experimental_passages/donor_inheritance/rank_cases.py) reads those tables and certifies full-rank minors. [verify_independent.py](in_vitro/05_experimental_passages/donor_inheritance/verify_independent.py) checks sensitivities with a separate implementation.
- **Writes** The design and sensitivities under `hist/` and `gaussian/`, plus `rank_summary.tsv`, `rank_minors.tsv`, `parameter_summary.tsv`, and `verification.tsv` in `oxygen/results/structural_identifiability/20261118_culture_experimental_inheritance_01/`.
- **Results** The eight specifications cross two starting-proxy models, two donor-transfer rules, and two dead-cell visibility rules. The full OLS column ranks in C1–C8 order are **16/16, 18/18, 40/40, 42/42, 20/20, 22/22, 44/44, and 46/46**. All 14 mechanism parameters are locally identifiable at the tested point in each specification. C1 uses fixed empirical karyotype proxies, exact transfer of donor karyotype proportions, and two dead-cell visibility fractions `eta_H` and `eta_M` shared by 2N and 4N.

### 5. Fit C1 to the 221 post-seeding counts

- **Inputs and setup** [count_model.py](../../synthetic_identifiability/experimental/counts/count_model.py) reads the same raw snapshot, audit mapping, and calendar map, plus the C1 design `hist/passage_design.tsv` and `hist/dp_proxy_histograms.tsv`. The compiled [count_kernel.cpp](../../synthetic_identifiability/experimental/counts/count_kernel.cpp) evaluates the log-count objective. The 114 fixed day-0 records initialize the model and the 221 positive-day records contribute residuals.
- **Fit and profile** [fit_counts.py](../../synthetic_identifiability/experimental/counts/fit_counts.py) in `real` mode starts from `main_constrained_best.tsv`, `main_constrained_refinement.tsv`, and `profile_reference_refits.tsv` in `oxygen/results/synthetic_identifiability/20261118_experimental_OLS_practical_01/`. It optimizes the 14 mechanism parameters and two shared visibility fractions jointly, subject to the parameter bounds and donor supply constraints. [profile_counts.py](../../synthetic_identifiability/experimental/counts/profile_counts.py) profiles parameters by fixing one coordinate and reoptimizing the others. [edge_counts.py](../../synthetic_identifiability/experimental/counts/edge_counts.py) and [cross_profiles.py](../../synthetic_identifiability/experimental/counts/cross_profiles.py) generate additional profile points.
- **Writes** `count_inputs.tsv`, `input_audit.tsv`, `real_fits.tsv`, `real_best.tsv`, `real_residuals.tsv`, and the profile tables, including `profiles.tsv` and `profile_intervals.tsv`, in `oxygen/results/synthetic_identifiability/20261118_experimental_counts_practical_01/`. The next structural step uses the fit and selected profile parameter vectors as points at which to check rank.

### 6. Certify the 221-row C1 count Jacobian directly

- **Select points** [direct_count_rank/prepare_points.py](in_vitro/05_experimental_passages/direct_count_rank/prepare_points.py) reads the original OLS witness in `hist/mechanism_witness.tsv`, the C1 `count_inputs.tsv`, `real_best.tsv`, and the five saved profile support tables `profiles_{direct,supplement,edges,cross,polish}.tsv`. It writes `exact_points.tsv`, `point_selection.tsv`, and `source_provenance.tsv`. The selection includes the original OLS witness, a strictly feasible point near the best count fit, two checks at and near the `mu_hp` lower bound, and profile points.
- **Compute and verify** [certify_counts.py](in_vitro/05_experimental_passages/direct_count_rank/certify_counts.py) reads the point tables, the audited count inputs, the donor design, and the saved C1 sensitivities used to check agreement at the OLS witness. It propagates interval derivatives for the 221 post-seeding log counts and certifies an invertible 16-by-16 minor. [summarize_stage1.py](in_vitro/05_experimental_passages/direct_count_rank/summarize_stage1.py) verifies the source checksums and independently checks the saved minors.
- **Writes** `observation_audit.tsv`, `certificates.tsv`, `certified_minors.tsv`, `certificate_validation.tsv`, and `REPORT.md` in `oxygen/results/structural_identifiability/20261118_culture_experimental_counts_local_01/`. The **221-by-16 Jacobian has rank 16 at all 18 feasible points**. The archived numerical best fit is a separate analytic check because it slightly violates one donor supply constraint.

With the day-0 counts fixed, a slope is a constant plus a linear combination of the post-seeding log counts with fixed weights. Consequently `J_OLS = W J_log_counts`, and full OLS column rank implies full count column rank for the same specification and parameter point. This gives the count ranks for C2–C8 at their OLS test points. The direct C1 calculation extends the structural check to the additional feasible points selected around the count fit and profiles.

## Where to start reproducing the results

### Replay the archived evidence

Start here if the local result archives are available. Set `workspace` to an existing isolated copy of the repository containing the local `oxygen/results` inputs. Provide the original CloneID snapshot at `../soft_couping_org/data/InVitroData/cloneid_passaging_sum159_snapshot_20260731.tsv` relative to that copy. Use a Python environment with NumPy, SciPy, and python-flint 0.8.0, and make `Rscript` available.

Run this from the source repository. [replay.py](replay.py) stages sources at the historical paths recorded in [replay_layout.tsv](replay_layout.tsv) and checks and rebases the C1 input paths inside the isolated copy. The historical scripts use their own location to find inputs and outputs.

```sh
workspace=/path/to/isolated-checkout
python3 oxygen/code/O2_supply_demand_MAP/analysis/structural_identifiability/structural_first/replay.py \
  --workspace "$workspace" --rebase-count-provenance
```

From the staged copy, rebuild the OLS observations from the saved audit mapping and verify the rank evidence. Run the direct C1 verifier before regenerating the C1–C8 summary because the C1 provenance records that summary's checksum.

```sh
cd "$workspace"
structural_results=oxygen/results/structural_identifiability
ols_run="$structural_results/20261118_culture_experimental_OLS_step02"
scratch_ols="$PWD/verification/ols"
mkdir -p "$scratch_ols"
Rscript "$ols_run/reproduce_ols.R" "$scratch_ols"
Rscript "$ols_run/build_calendar_map.R" "$scratch_ols"
Rscript "$ols_run/build_passage_followup.R" "$scratch_ols"

python3 "$structural_results/20261118_culture_experimental_OLS_rank_01/verify_results.py"
python3 "$structural_results/20261118_culture_experimental_fixed_inoculum_01/verify_fixed.py"
python3 "$structural_results/20261118_culture_experimental_counts_local_01/summarize_stage1.py"
python3 "$structural_results/20261118_culture_experimental_inheritance_01/rank_cases.py"
```

The R scripts write fresh tables to `verification/ols/`. The Python scripts read the saved inputs in the isolated result directories and verify or regenerate the rank evidence there.

### Rebuild from an earlier stage

- **Start at the original records.** Run the staged `audit_counts.py` from Step 1 with access to the contract and frozen seed144 inputs. Use its `--snapshot`, `--prepared-counts`, and `--dead-annotations` options if those source files are in different locations. Follow with the three R scripts in Step 2. Use a separate staged copy for regenerated inputs so the saved evidence still has its original provenance.
- **Recompute the OLS certificates.** Start with `certify_experimental.py --phase ols --visibility shared --out <A-run>/ols_shared` and repeat with `split`. For A's count phase, pass the corresponding OLS `witness_inocula.tsv` via `--witness`. Case B reads the corresponding A runs and accepts the same phase, visibility, and output options. For C, run `inherit.py --root hist --out <C-run>/hist` and its `gaussian` equivalent before `rank_cases.py`. Here `<A-run>` and `<C-run>` represent the staged result directories above.
- **Recompute the direct C1 Jacobian.** Use the archived `exact_points.tsv` and `point_selection.tsv`. These commands are for macOS and start at the root of the staged copy. The certification writes new minors to that copy's C1 archive.

```sh
count_run=oxygen/results/synthetic_identifiability/20261118_experimental_counts_practical_01
c1_run=oxygen/results/structural_identifiability/20261118_culture_experimental_counts_local_01
mkdir -p verification
export CULTURE_COUNT_LIB="$PWD/verification/culture_count_practical.dylib"
c++ -std=c++17 -O3 -dynamiclib "$count_run/count_kernel.cpp" -o "$CULTURE_COUNT_LIB"
python3 "$c1_run/certify_counts.py" --bits 384
python3 "$c1_run/summarize_stage1.py"
```

To select new points from a fresh count fit, first run the staged `count_model.py`, `fit_counts.py real`, and the profile scripts in Step 5. Use fresh fit and profile outputs, since the profile scripts can reuse saved work. Then run `prepare_points.py` to write the new point tables and provenance before certifying. The refit starts from the archived OLS fits, and the final report reads the practical comparison table `parameter_comparison.tsv` alongside the structural certificates.

The local documentation archive is at `oxygen/results/analysis_documentation/`. Its `index.tsv` maps former `docs/` locations to exact archived documents. For a new preparation, invoke the relevant versioned preparation entry with the repository root containing those inputs.
