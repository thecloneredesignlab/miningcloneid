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

## Replay on archived inputs

Historical scripts use each study directory as the working location for both the source and its recorded observation design. [replay_layout.tsv](replay_layout.tsv) maps the reorganized files to those study paths. First create a separate working copy with the local `oxygen/results` inputs. The experimental count replay also uses the original CloneID snapshot under the sibling `soft_couping_org` repository. In that separate parent directory, provide the same relative snapshot path.

Run the following from the source repository, substituting the path of that working copy.

```sh
python3 oxygen/code/O2_supply_demand_MAP/analysis/structural_identifiability/structural_first/replay.py \
  --workspace /path/to/isolated-checkout --rebase-count-provenance
```

The staging command restores the historical source paths in the isolated copy and checks archived plan inputs. From the isolated root, run the entry scripts at their study paths listed in the layout. For example, use `Rscript` with `.../20261118_culture_experimental_OLS_step02/reproduce_ols.R` and a scratch output directory. For the donor inheritance result use `python3` with `.../20261118_culture_experimental_inheritance_01/rank_cases.py`. Use separate processes for cases that import different versions of the numerical modules. Python rank verification requires NumPy and python-flint 0.8.0, and the C1 replay imports SciPy through the fitting code.

The local documentation archive is at `oxygen/results/analysis_documentation/`. Its `index.tsv` maps former `docs/` locations to exact archived documents. For a new preparation, invoke the relevant versioned preparation entry with the repository root containing those inputs.
