# Remove Migrated PKPD Live/Dead Model From miningcloneid

Date: 2026-06-29

## Goal

Remove the migrated gemcitabine PK/PD live/dead in-vitro fitting and visualization code from this miningcloneid repository after confirming the canonical copy exists in `Gemcitabine-model`.

The cleanup should be active on these branches:

- `pujan_branch_edits`
- `doxoFits`
- `soft_coupling`
- `feat-O2-G-resource`
- `main`

The removal must not delete shared dependencies or data used by other miningcloneid analyses.

## Canonical Destination To Verify First

Before deleting anything from miningcloneid, verify the migrated copy from the committed `Gemcitabine-model` tree, not from untracked files in its working tree. As of this review, both `HEAD` and `origin/main` in `/Users/4470246/Repositories/Gemcitabine-model` resolve to:

```text
d8e85a7407378e07ad929b5e0009c6e7644d9256
```

Use `git ls-tree` and `git rev-parse <ref>:<path>` against the chosen target commit. Do not rely on `find` output from a dirty working tree or on untracked timestamped validation outputs.

Verify the migrated copy in:

- `/Users/4470246/Repositories/Gemcitabine-model/Code/in-vitro/pkpd_live_dead_model/`
- `/Users/4470246/Repositories/Gemcitabine-model/Data/in-vitro/pkpd_live_dead_model/`
- `/Users/4470246/Repositories/Gemcitabine-model/Data/GemDelayKillTerm/processed/`
- `/Users/4470246/Repositories/Gemcitabine-model/docs/pkpd_live_dead_model_fitting_migration_status.md`

Minimum expected target files:

- `Code/in-vitro/pkpd_live_dead_model/plot_invitro_fit_outputs.py`
- `Code/in-vitro/pkpd_live_dead_model/run_invitro_fit.py`
- `Code/in-vitro/pkpd_live_dead_model/src/invitro_fitting.py`
- `Code/in-vitro/pkpd_live_dead_model/tests/`
- `Data/in-vitro/pkpd_live_dead_model/raw/Gemcitabine_PlateMap_20240111.xlsx`
- `Data/in-vitro/pkpd_live_dead_model/raw/drugKinetics/GemcitabineExposure_PKPD.xlsx`
- `Data/GemDelayKillTerm/processed/counts_by_well_time.parquet`
- `Data/GemDelayKillTerm/processed/counts_by_well_time_wellAggregated.parquet`
- `Data/GemDelayKillTerm/plotTimSeries.R`, if retiring the old miningcloneid `plotTimSeries.R`
- `Data/M00_GemcitabinePKPD_101823/nM1000_*.txt`, if deleting the old raw PK text files
- `Data/in-vitro/pkpd_live_dead_model/invitro_fitting_outputs/alsoGoodFit_20260514T093906/`
- `Data/in-vitro/pkpd_live_dead_model/invitro_fitting_outputs/bestFitSoFar_20260513T164159/`
- `Data/in-vitro/pkpd_live_dead_model/input_manifest.tsv`

Record the Gemcitabine-model commit hash and the scoped dirty status in the miningcloneid removal commit message.

## Source-To-Target Migration Manifest

Before deleting a source path, create a manifest for the exact source ref and target ref. The manifest should have one row per source asset:

```text
source_ref	source_path	source_blob	target_ref	target_path	target_blob	disposition	validation_evidence
```

Allowed `disposition` values:

- `identical`: source and target blobs match.
- `migrated-with-diff`: target is the intentional maintained version and validation has passed.
- `retired`: no target copy is needed; record the reason.

Required rows for branches where these files exist:

- `code/invitro_fitting.py` -> `Code/in-vitro/pkpd_live_dead_model/src/invitro_fitting.py` (`migrated-with-diff`; require target tests/validation).
- `code/plot_invitro_fit_outputs.py` -> `Code/in-vitro/pkpd_live_dead_model/plot_invitro_fit_outputs.py` (`migrated-with-diff`; require saved-summary reproduction validation).
- `code/invitro_fitting_outputs/...` -> `Data/in-vitro/pkpd_live_dead_model/invitro_fitting_outputs/...` (`identical` or `migrated-with-diff`, depending on regenerated outputs).
- `data/InVitroData_Gemcitabine/Gemcitabine_PlateMap_20240111.xlsx` -> `Data/in-vitro/pkpd_live_dead_model/raw/Gemcitabine_PlateMap_20240111.xlsx`.
- `data/InVitroData_Gemcitabine/drugKinetics/GemcitabineExposure_PKPD.xlsx` -> `Data/in-vitro/pkpd_live_dead_model/raw/drugKinetics/GemcitabineExposure_PKPD.xlsx`.
- `data/InVitroData_Gemcitabine/processed/counts_by_well_time*.parquet` -> `Data/GemDelayKillTerm/processed/counts_by_well_time*.parquet`.
- `data/InVitroData_Gemcitabine/plotTimSeries.R` -> `Data/GemDelayKillTerm/plotTimSeries.R` (`migrated-with-diff`, if retiring the old copy).
- `data/InVitroData_Gemcitabine/drugKinetics/nM1000_*.txt` -> `Data/M00_GemcitabinePKPD_101823/nM1000_*.txt`.
- `data/InVitroData_Gemcitabine/drugKinetics/README.md` -> either a confirmed target documentation path or `retired` with an explicit reason. No target equivalent was confirmed during this review.

## Current Branch Inventory

Focused scans were run from `/Users/4470246/Downloads/miningcloneid` for:

```bash
git ls-tree -r --name-only <branch> | rg "invitro_fitting|InVitroData_Gemcitabine|plot_invitro|invitro_fitting_outputs|GemcitabineExposure|Gemcitabine_PlateMap|counts_by_well|drugKinetics|test_invitro_fitting_extracellular_decay|joint_fit_summary|dfdctp|alsoGoodFit|bestFitSoFar"
git grep -n -e "invitro_fitting" -e "InVitroData_Gemcitabine" -e "GemcitabineExposure_PKPD" -e "counts_by_well_time" -e "Gemcitabine_PlateMap_20240111" <branch> -- .
```

Initial local tracked candidates observed before exact-ref freezing:

| Branch | Tracked migrated PK/PD assets found |
| --- | --- |
| `pujan_branch_edits` | `code/invitro_fitting.py`, `code/__pycache__/invitro_fitting.cpython-39.pyc`, `data/InVitroData_Gemcitabine/`, `oxygen/tests/test_invitro_fitting_extracellular_decay.py` |
| `main` | `data/InVitroData_Gemcitabine/` |
| `doxoFits` | None found by the focused scan |
| `soft_coupling` | None found by the focused scan |
| `feat-O2-G-resource` | None found by the focused scan |

This table is only an initial local snapshot. It is not sufficient as an execution manifest because local and remote branch tips may differ. In particular, `origin/pujan_branch_edits` also tracks `code/plot_invitro_fit_outputs.py` and `code/invitro_fitting_outputs/...`, which were not present in the local `pujan_branch_edits` snapshot used for the first inventory.

The current `soft_coupling` working tree also has untracked generated artifacts such as `code/invitro_fitting_outputs/`, `code/__pycache__/`, and gemcitabine response PNGs. Those are local working-tree artifacts, not tracked branch content. Do not include them in a tracked cleanup commit unless they are intentionally being ignored or documented.

Important nuance: `data/InVitroData_Gemcitabine/` is not only passive data. On `pujan_branch_edits` and `main`, it also contains `data/InVitroData_Gemcitabine/plotTimSeries.R`, which reads `processed/counts_by_well_time_wellAggregated.parquet` and writes plots from inside that folder. A grep excluding `code/invitro_fitting.py` and excluding files inside `data/InVitroData_Gemcitabine/` found no external tracked consumers of the folder, but the folder's own plotting script is a separate code asset. The folder also includes raw PK text files and `drugKinetics/README.md`. Do not delete the whole folder unless every file has a source-to-target manifest row and `plotTimSeries.R` plus `drugKinetics/README.md` are migrated or explicitly retired.

## Removal Boundary

Remove these if present and tracked on a target branch:

- `code/invitro_fitting.py`
- `code/plot_invitro_fit_outputs.py`
- `code/invitro_fitting_outputs/`
- `code/__pycache__/invitro_fitting*.pyc`
- `data/InVitroData_Gemcitabine/`, only after completing the source-to-target manifest and deciding what to do with `plotTimSeries.R`, `drugKinetics/README.md`, raw PK text files, and processed parquet files
- `oxygen/tests/test_invitro_fitting_extracellular_decay.py`, because it imports the removed `code/invitro_fitting.py` directly
- docs that only describe running the migrated PK/PD live/dead model from miningcloneid, if any are found by grep

Do not remove these without a separate dependency audit:

- `data/InVivoData_Gemcitabine/`
- `data/DrugResponseData_PloidyJumps/` gemcitabine or doxorubicin response files
- `oxygen/code/in-vitro-utils/`
- `oxygen/code/O2_supply_demand_MAP/`
- `oxygen/code/O2G_supply_demand_MAP/`
- `oxygen/data/O2_supply_demand/parameter_table_invitro*.csv`
- `oxygen/code/*/util/*fit_invitro_backend.R`
- `oxygen/code/*/vis/viz_invitro_model_*_results.R`
- generic Python/R dependencies such as `numpy`, `pandas`, `scipy`, `matplotlib`, `pyarrow`, or `openpyxl`, because other code may use them

## Implementation Steps

1. Work from a clean tracked state.

   Because the active checkout currently has many unrelated untracked files, either use clean temporary worktrees per branch or commit from a checkout where `git status --short` has no tracked modifications.

2. Freeze exact refs.

   Fetch and record the exact source refs before building manifests:

   ```bash
   git fetch origin
   git rev-parse pujan_branch_edits origin/pujan_branch_edits doxoFits origin/doxoFits soft_coupling origin/soft_coupling feat-O2-G-resource origin/feat-O2-G-resource main origin/main
   git rev-list --left-right --count <local-branch>...origin/<branch>
   ```

   Decide whether the cleanup targets local branch tips or `origin/*` tips. Build all deletion manifests from those exact SHAs, not from ambiguous branch names.

3. Verify the target migration.

   In `Gemcitabine-model`, verify the canonical code/data paths above from the committed tree. Also run or inspect the latest migration validation status. Do not proceed with source deletion if the migrated copy is incomplete or if a source file lacks an `identical`, `migrated-with-diff`, or `retired` manifest row.

4. Build a branch-specific removal manifest.

   For each target branch/ref, run the inventory commands above and save the output in the terminal log or commit notes. Only delete tracked files that are present on that exact ref.

5. Remove migrated files.

   Use branch-specific `git rm --ignore-unmatch` commands for the removable paths. Expected branch-specific effects:

   - `pujan_branch_edits` or `origin/pujan_branch_edits`: remove old Python module, stale pyc, plotter/output artifacts if present on the chosen ref, and the Python test tied to the removed module. Remove `data/InVitroData_Gemcitabine/` only if its full per-file manifest is complete.
   - `main` or `origin/main`: remove `data/InVitroData_Gemcitabine/` only if its full per-file manifest is complete.
   - `doxoFits`, `soft_coupling`, `feat-O2-G-resource`: likely no tracked removals unless a fresh inventory finds new matching files.

6. Re-scan references after deletion.

   Required clean grep:

   ```bash
   git grep -n -e "invitro_fitting" -e "InVitroData_Gemcitabine" -e "GemcitabineExposure_PKPD" -e "counts_by_well_time" -e "Gemcitabine_PlateMap_20240111" -- .
   ```

   Acceptable remaining matches are only migration/removal documentation. Any code or test match means either the file should also be removed or the reference should be rewritten to point to `Gemcitabine-model`.

7. Keep shared dependencies.

   Do not edit requirements, environment files, or shared data directories unless a repo-wide search proves they are only used by the removed migrated code. If a dependency looks unused, leave it for a later dedicated dependency cleanup.

8. Run lightweight validation per branch.

   Minimum checks:

   ```bash
   git grep -n -e "invitro_fitting" -e "InVitroData_Gemcitabine" -- .
   python3 -m compileall code
   ```

   On branches with oxygen tests, also run the branch's normal oxygen smoke/unit test entrypoint if the required R/Python dependencies are available:

   ```bash
   Rscript oxygen/tests/run_unit_tests.R
   ```

   If tests are too expensive or dependencies are unavailable, record the exact reason in the final commit notes.

9. Commit and propagate.

   Preferred sequence:

   - Make branch-specific cleanup commits for refs with different manifests.
   - Cherry-pick only when two branches have the same deletion manifest and the same validation outcome.
   - Treat `doxoFits`, `soft_coupling`, and `feat-O2-G-resource` as no-op branches if refreshed exact-ref inventories still find no tracked migrated assets.
   - If a branch has no tracked migrated assets, leave it unchanged and record that the branch was already clean.

## Acceptance Criteria

- `Gemcitabine-model` has the canonical migrated code and data verified from a committed tree, with target commit hash recorded.
- Every deleted source file has a source-to-target manifest row with `identical`, `migrated-with-diff`, or `retired` disposition.
- No target miningcloneid branch tracks the migrated PK/PD live/dead fitting or visualization code.
- No target miningcloneid branch tracks `data/InVitroData_Gemcitabine/` unless the branch still intentionally owns `plotTimSeries.R`, `drugKinetics/README.md`, or another non-migrated user of that data.
- No code/test references remain to removed paths.
- Shared oxygen, in-vivo gemcitabine, drug-response, and generic dependency files remain untouched.
- Branches with no tracked migrated files are explicitly reported as already clean.
