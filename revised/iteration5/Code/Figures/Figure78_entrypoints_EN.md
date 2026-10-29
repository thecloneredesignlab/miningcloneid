# Figure7/8 unified entry points

The same current code and private model layout is used locally and on HPC.
`Code/config/manifests/figure78_entrypoints.tsv` is the complete entry inventory.

Figure7 update (2026-10-24): generation now dispatches the native-fit passage /
last-five cycling workflow in `Figure7_current_pipeline_EN.md`. Its320 conditions
expand to660 grid workers; `submit_Figure7_current.sh` submits the full array and
dependent numerical/grid collection. Figure7 generation and HPC drawing default
to current products only; historical products require explicit selection and are
not regenerated here. Figure8 is unchanged. During integration, `static
--source-only=yes` explicitly does not certify old numerical data; final complete
data/provenance checks remain mandatory.

## Paths and source policy

- Model: workspace-relative `Code/vendor/figure78_model`; all 12 SHA-256 values
  must match `figure78_model_sources.tsv`. A shared/external model override fails.
- Original fits: workspace-relative sibling `soft_couping_org/oxygen/results/`
  and the approved final joint-fit campaign in `figure78_inputs.json`. An explicit
  `--fit-root` / `FIGURE_JOINT_RESULT_ROOT` may relocate these original fits.
- Existing results: `data/Figures/Figure7` and `Figure8`, including their immutable
  provenance; empirical source data: `data/share_bundle_figures_3_5`.
- Old version directories and archives are not runtime roots or fallbacks.
  Provenance snapshots inside current data remain identity evidence, not numerical
  implementations. Their original receipts and RDS fingerprints are not rewritten.

R entries resolve paths themselves. Shell/HPC entries automatically load
`Code/config/figure78_hpc_paths.env`; no login profile is modified. The shared HPC
model is never patched. `figure78_Rscript.sh` uses Rscript from PATH locally, or an
explicit `FIGURE78_NATIVE_RSCRIPT`; HPC uses the approved R 4.4.2 Apptainer image.
`FIGURE78_CONTAINER` can identify another separately approved environment.

## Entry mapping

| Purpose | Current entry | Behavior |
|---|---|---|
| Generate Figure7 | `run_Figure7.R`, `data_Figure7.R` | Current preparation, tasks, aggregation, validation, plotting |
| Generate Figure8 | `run_Figure8.R`, `data_Figure8.R` | Same current forward pipeline |
| Draw retained results | `draw_Figure7.R`, `draw_Figure8.R`, individual wrappers | No simulation or archive access; 26 products |
| Validate retained/full result receipts | `validate_Figure78.R` | Metadata and identity only; may write audit receipts |
| Static audit | `figure78_entry_tools.py static` | Python AST, R parse-only, shell syntax, paths and hashes |
| Host/HPC dispatch | `Code/hpc/run_Figure78.sh` | Explicit figure and stage; never submits |
| One prepared task | `Code/hpc/run_Figure78_task.sh` | Explicit data root, full manifest and Slurm index required |
| Future approved submission | `Code/hpc/submit_Figure78.sh` | Explicit `--submit=yes`; xxlarge; full manifest, no concurrency cap |

Generation and drawing remain separate. Retained current results have no new
`prepared/task_manifest.tsv`; they are not resumable prepared simulations.
Generating data requires a separately approved empty output root or the existing
explicit rebuild protection, never an implicit overwrite of retained tables.

Old data/finalize/audit/HPC filenames are compatibility guards, not alternate
workers. Explicit current `--stage=...` calls may use generation aliases. Legacy
`--action`, task IDs, configs, positional arguments and version directories fail
with migration instructions instead of being silently mapped to different tasks.
Old validation aliases require `--mode=retained` (R/Python) or `--mode=static`
(Python). Old Slurm spool adapters fail; they cannot replay a retired array index.
Retired implementations/configs are backed up outside the active Code tree and are
not loaded by any current entry.

## Safe static command

From iteration5, with Python 3.6+ and an existing approved R environment:

```sh
python3 Code/Figures/figure78_entry_tools.py static --workspace=. --report=audit/figure78_static/checks.json
```

Only the audit report is written. R files are parsed, never sourced or evaluated.
No `plan`, `prepare`, simulation, plot, numerical acceptance or Slurm command is
executed. Static success does not establish numerical or graphics runtime success.

## Commands requiring separate authorization

```sh
bash Code/hpc/run_Figure78.sh --figure=8 --stage=retained-validate
bash Code/hpc/run_Figure78.sh --figure=7 --stage=plot
```

The first authenticates result receipts and may write audit records; the second
redraws all Figure7 analyses. Neither is part of the static deployment acceptance.
The cold-start, independent numerical acceptance and render-regression tools also
execute scientific/render operations and require their own explicit authorization.
Existing calculation definitions, scientific configs, kernels, random streams,
plot builders, styles and publication filenames are unchanged.
