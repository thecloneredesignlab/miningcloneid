# Figure7 current production workflow

`run_Figure7.R` is the complete entry; `data_Figure7.R` calls it. Figure7-only
dispatch is in `util/analysis/figure7_current_pipeline.R`, using the independently
accepted `figure7_passaged.R/.cpp`. Figure8 and historical `_ver2` are not generated
or written by this workflow.

## Scientific contract

Use complete C01 seed94 / C02 seed80 joint winners:26 in-vivo and33 in-vitro
parameters. Culture configuration comes from `fit_result.rds$ctx$invitro_cfg`.
Scenario oxygen, native p_misseg=.1/.3, pure44/88-chromosome roots, designed blockade
and the prediction grid are explicit settings; effective missegregation retains
its native state/oxygen dependence. The private12-source lock and HPC shared
model are not modified.

Replay the original experimental DAG, then cycle each control/O1/O2 route's last
five complete measured tuples in order. Select the positive observation day
nearest the target over the full candidate window, then rescale composition to
the next actual inoculum. No standardized inoculum, threshold-triggered random
sampling, insufficient-stock stop, N<1 stop or unreached-target stop is added.
The1000-day selected-elapsed clock includes the prefix; final selection precedes
truncation. Conditional reseeding can amplify expected numbers and is not proof
of physical inventory feasibility. Detailed equations are in
`manuscript/figure_legends/Figure7_calculation_methods_EN.md`.

A/B are matched daily viable low-daughter shares from high/low mothers, both
divided by HL+LL. They complement when production is positive; day0 and zero
production are NA. Equal-route medians/type-7 IQRs describe three protocol
scenarios, not confidence intervals or200 random repeats. supp_fig7_8 retains
non-inherited first-generation D and permanent ancestry diagnostics at current
low-state proliferation rates.

## Stages and storage

Use `--stage=plan|prepare|simulate|aggregate|validate|plot|all`, `--profile=full`,
an explicit isolated `--data-root=/absolute/campaign/data`, and a separate
`--output-root` for drawing. Prepare/rebuild refuses implicit production overwrite.
Workers accept `--task-index=N` or `--task-ids=N,M` through the same entry.

The320 conditions become660 workers:170 culture conditions at137/139/141N and150
unchanged continuous conditions. Each evaluates three beta settings; culture
uses three routes per beta. Atomic task and route outputs carry SHA receipts
and exact source/input fingerprints. Only validated matching caches are reused.
Unregistered interrupted files must first be preserved in a recoverable audit
location. Each worker has a physically distinct model/compile runtime.

Both adjacent grid refinements must pass: shares≤.5 percentage points, mean
ploidy≤.01N, log size≤.01, full-state CDF≤.01, upper lost-birth fraction≤.001,
and identical passage calendars. A failed certificate blocks publication and
does not authorize automatic grid expansion. Numerical completion and figure
publication are separate states.

## HPC and acceptance

```sh
bash Code/hpc/submit_Figure7_current.sh --data-root=/absolute/prepared/campaign/data --submit=yes
```

Use red partition / xxlarge QoS and the full array index range without a `%`
concurrency cap. Dependent postprocessing only aggregates/validates; it neither
draws nor promotes. Retry explicit failed/missing indices with `--task-ids=...`,
preserving successful matching receipts.

`test_Figure7_current.R` checks formal dispatch, native selection, independent
Euler quantities and selected tasks. `figure78_entry_tools.py static
--source-only=yes` checks syntax, entries, private lock and original inputs only;
it explicitly does not certify retained numerical data against changed code.
Final normal static/render authentication still requires matching completion,
numerical, table and source receipts; source-only cannot replace those gates.

Draw current Figure7 and supplements1–8 only, preserving existing style. Stage
and inspect PDF/PNG outputs before allowlisted promotion. Preserve historical
tables with their original authentication chain, never re-signing them under
the current campaign fingerprint. Figure8 and every `_ver2` data, provenance and
image file remain byte-identical.

`validate_Figure7_current_figures.py` inspects staged current products without
simulating or drawing: nine products and16 matrix components,50 PDF/PNG files
per mirror and86 PDF pages. It exports every page raster for human review.
Numerical results intentionally change, so old-pixel identity is not an
acceptance criterion; preserve dimensions/page counts/style and inspect the
actual data binding, labels and readability.

Small-entry authorization flags in `figure7_passaged.json` govern the standalone
acceptance program only. Formal dispatch reuses its scientific protocol/grid
limits, not that bounded entry to run production. A formal campaign needs an
explicit separate root and this execution's approved scope/backup records.
