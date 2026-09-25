# Stage 1 regression: in-vitro-utils migration

Stage 1 moved the complete former `oxygen/code/in-vitro-utils/` implementation
into `O2_supply_demand_MAP/util/` and `O2_supply_demand_MAP/vis/invitro/`.

## Public API

- All 58 historical `ivt_*` functions remain available through
  `util/o2_supply_demand_map_invitro_utils.R`.
- Function names and formal arguments were compared before and after migration.
- The optimizer specification still contains the same 20 fitted parameters.
- Loader tests were run from the repository root and from an unrelated working
  directory.

## Seed10 golden comparison

- Objective: `3.8525352626059366`.
- Total, growth, ploidy, and flow objective components were frozen.
- All ten scientific output tables have identical row counts and byte content.
- The 22 API/table/objective comparisons pass at tolerance `1e-12`.
- Regression output: `/tmp/oxygen_stage1_invitro`.

## Visualization

The migrated plot constructors reproduced the Stage 0 isolated-output fixture:

`/tmp/o2_invitro_viz_seed10_baseline.AoMnnw`

All eight PDF and eight PNG names are present, and all eight fixed-seed PNG
files are byte-identical.

## Stage gate

- R parse: 126/126 files.
- Shell syntax: 24/24 scripts.
- In-vitro migration/full test gate: 66/66 tests and 352 expectations.
- References to the removed `oxygen/code/in-vitro-utils/` runtime path: zero.
- Protected model/optimizer tree: 18/18 files unchanged.
- Gate logs: `/tmp/oxygen_phase1_full_gate/_gate_logs/`.
