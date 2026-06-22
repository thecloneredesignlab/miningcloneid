# Optimization Data Streams Figure Implementation Plan

## Purpose

Build an automated figure that shows the measured data streams used by the in vivo and in vitro O2 supply-demand optimizations. The figure should make three points clear at a glance:

- the same isogenic SUM159 2N and 4N starting cell-line pair underlies both the in vitro and in vivo experiments;
- time and growth history are the organizing structure in both settings;
- oxygen exposure is known and experimentally controlled in vitro, but is latent and unmeasured in vivo.

The figure must be generated from the actual source data, not from static report text or hand-entered counts. If the input data change, the normalized tables and figure contents should change with them.

## Inputs Reviewed

This plan is based on the two data-stream reports:

- `docs/in_vivo_optimization_data_streams_report.md`
- `docs/in_vitro_optimization_data_streams_report.md`

It also incorporates six parallel review passes:

- two data-extraction reviews;
- two visual-design reviews;
- one implementation-architecture review;
- one validation/QC review.

The main consensus from those reviews is:

- implement the automated figure in R because the existing loaders, RDS objects, config semantics, and plotting stack are R-based;
- use one shared absolute-day x-axis for the main overview, with no axis break in the primary panel;
- show in vitro known oxygen as colored passage blocks;
- show in vivo oxygen as "not measured / latent", not as a colored oxygen timeline;
- emit intermediate TSVs and a manifest so the figure is auditable.

## Proposed Script Location

Add a new runnable script:

- `oxygen/code/O2_supply_demand_MAP/analysis/plot_optimization_data_streams.R`

Add shared helpers only if the script becomes too large:

- `oxygen/code/O2_supply_demand_MAP/util/o2_supply_demand_map_data_streams.R`

Default output directory:

- `figures/`

Expected primary outputs:

- `figures/optimization_data_streams_overview.pdf`
- `figures/optimization_data_streams_overview.png`
- `figures/in_vivo_optimization_data_streams.pdf`
- `figures/in_vivo_optimization_data_streams.png`
- `figures/in_vitro_optimization_data_streams.pdf`
- `figures/in_vitro_optimization_data_streams.png`

## Implementation Language and Dependencies

Use R.

Required packages should stay close to the existing workflow:

- `ggplot2`
- `dplyr`
- `tidyr`
- `readxl`
- `yaml`
- `patchwork`
- `scales`

Reuse existing loaders and semantics rather than duplicating them:

- `oxygen/code/O2_supply_demand_MAP/util/o2_supply_demand_map_shared.R`
- `oxygen/code/in-vitro-utils/io.R`
- `oxygen/code/in-vitro-utils/lineage_adapter.R`
- `oxygen/code/in-vitro-utils/summaries.R`, if useful for observed-data summaries

## Command Interface

The script should accept:

- `--config=oxygen/config/O2_supply_demand.yaml`
- `--out_dir=figures`
- `--data_dir=...`, optional override for in vivo data
- `--fit_objects_dir=...`, optional override for in vitro RDS objects
- `--flow_density_path=...`, optional override for flow density
- `--fit_dir=...`, optional completed-fit directory; if present, prefer its `fit_config.rds` for active config flags
- `--include_source_only=TRUE`, default `TRUE`
- `--strict_snapshot=FALSE`, default `FALSE`

`--strict_snapshot` should be reserved for regression testing against the current repository snapshot. The normal plotting path should recompute counts and update automatically when source data change.

## Figure Concept

### Main Overview

Use one shared absolute x-axis in days. The current data suggest a natural range of day 0 to about day 160, but this should be computed dynamically from the maximum of:

- in vivo final observed harvest days;
- in vitro terminal lineage cumulative days.

The main layout should use four stacked lanes:

1. in vitro 2N;
2. in vitro 4N;
3. in vivo 2N;
4. in vivo 4N.

This makes the matched isogenic 2N and 4N pair visible in both environments, and makes the absolute time-frame comparison immediate.

### In Vitro Lanes

In vitro should be organized as cumulative passage time:

- each passage segment is a rectangle from `cumulative_start_day` to `cumulative_end_day`;
- rectangle fill encodes known target O2 percentage;
- oxygen severity should use a stable monotone scale, with 20.5% as least severe and 0% as most severe;
- passage or depth labels can be shown selectively to avoid clutter;
- control and oxygen-deprived lineages should be distinguishable within each 2N/4N lane.

Attach measured data streams to the passage timeline:

- growth: interval glyphs spanning the passage duration, because growth is measured per passage interval;
- karyotype: endpoint glyphs at the relevant passage end, with cell count encoded by size or a small distribution summary;
- flow: endpoint glyphs or compact density mini-panels at the relevant passage end.

Do not plot in vitro growth as a continuous daily cell-count trajectory in this source-data figure. The measured stream is passage-level growth, not a daily time course.

### In Vivo Lanes

In vivo should be organized as measured tumor burden trajectories:

- draw observed tumor-burden points and connecting lines from day 0 to the last observed day;
- facet or lane by starting cohort, 2N versus 4N;
- use subdued line styling for source-only observations and stronger styling for observations included under the active fit config;
- mark each harvest day explicitly.

Attach terminal data streams only at harvest:

- terminal single-cell ploidy/chromosome-count data from `all_ploidy.csv`;
- terminal histology necrosis fraction from the mapping CSV.

Show in vivo oxygen as unmeasured or latent:

- use a gray hatched or outlined "O2 not directly measured" track;
- do not color in vivo time intervals by oxygen unless adding a separate, clearly labeled model-prediction panel later.

### Optional Detail Panels

The main overview should prioritize time placement. If the overview becomes crowded, add detail panels below or in a companion page:

- in vitro flow-density distributions for the measured passages;
- in vitro karyotype distributions by passage;
- in vivo terminal ploidy distributions by harvest;
- in vivo necrosis percentages by harvest.

These detail panels should reuse the same event IDs and day positions from the normalized data tables.

## Data Extraction Work Packages

### Work Package 1: Config and Path Resolution

Build a `config_snapshot` table with:

- resolved `config_path`;
- resolved `data_dir`;
- resolved `fit_objects_dir`;
- resolved `flow_density_path`;
- `dose_zero_only`;
- `paired_only`;
- `burden_exclude_day0`;
- `truncate_at_treatment`;
- `ploidy_at_harvest`;
- `use_necrosis_loss`;
- `necrosis_mapping_csv`;
- `start_with`;
- timestamp and git SHA if available.

If `--fit_dir` is supplied and contains `fit_config.rds`, use that config for active inclusion flags. Otherwise, read the YAML config.

### Work Package 2: In Vivo Normalized Tables

Read all in vivo sources dynamically:

- burden workbook from `dt_Gem_VT_20260209_v5.xlsx`;
- terminal ploidy table from `all_ploidy.csv` or `.tsv`, preserving the existing tab-delimited `.csv` fallback;
- necrosis mapping CSV from `necrosis_mapping_csv`.

Build:

- `invivo_harvest_catalog.tsv`: one row per burden harvest with cohort, dose, treatment day, first observed day, harvest day, stream-presence flags, and active-fit inclusion flags.
- `invivo_burden_long.tsv`: one row per finite `Day_*` value with numeric day, burden, source presence, and `used_in_burden_loss`.
- `invivo_ploidy_cells.tsv`: one row per terminal cell with harvest, cell ID, raw ploidy, raw total chromosomes, selected endpoint value, endpoint mode, harvest day, and `used_in_endpoint_loss`.
- `invivo_necrosis_endpoint.tsv`: one row per mapped harvest with necrosis fraction, source count, harvest day, and `used_in_necrosis_loss`.

Implementation notes:

- derive cohort from `harvest` containing `2N` or `4N`;
- derive ploidy harvest by stripping `.sps.cbs` from `file`;
- average finite mapped necrosis values by `dt_harvest`;
- keep day 0 in the display table even when `burden_exclude_day0` excludes it from the burden likelihood;
- distinguish source availability from active objective inclusion.

### Work Package 3: In Vitro Normalized Tables

Load in vitro sources through existing fit-object helpers:

- `fit_data.Rds`;
- `jobs_2N.Rds`;
- `jobs_4N.Rds`;
- `g0g1_ploidy_density_grid.csv`, attached by the canonical flow loader.

Build:

- `invitro_segment_catalog.tsv`: one row per job segment with cohort, segment ID, parent segment ID, depth, oxygen percent, and data IDs.
- `invitro_lineage_timeline.tsv`: one row per segment per terminal path, including terminal key, lineage label, passage index, cumulative start day, cumulative end day, duration source, and oxygen percent.
- `invitro_passage_observations.tsv`: one row per passage ID with growth, passage duration, initial cells, final cells, karyotype count, flow-grid count, and linked segment metadata.
- `invitro_kary_cells.tsv`: one row per observed karyotype cell, attached to segment and cumulative passage time.
- `invitro_flow_density_grid.tsv`: one row per flow grid point, with canonical sample name, passage ID, ploidy, density, log density, and cumulative passage time.

Implementation notes:

- compute segment duration as the median finite observed passage duration across the segment's data IDs;
- fall back to 14 days only when no finite duration is available, matching the in vitro adapter behavior;
- derive terminal paths from the parent-child job graph rather than hard-coded path names;
- treat oxygen values as percent O2, not fractions;
- duplicate shared segments per terminal lineage path when building the timeline, because the figure needs terminal lineage histories.

### Work Package 4: Unified Measurement Events Table

Build a common `measurement_events.tsv` table to drive legends and event layers:

- `context`: `in_vivo` or `in_vitro`;
- `cohort`: `2N` or `4N`;
- `stream`: burden, ploidy, necrosis, growth, karyotype, flow, oxygen_design, oxygen_latent;
- `entity_id`: harvest, segment, passage, or sample identifier;
- `time_start_day`;
- `time_end_day`;
- `measurement_level`: interval, endpoint, distribution, trajectory, design;
- `used_in_objective`;
- `source_file`;
- `n_observations`;
- `display_label`.

This table should be generated from normalized source tables, not hand-authored.

## Plotting Work Packages

### Work Package 5: In Vitro Panel

Use `geom_rect` for oxygen passage blocks:

- `xmin = cumulative_start_day`;
- `xmax = cumulative_end_day`;
- `ymin/ymax` define lane position;
- fill = oxygen percent severity.

Overlay source-data glyphs:

- growth intervals as bracket/line glyphs across segment duration;
- karyotype observations as endpoint markers sized by cell count;
- flow observations as endpoint markers shaped differently from karyotype or as compact density sparklines if feasible.

Faceting or lane ordering:

- top-level cohort rows: 2N and 4N;
- within each cohort, separate control and oxygen-deprived terminal paths;
- preserve consistent ordering between 2N and 4N.

### Work Package 6: In Vivo Panel

Use measured burden trajectories as the central visual element:

- `geom_line` and `geom_point` by harvest;
- x = numeric day;
- y = tumor burden, likely log-scaled;
- facet or lane by cohort.

Overlay terminal endpoint glyphs:

- ploidy/sequencing glyph at harvest day;
- histology necrosis glyph at harvest day;
- optional size or color intensity encodes number of cells or percent necrosis.

Add an in vivo oxygen-status track:

- gray, hatched, or outlined;
- labeled "O2 not directly measured; latent in model";
- no oxygen severity colors.

### Work Package 7: Combined Overview

Use `patchwork` to stack:

- title/header strip explaining the common isogenic 2N/4N pair;
- in vitro source-data panel;
- in vivo source-data panel;
- shared legend area.

Use a common x-axis limit across panels:

- minimum day 0;
- maximum day computed from all in vivo and in vitro timelines, rounded up to a clean tick.

The main overview should not use an axis break. A secondary detail figure may use axis breaks if needed, but the overview must preserve the absolute timeframe comparison.

## Visual Encoding Requirements

Use separate encodings for separate concepts:

- cohort: lane/facet labels, not oxygen color;
- context: panel labels, in vitro versus in vivo;
- oxygen: in vitro block fill only;
- measured stream type: glyph shape or small icons;
- active objective inclusion: opacity, outline, or line weight;
- source-only observations: lighter styling, still visible.

Avoid these failure modes:

- do not color in vivo time by oxygen severity;
- do not plot in vivo ploidy or necrosis as time courses;
- do not plot in vitro growth as daily continuous cell counts;
- do not let the oxygen color scale encode cohort identity;
- do not hide observations that exist in source data but are excluded by the active config.

## Validation and QC

### Default Runtime Checks

Treat these as hard errors:

- missing required files;
- missing required columns;
- nonnumeric day columns after parsing `Day_*`;
- no finite burden observations;
- ploidy table cannot be read with the existing CSV/tab fallback;
- necrosis values outside `[0, 100]`;
- job `data_ids` missing from `fit_data`;
- duplicate or cyclic in vitro parent-child job graph;
- ambiguous or failed flow attachment for samples that should match;
- oxygen values not finite or outside plausible percent range.

Treat count changes as warnings by default, not errors, so the figure can update when data change.

### Optional Strict Snapshot Checks

Under `--strict_snapshot=TRUE`, compare against the current source snapshot and fail on drift unless the expected snapshot is intentionally refreshed.

Current in vivo snapshot diagnostics:

- 40 burden rows;
- 31 `Day_*` columns;
- 14,125 terminal ploidy cells across 16 harvests;
- 7 mapped finite necrosis harvests;
- 10 untreated burden rows;
- 8 untreated rows with terminal ploidy under the paired filter;
- 6 paired untreated rows with necrosis.

Current in vitro snapshot diagnostics:

- 131 fit-data entries;
- 75 job segments;
- 6 terminal paths;
- terminal cumulative durations from 57 to 151 days;
- 114 growth entries;
- 12 karyotype passages with 220 observed cells;
- 20 flow samples with 4,000 total grid rows;
- oxygen values: 0, 0.1, 0.2, 0.3, 0.5, 1, 2, 20.5.

### Figure Smoke Tests

After rendering, verify:

- output PDF and PNG exist and are nontrivial size;
- both 2N and 4N appear in both contexts;
- in vitro timeline reaches the maximum terminal cumulative day;
- every in vitro passage segment has an oxygen block;
- in vivo burden trajectories cover day 0 through harvest;
- in vivo endpoint glyphs occur only at harvest days;
- in vitro growth, karyotype, and flow layers have nonzero rows when source data contain them;
- source-only and active-fit observations are visually distinguishable.

## Subtask Breakdown

### Subtask 1: Data Loader and Normalized Tables

Implement path/config resolution and source-driven extraction.

Deliverables:

- `config_snapshot.tsv`;
- all in vivo normalized TSVs;
- all in vitro normalized TSVs;
- `measurement_events.tsv`.

Acceptance criteria:

- rerunning after source changes updates tables automatically;
- no figure layer reads from static report prose;
- all joins expose source-only and active-in-objective flags.

### Subtask 2: In Vitro Timeline Panel

Implement the passage-lineage timeline with known oxygen exposure.

Deliverables:

- `in_vitro_optimization_data_streams.pdf`;
- `in_vitro_optimization_data_streams.png`.

Acceptance criteria:

- oxygen blocks cover every passage;
- cumulative day scale is proportional to passage duration;
- growth, karyotype, and flow observations are placed on the relevant passage intervals/endpoints;
- 2N and 4N rows are visually paired.

### Subtask 3: In Vivo Timeline Panel

Implement measured tumor-burden trajectories with terminal endpoint data.

Deliverables:

- `in_vivo_optimization_data_streams.pdf`;
- `in_vivo_optimization_data_streams.png`.

Acceptance criteria:

- burden trajectories are plotted over observed days;
- ploidy and necrosis appear only at harvest endpoints;
- active objective subset is distinguishable from source-only data;
- oxygen is labeled as latent/unmeasured rather than color-coded.

### Subtask 4: Combined Figure Assembly

Assemble the shared-axis overview using the two panels.

Deliverables:

- `optimization_data_streams_overview.pdf`;
- `optimization_data_streams_overview.png`.

Acceptance criteria:

- a single absolute-day x-axis governs both settings;
- the same isogenic 2N/4N pair is clear in labels or header;
- in vitro and in vivo timeframes can be compared without reading the caption;
- legends separate oxygen, cohort, stream type, and inclusion status.

### Subtask 5: Validation and Documentation

Add runtime manifest, smoke checks, and a short usage note.

Deliverables:

- `optimization_data_streams_manifest.tsv`;
- optional `optimization_data_streams_qc.tsv`;
- usage example in either the script header or workflow README.

Acceptance criteria:

- schema/join failures are caught early;
- count drift is reported;
- strict snapshot mode can be used for regression checks.

## Suggested First Implementation Pass

Start with a robust overview rather than highly detailed mini-distribution glyphs.

First pass:

- in vitro oxygen blocks plus event glyphs for growth, karyotype, and flow;
- in vivo burden curves plus endpoint glyphs for ploidy and necrosis;
- combined shared-day layout;
- normalized TSV outputs and manifest.

Second pass:

- add compact density sparklines or small multiples for flow;
- add richer terminal ploidy/karyotype distribution summaries;
- tune legends and manuscript-ready styling;
- add optional detail figure pages.

This staging keeps the first implementation focused on the central biological message: the same isogenic 2N/4N cell-line pair is observed across two long experimental contexts, but in vitro has a known progressive oxygen history while in vivo oxygen is latent and inferred from growth plus terminal biological states.
