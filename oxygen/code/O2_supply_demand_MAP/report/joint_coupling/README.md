# Joint coupling report

This directory contains the report-only layer for the joint soft-coupling and
ploidy-category workflow. It consumes materialized analysis tables and existing
visualization files; it does not recalculate statistics, simulate trajectories,
or draw figures.

## Files

| File | Function |
|---|---|
| `render_joint_coupling_report.R` | Validates the report inputs, copies portable PNG/PDF assets, and renders `joint_coupling_analysis_report.html` with a responsive left navigation bar. |
| `joint_coupling_figure_catalog.tsv` | Defines the authoritative order and report text for all 29 results: major section, subsection ID, figure title, legend explanation, interpretation, and limitation. |
| `README.md` | Documents the report contract and generated files. |

## Report structure

Each of the 29 results is rendered as its own numbered subsection and contains
exactly one numbered figure, a descriptive title, a plot-reading legend, a
result interpretation, and a result-specific limitation. The seven result
groups are overview, within-pair stability, between-pair stability, biological
processes, ploidy categories, Cat-by-Class association, and robustness. The
left navigation contains both group links and one link for every result.
The navigation is viewport-aware: the result group containing the active
subsection expands automatically, other result groups collapse, the active
figure link is highlighted, and the sidebar scroll position follows that link.
Each result group also has an accessible manual toggle with synchronized
`aria-expanded` state.

Before the numbered results, Section 2.1 defines both classification systems.
It states the in-vivo/in-vitro ratio direction and inclusive ClassB boundaries,
explains the biological reading of ClassA/B/C, and distinguishes these ratio
classes from the CatA/B/C/U in-vivo ploidy-trajectory categories. The Cat rules
include the operational high-ploidy, terminal-low, drop, plateau, and
two-transition BIC criteria, plus the seed-level 2N/4N combination rules.

The HTML shell follows the project's fit-report visual system in
`report/render_fit_report.R` and `report/fit_result_report.Rmd`: a light-gray
page background, white rounded report cards, a dark-blue gradient sidebar
header, hierarchical blue navigation, restrained blue/orange callouts, and the
same 1100/900-pixel responsive breakpoints. The FixO2 report supplies the
secondary callout and figure-caption conventions. Shared argument, escaping,
and numeric-display helpers are sourced from `util/`.

The renderer rejects missing figures, undocumented figures, duplicate figure
stems, missing PDF companions, incomplete catalog text, and any figure count
other than 29. This prevents the HTML from silently omitting a result when the
visualization suite changes.

## Generated output contract

The report output directory contains:

| Generated file | Function |
|---|---|
| `joint_coupling_analysis_report.html` | Self-contained report markup and styling with relative local asset links. |
| `figures/*.png` and `figures/*.pdf` | Report-local raster previews and vector downloads for every numbered figure. |
| `report_figure_catalog.tsv` | Validated figure order and all displayed report metadata, including report-local asset paths. |
| `chart_map.tsv` | Visualization provenance enriched with report figure number, section, title, legend, interpretation, and limitation. |
| `report_manifest.tsv` | Inventory of the HTML, metadata files, and all 58 figure assets. |

Sample sizes, the ClassB threshold, in-vitro anchor, and Cat pair counts are
read from materialized config and quality tables rather than hard-coded. The
Cat definition text additionally consumes
`ploidy_coupling/ploidy_category_definition.tsv`, so its displayed numeric
thresholds remain synchronized with the materialized analysis contract. The
HTML remains portable when the entire report directory is copied elsewhere.
