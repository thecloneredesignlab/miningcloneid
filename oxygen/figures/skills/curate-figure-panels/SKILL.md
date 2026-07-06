---
name: curate-figure-panels
description: Curate first-pass scientific figure panels from existing reports, rendered plots, and code outputs. Use when Codex must choose among multiple candidate plots based on a figure caption, claim graph, panel backlog, or manuscript figure plan; save one selected candidate per applicable panel; avoid new analysis unless requested; and document selection rationale, provenance, caveats, and skipped alternatives.
---

# Curate Figure Panels

## Overview

Use this skill to turn a figure plan into a concrete first-pass panel set using only existing artifacts. The output should make the available evidence visible, preserve provenance, and record why each panel was chosen over plausible alternatives.

## Workflow

### 1. Establish Panel Contracts

Read the figure caption, claim graph, panel backlog, evidence plan, or manuscript outline that defines what each panel should convey.

For each panel, record a short contract:

- claim or message the reader should take away;
- expected visual form;
- required source type, such as fit report, rendered figure, table, or code-generated plot;
- whether the panel can be satisfied by existing artifacts.

Do not force a panel if no existing code/report/asset supports it. Mark it deferred and state what new work would be required.

### 2. Inventory Existing Evidence

Search for existing assets before creating anything:

- rendered figures: `*.png`, `*.pdf`, `*.svg`;
- report files: `*.html`, especially self-contained HTML reports with embedded images;
- figure/code manifests, evidence plans, and report-generation scripts;
- archives containing generated report outputs.

Prefer `rg --files`, `find`, `unzip -l`, and targeted parsing over broad text dumps from self-contained HTML reports.

For HTML reports, list section IDs and titles first. If images are embedded as base64, extract only the needed section images to a temporary candidate folder before final selection.

### 3. Generate Candidate Set

For each panel contract, identify two or more candidates when available. Include obvious alternatives even if they are later rejected, such as:

- broad overview versus split subpanels;
- compact summary versus detailed diagnostic plot;
- UMAP versus t-SNE versus PCA embeddings;
- parameter-ratio plot versus mechanism-function plot;
- one representative fit versus multi-run summary.

Keep candidate extraction mechanical. Do not modify the analysis or regenerate new results unless explicitly asked.

### 4. Select One Candidate Per Applicable Panel

Choose the candidate that best satisfies the panel contract and is most readable as a standalone panel.

Use these criteria:

- directness: does the plot actually show the claim, not just a related quantity;
- provenance: can the source report/code path be named;
- readability: are axes, legends, and labels visible;
- scope: does the panel answer one panel-sized question rather than an entire figure;
- consistency: does it use the same run, cluster definition, or preprocessing as adjacent panels;
- caveats: does the panel require a provenance warning or later regeneration.

Do not let an existing caption overrule the plot content. Captions and panel choices should describe what the selected artifact actually shows.

### 5. Save Panels

Save only selected panels to the requested output folder. Use stable, descriptive names:

```text
fig<figure><panel>_<short_content_description>.png
```

Examples:

```text
fig3a_invitro_growth_ploidy_burden_fit.png
fig4g_landscape_cluster_driver_violins.png
fig6b_objective_distribution_by_curve_class.png
```

Do not assemble multi-panel figures unless requested.

### 6. Write Selection Log

Create a concise log in the output folder. For each saved panel, include:

- selected source file/report and section ID if applicable;
- candidates considered;
- decision rationale;
- caveats, especially run provenance, clipped labels, stale captions, or mismatch between figure plan and existing assets.

Also list deferred panels and explain why they were not generated in this iteration.

### 7. Validate

Run lightweight checks:

- confirm every saved panel is a valid image, for example with `file`;
- confirm the log exists and is readable;
- inspect a few key panels visually when judgment was involved;
- check git status and report only changes relevant to the request.

If a report dependency, archive, or source figure is missing, state that explicitly. Do not treat a missing input as a successful panel.

## Caption Guidance

When asked to write figure captions from selected panels:

- write from the actual saved panels and selection log, not from stale or aspirational caption drafts;
- distinguish model-implied quantities from measurements;
- state source/provenance succinctly when useful;
- avoid claiming evidence that the selected panel does not show;
- mark provisional panels when run provenance is not yet canonical;
- use LaTeX figure blocks only if requested, and omit `\includegraphics` when panels are not yet assembled into full figures.
