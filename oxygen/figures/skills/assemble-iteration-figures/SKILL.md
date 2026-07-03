---
name: assemble-iteration-figures
description: Assemble manuscript figure PNGs from an iteration panel folder, usually `oxygen/figures/iterationN`, using the repository assembly script. Use when Codex is asked to run or rerun figure assembly, add or remove a panel from an assembled figure, update `assembled_fig*.png` outputs, validate panel naming, or align assembled figures with `figureCaptions.txt`, panel-selection logs, or manuscript `.tex` figure blocks.
---

# Assemble Iteration Figures

## Overview

Use this skill to make assembled figure files reproducible after panel-level changes. The goal is to run the existing assembly path, validate what was assembled, and keep captions/logs synchronized with the actual panel set.

## Workflow

### 1. Identify Inputs And Outputs

Confirm the iteration folder and assembly script before running anything.

Default paths in this repository:

```text
oxygen/figures/iteration1/
oxygen/figures/assemble_iteration_panels.py
oxygen/figures/assembled_fig*.png
```

Use the requested iteration folder when the user names one. If no folder is named, infer it from the current figure work and state the assumption.

### 2. Inventory Panels

List panel images in the iteration folder. Treat filenames matching this pattern as assembly inputs:

```text
fig<figure-number><panel-letter>_<description>.png
```

Examples:

```text
fig3b_invitro_predicted_ploidy_distribution.png
fig4f_landscape_tsne_clusters.png
fig6a_pooled_fixed_o2_curve_class_examples.png
```

Compare the file list with `figureCaptions.txt`, `panel_selection_log.md`, and relevant `.tex` figure blocks when present. Missing planned panels are not automatically an error if the caption intentionally marks them as missing or "to be generated."

### 3. Preflight Guardrails

Before assembly:

- verify the assembly script exists;
- verify the iteration folder exists;
- identify newly added or removed panel files;
- check whether captions still reference stale output paths such as `oxygen/figures/iteration1/assembled_fig*.png` when the script now writes to `oxygen/figures/assembled_fig*.png`;
- preserve user changes and do not delete panel files unless explicitly asked.

If the requested source panel is missing but a report/archive source is known, restore it first using a report/panel curation workflow, then assemble.

### 4. Run Assembly

Run the repository script rather than manually composing figures:

```bash
python3 oxygen/figures/assemble_iteration_panels.py oxygen/figures/iteration1
```

Use the requested iteration folder in place of `iteration1` when applicable.

Read the script output. It should state which figures and panels were assembled and where outputs were written.

### 5. Validate Outputs

Do not skip validation silently. At minimum:

- confirm the expected assembled files exist;
- run `file` on changed `assembled_fig*.png` outputs;
- inspect layout-sensitive figures visually when panels were added, removed, or reordered;
- confirm the assembler output panel list matches the intended figure content;
- check `git status --short` for relevant changed/untracked files.

If ImageMagick or another optional image validator is unavailable, use the script's own validation output plus `file` and visual inspection.

### 6. Update Captions And Logs

When panel composition changes, update all applicable text:

- `oxygen/figures/iteration*/figureCaptions.txt`;
- `oxygen/figures/iteration*/panel_selection_log.md`;
- manuscript `.tex` figure captions or in-text panel references.

Caption what the assembled figure actually shows. Do not let aspirational captions claim evidence that is not present. If a panel is intentionally missing, keep the placeholder if the manuscript plan uses missing panels as work markers.

### 7. Report Back

Summarize:

- which panel files were added, removed, or restored;
- which assembled figure files were regenerated;
- the assembler output panel list;
- validation performed;
- the updated caption when the user asks for it.
