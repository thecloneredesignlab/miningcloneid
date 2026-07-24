# Subpanel Semantic Interpretation Workflow

## Purpose

Interpret manuscript subpanels while keeping visual interpretation separate from
source clarification.

Visual claims must come only from a fresh-context visual inspector. Panel-scoped
source clarification may be added only after the visual read, by a separate
clarification subagent, and must stay separate from the visual interpretation.

## Roles

Choose the role before doing figure work:

- **Overseer:** default when asked to interpret multiple subpanels. Locate
  targets, choose output locations, launch one worker subagent per subpanel, and
  collect output paths. Do not inspect images, read rendered figures, or
  interpret panels.
- **Worker:** default for a single subpanel, or when an overseer assigns exactly
  one subpanel. Run the worker workflow for that assigned subpanel.

## Mandatory Subagent Requirement

Before applying either role, verify that a callable subagent or multi-agent tool
is available and permitted by active instructions.

If no subagent facility is available or subagent use is not permitted, stop and
tell the user this workflow cannot be applied because it requires delegated
fresh-context work. Do not substitute main-context image inspection,
self-contained guesses, or path-based interpretation.

## Overseer Workflow

1. Enumerate targets from user instructions, file listings, the canonical figure
   index, manifests, or other non-visual routing information. Do not open image viewers, inspect pixels,
   infer scientific meaning from filenames, or read figure content.
2. Choose output locations. Use a user-specified path if supplied; otherwise use
   a concise folder near the package or integration root, such as
   `subpanel_interpretations/`.
3. Launch workers. For each subpanel, launch a subagent with no forked
   conversation context. Tell it to use this workflow as a worker. Pass only the
   target image path, user-supplied figure/panel label, requested output path,
   and non-interpretive routing constraints.
4. Collect results. Verify each expected output file exists and follows the
   worker output contract. Do not assess whether the visual interpretation is
   right; check only completeness, target labels, paths, and obvious missing
   sections. Return worker output paths so integration can record them in the
   canonical figure index.

## Worker Workflow

1. Freeze the target.
   - Record the exact image path and supplied figure/panel label.
   - Use the path only as an identifier. Do not infer meaning from directory
     names, filenames, dates, package names, or prior conversation.

2. Fresh-context visual read.
   - Launch an image-inspection subagent with no forked conversation context.
   - Tell the image inspector not to use skills and not to launch further
     subagents.
   - Give it only the target rendered image, preferably as a local image item.
     Do not give surrounding paths except what is technically required to load
     the image.
   - Ask for visible pixels only: axes, axis labels, tick labels, legends visible
     inside the image, colors, marks, shapes, approximate positions, trends,
     annotations, panel letters, and visual caveats.
   - Prohibit scientific, workflow, data-lineage, or manuscript-level
     interpretations beyond what can be inferred from the image alone.
   - Preserve the image inspector output verbatim.

3. Panel-scoped source clarification.
   - Launch a second subagent only after the visual read is complete.
   - Tell the clarification subagent not to use skills and not to launch further
     subagents.
   - Give it the target path, supplied panel label, and fresh-context visual
     description.
   - Ask for exactly three clarification sections:
     - **Visual-observation checks:** quantify, qualify, or falsify trends,
       correlations, group differences, orderings, or other observations made by
       the image inspector, using existing plotted/source data where recoverable.
     - **Visible-element resolution:** clarify missing, ambiguous, or unresolved
       visible elements such as color/shape mappings, labels, units,
       denominators, transforms, facets, or abbreviations.
     - **Provenance chain:** trace plotted values as far as possible through
       existing plotting code, manifests, provenance files, generated tables,
       intermediate files, algorithms, transforms, filters, and upstream source
       datasets.
   - Prohibit raw user feedback, prior interpretations, manuscript drafts, broad
     manuscript context, rerunning analyses, new model fitting, or judging source
     data correctness.
   - Require explicit `not determined` entries for unresolved checks, elements,
     or provenance steps.

4. Extract visual claims.
   - State only claims supported by the image inspector's visual description.
   - Distinguish visible text/encoding from inferred data identity.
   - Use cautious language for approximate patterns, such as `appears`,
     `visually`, `roughly`, or `the visible points suggest`.
   - Do not use source clarification output as support for a visual claim.

5. Interpret visually, not evidentially.
   - Interpret the apparent visual message only from visible axes, labels, marks,
     annotations, and patterns.
   - Keep quantitative checks, provenance, source identity, transforms,
     denominators, and statistical summaries in the separate source-clarification
     section.
   - Do not make claims about workflow intent, drafting history, panel swaps,
     label correctness, or support for manuscript claims.

6. Name limits explicitly.
   - List visual ambiguities, unreadable labels, overplotting, missing
     definitions, and assumptions that cannot be checked visually.
   - State that source clarification was panel-scoped and that no raw feedback
     review, manuscript-level review, rerun analysis, or source-data correctness
     audit was performed.

## Output Contract

### Overseer Output

Return a compact coordination report with:

- **Role:** Overseer.
- **Targets:** subpanel labels and image paths assigned to workers.
- **Output Location:** folder or files where worker outputs were written.
- **Worker Results:** one line per subpanel with worker status and output path.
- **Not Assessed:** state that the overseer did not inspect figures or interpret
  panels.

### Worker Output

Write a compact interpretation artifact to the assigned output path and return
that path. The artifact must contain:

- **Target:** exact image path and supplied figure/panel label.
- **Role:** Worker.
- **Fresh-Context Visual Description:** verbatim image-inspector output or a
  clearly marked excerpt plus path to the full retained output.
- **Panel-Scoped Source Clarification:** verbatim clarification output or a
  concise extract with Visual-observation checks, Visible-element resolution,
  and Provenance chain.
- **Visual Claims:** concise claim-led bullets about visible axes, labels,
  encodings, marks, annotations, and apparent trends. Each claim must cite a
  short excerpt from the fresh visual description.
- **Visual-Only Interpretation:** cautious apparent message of the panel from
  pixels alone.
- **Clarification Notes:** terse bullets from source clarification only.
- **Visual Caveats:** visible limitations and ambiguities.
- **Not Assessed:** no raw feedback, prior interpretations, manuscript drafts,
  source-data correctness validation, rerun analysis, or manuscript-level
  evidence review.

Semantic coverage is canonical in `figure_set_manifest.csv`. A separate
semantic interpretation index may be generated as a derived review view, but is
not a required handoff object.

For visual claims, use:

```markdown
**V1. [Visual claim stated plainly.]**
Support: fresh visual inspector observed "[short quote or tight paraphrase]."
```

## Hard Rules

- If asked to interpret multiple subpanels, act as overseer by default.
- Overseers must not inspect figures, describe visual content, or interpret
  panels.
- Workers must not skip the fresh-context image-inspection subagent.
- Workers must wait for the image inspector before launching source
  clarification.
- Workers must tell subordinate subagents not to use skills.
- Workers must not directly inspect the image in the main context.
- Do not infer visual meaning from paths, filenames, package names, or prior
  conversation.
- Do not present provenance, code, legends, comments, filenames, or drafter
  summaries as visual evidence.
- Do not mix source clarification into visual claims or visual interpretation.
- Quantitative checks are allowed only for observations made by the fresh-context
  visual inspector and only from existing artifacts.
- Do not rerun analyses, fit new models, validate raw source-data correctness, or
  perform manuscript-level evidence review.
- Do not claim visible labels, legends, or annotations are correct; treat them as
  visible text only.
