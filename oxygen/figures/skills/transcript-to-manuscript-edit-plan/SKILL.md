---
name: transcript-to-manuscript-edit-plan
description: Convert a transcript or free-form review comments about a manuscript into a line-mapped implementation plan. Use when Codex is given spoken/transcribed comments, reviewer notes, scientist feedback, or pasted annotations about a `.tex` manuscript or figure captions and must map each comment to manuscript lines, separate conceptual edits from wording edits, identify figure/caption changes, define shorthand terms succinctly, and prepare an actionable edit plan before modifying the manuscript.
---

# Transcript To Manuscript Edit Plan

## Overview

Use this skill to turn messy manuscript feedback into a concrete edit plan. The output should make clear where each comment applies, what action it implies, and which edits are safe to implement now versus dependent on missing figures or later Methods/Supplement text.

## Workflow

### 1. Locate The Manuscript And Context Files

Find the current manuscript source first. Prefer the file named by the user. If the user refers to "the `.tex` file" without a path, search likely locations with `rg --files -g '*.tex'` and inspect current candidates.

Also inspect nearby planning files when relevant:

- `oxygen/figures/iteration*/figureCaptions.txt`;
- `oxygen/figures/iteration*/panel_selection_log.md`;
- figure plans, caption drafts, claim notes, and evidence plans named by the user.

Use `nl -ba <file>` or targeted `rg -n` output so line references are stable.

### 2. Split The Transcript Into Atomic Comments

Break the transcript into small items. For each item, classify it as one of:

- manuscript text change;
- figure caption change;
- figure panel/backlog change;
- citation or figure-numbering correction;
- terminology/definition problem;
- duplicate or superseded comment;
- global style rule.

Keep the speaker's scientific intent, but clean obvious transcription artifacts only when the intended term is clear from context. If a correction could change the science, flag it instead of silently rewriting it.

### 3. Map Comments To Lines

For each atomic comment, identify the manuscript anchor:

- exact line number when the quoted text exists;
- section or figure block line when the comment applies to a paragraph/caption;
- `global` when it applies across the manuscript;
- `missing/new` when the requested content is absent and needs insertion.

When the transcript contains repeated passages, map only the latest or most specific version unless the repetition contains distinct instructions.

Use concise mapping tables when helpful:

```text
Transcript item | Manuscript anchor | Required action
```

### 4. Convert Feedback Into An Implementation Plan

Organize the plan in the order edits should be made, usually:

1. section titles and paragraph framing;
2. main-text scientific claims;
3. figure citations and numbering;
4. figure captions;
5. panel backlog or assembly changes;
6. terminology cleanup;
7. validation.

For each task, state:

- target file and line/section;
- specific change to make;
- reason from the transcript;
- dependencies, such as a missing panel or Methods section to be written later.

### 5. Apply Manuscript-Specific Guardrails

Use these defaults for this project unless the user says otherwise:

- do not include provenance phrases like "copied from an HTML report" in manuscript captions;
- keep missing panels present as manuscript markers when the user wants them to signal future figure generation;
- keep technical details succinct in Results and point to Methods/Supplement when detailed computation will be described later;
- replace shorthand terms such as "fixed-O2 analysis", "feature AUC", "mode", or "parameter landscape" with one or two sentences stating what was computed, what input was used, what output was produced, and why it matters biologically;
- distinguish measured oxygen from model-implied in vivo effective O2/resource stress;
- do not collapse biologically distinct claims just because they share a figure panel.

### 6. Validate The Plan Or Edits

If only drafting a plan:

- confirm every transcript item is either mapped, marked duplicate, or marked unmapped;
- confirm no requested comment has been silently dropped;
- identify open questions only when they block implementation.

If implementing edits:

- re-read the relevant `.tex` sections after editing;
- run targeted `rg -n` checks for stale figure numbers, stale phrases, and undefined shorthand;
- if figures were changed, run the figure assembly workflow and validate outputs;
- report any validation that could not be run.

### 7. Report Back

Return a short summary that separates:

- line-mapped transcript coverage;
- planned or completed manuscript edits;
- missing panels or evidence still needed;
- validation performed.
