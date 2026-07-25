# LTee hypoxia manuscript handoff

## Branch and base

- Handoff branch: `manusctip_draft_ai`
- Parent branch: `hypoxia_ltee_figures`
- Parent/base commit:
  `5a1a4aa689aeaf2f726baa98261d9f4bf2cb0a5a`
- Initial handoff commits:
  - `38529fc57dc555472fdfd520c57510bbaeee1924` — Gate B drafting
    package, workflow skills, and support files.
  - `c9f126d992d333f5443c4e39d7f04cc2cba44653` — frozen
    manuscript-context copy.

This branch packages the completed Gate B figure-drafting state. It does not
claim to contain the scientist's newest parallel manuscript edits. The root
`ltee_hypoxia_model.tex` is deliberately not tracked at the branch tip. The
read-only Gate B manuscript snapshot is:

`agent-dev/manuscript_work_packages/ltee_hypoxia_model/drafting/source_context/ltee_hypoxia_model.tex`

Do not edit that frozen copy as the production manuscript. During integration,
the scientist's latest manuscript should remain the production source.

## Current Gate B status

The Gate B workflow is complete and its feedback state is signed off against
the approved Gate A architecture. This is a technical workflow milestone, not
scientist approval to polish, insert, or publish the figures. The scientist
must review `drafting/review_report.html` and explicitly approve the set for
production work.

The five-figure review set is:

- Figure 1: matched experimental design and motivating chromosome trajectories.
- Figure 2: integrated resource-state and chromosome-state model schematic.
- Figure 3: in-vitro fit adequacy and fitted/model-implied functions.
- Figure 4: in-vivo fit adequacy, latent effective oxygen, fixed-O2 separators,
  the preserved pooled in-vivo/in-vitro embedding, and all three formal in-vivo
  regions.
- Figure 5: joint-fit adequacy, all-six parameter ratios, and matched
  context-specific functions.
- Supplement: objective distributions and stored optimizer diagnostics.

The package contains six final PNG/PDF pairs, 98 frozen generator inputs,
file-level provenance, prior-code fidelity records, consolidated legends, and
a self-contained HTML review report. The final validator passed with all 31
drafted PNGs represented in the report and all 231 provenance records verified.

These are manuscript-ready review drafts, not a production manuscript
assembly. Fit checks are in-sample; they are not held-out validation or
inferential uncertainty.

## Start here

Read these files in order:

1. `MANUSCRIPT_HANDOFF.md`
2. `agent-dev/manuscript_work_packages/ltee_hypoxia_model/drafting/review_report.html`
3. `agent-dev/manuscript_work_packages/ltee_hypoxia_model/drafting/drafting_panels.md`
4. `agent-dev/manuscript_work_packages/ltee_hypoxia_model/drafting/validation_report.md`
5. `agent-dev/manuscript_work_packages/ltee_hypoxia_model/drafting/final_figures/recommended/legend.md`
6. `agent-dev/manuscript_work_packages/ltee_hypoxia_model/drafting/feedback_intake.md`
7. `agent-dev/manuscript_work_packages/ltee_hypoxia_model/drafting/feedback_manager_context.md`

For detailed provenance:

- `agent-dev/manuscript_work_packages/ltee_hypoxia_model/drafting/source_manifest.csv`
- `agent-dev/manuscript_work_packages/ltee_hypoxia_model/drafting/source_tables/frozen_input_manifest.csv`
- `agent-dev/manuscript_work_packages/ltee_hypoxia_model/drafting/prior_code_fidelity.csv`
- `agent-dev/manuscript_work_packages/ltee_hypoxia_model/drafting/prior_panel_disposition.csv`

The recommended figures are under:

`agent-dev/manuscript_work_packages/ltee_hypoxia_model/drafting/final_figures/recommended/`

## Validation

Run from the repository root:

```bash
scripts/agentRrunner.sh --check
scripts/agentRrunner.sh \
  agent-dev/manuscript_work_packages/ltee_hypoxia_model/drafting/scripts/validate_drafting_package.R
```

Expected package checks:

- six 300-dpi PNG/PDF pairs;
- one page per final PDF;
- every final artifact within 7.1 × 9.0 inches;
- 31 of 31 drafted PNGs represented in `review_report.html`;
- no pending feedback or fidelity status;
- all provenance hashes valid;
- no Figure 3E or Figure 6 final artifact.

The validator requires R package `png` and system commands `pdfinfo` and
`sha256sum`.

## Regeneration

The active figure generators use only committed frozen inputs, except Figure 2,
which uses the tracked model implementation as its semantic source.

Run in this order:

```bash
scripts/agentRrunner.sh \
  agent-dev/manuscript_work_packages/ltee_hypoxia_model/drafting/scripts/make_figure1.R
scripts/agentRrunner.sh \
  agent-dev/manuscript_work_packages/ltee_hypoxia_model/drafting/scripts/make_figure2.R
scripts/agentRrunner.sh \
  agent-dev/manuscript_work_packages/ltee_hypoxia_model/drafting/scripts/make_figure3.R
scripts/agentRrunner.sh \
  agent-dev/manuscript_work_packages/ltee_hypoxia_model/drafting/scripts/make_figure4.R
scripts/agentRrunner.sh \
  agent-dev/manuscript_work_packages/ltee_hypoxia_model/drafting/scripts/make_figure5.R
scripts/agentRrunner.sh \
  agent-dev/manuscript_work_packages/ltee_hypoxia_model/drafting/scripts/make_fit_diagnostics.R
scripts/agentRrunner.sh \
  agent-dev/manuscript_work_packages/ltee_hypoxia_model/drafting/scripts/write_source_manifest.R
scripts/agentRrunner.sh \
  agent-dev/manuscript_work_packages/ltee_hypoxia_model/drafting/scripts/make_ltee_hypoxia_model_review_report.R
scripts/agentRrunner.sh \
  agent-dev/manuscript_work_packages/ltee_hypoxia_model/drafting/scripts/validate_drafting_package.R
```

R packages used across the generators are `cluster`, `dplyr`, `ggplot2`,
`magick`, `patchwork`, `png`, `scales`, and `tidyr`. Report generation also
uses the system `base64` command. Package versions are not pinned in an
`renv.lock`.

Do not run `materialize_frozen_inputs.R` as part of a clean build. It is a
provenance-refresh utility that reads ignored upstream result stores. Refreshing
frozen inputs changes the evidentiary basis and requires explicit scientific
review and updated manifests.

## Skills to use

The exact workflow skills used for this handoff are committed under:

`manuscript_redraft_skills_generalized_Bvfzfg/`

These are repository-vendored skills, not installed under `.agents/skills/`,
so an agent may not discover them automatically. Explicitly ask the agent to
read the relevant local `SKILL.md` completely before acting. Where a bundled
example uses `.agents/skills/...`, substitute the vendored path above.

Use the workflow in this order:

0. **Integrate first.**
   Merge the final remote tip of this branch into a branch based on the
   scientist's checkpointed parallel work, following the instructions below.
   The scientist's current manuscript is authoritative; the frozen `.tex` is
   provenance only. Rerun the Gate B validator after the merge.
1. **Feedback manager — review and approval.**
   Use `feedback-manager/SKILL.md` while the scientist reviews
   `drafting/review_report.html`. The manager owns feedback state, routing, and
   evidence-backed signoff; it does not substitute for scientific
   interpretation. Never edit `feedback_tracking/*.jsonl` directly.
2. **Figure drafting module — only if corrections are requested.**
   Use `manuscript-figure-workflow/SKILL.md`, selecting its Figure drafting
   module. Revalidate the package and obtain feedback-manager signoff again.
   Do not jump from requested corrections directly to polishing.
3. **Figure polishing module — after two explicit inputs.**
   Proceed only after the scientist has approved Gate B and supplied the target
   journal's dimensions and formatting constraints. Use the Figure polishing
   module in `manuscript-figure-workflow/SKILL.md`. This is a new authorization
   point.
4. **Subpanel semantic interpretation module.**
   Use the corresponding module in `manuscript-figure-workflow/SKILL.md` for
   manuscript-visible evidence panels. The current Gate B package does not yet
   contain these interpretation sidecars.
5. **Figure set integration module.**
   After polishing and semantic interpretation are complete, use the Figure
   set integration module in `manuscript-figure-workflow/SKILL.md` to create
   the canonical normalized figure package. Production integration requires
   explicit approval.
6. **Legends, Methods provenance, and claim graph.**
   After the integrated figure package is stable, use
   `manuscript-legend-writing/SKILL.md` and
   `method-table-provenance/SKILL.md`. Use
   `claim-graph-integration/SKILL.md` if an integrated claim/evidence map is
   wanted; there is no current LTee starting claim graph.
7. **Results text.**
   Use `results-text/SKILL.md` after semantic interpretations are available.
   If a claim graph is used, refresh it before drafting Results.
8. **Mechanical insertion and manuscript assembly — last.**
   Once substantive manuscript edits are explicitly authorized, insert
   approved sidecars into the scientist's authoritative manuscript and use
   `manuscript-assembly/SKILL.md` to assemble and validate the production
   package. Component skills produce sidecars; assembly does not independently
   authorize scientific prose changes.

Use `analysis/SKILL.md` only when a new analysis, sensitivity check, model
rerun, or major computation has been explicitly requested and approved. Gate B
does not authorize such work.

Separate authorization is required for journal-specific polishing,
figure-set/production integration, substantive manuscript edits, new or
long-running analysis, necrosis repair, Figure 3E, and Figure 6.

The existing Gate A/Gate B tracking root is:

`agent-dev/manuscript_work_packages/ltee_hypoxia_model/feedback_tracking/20260724_gate_a_approval`

When new feedback exists, serve it through the committed feedback-manager
scripts and return changed artifacts and validation evidence through the same
mechanism.

## Exclusions and unresolved work

The following were intentionally not done:

- no Figure 3E or restricted-model/negative-control ablation;
- no Figure 6 or its supplementary analytical grid;
- no optimizer refit, alternate-seed search, new parameter grid, or new
  embedding;
- no repair or reconstruction of the missing necrosis predictions;
- no production manuscript edit or figure insertion;
- no journal-specific polishing;
- no held-out validation, posterior interval, bootstrap interval, or
  identifiability analysis.

Open decisions for the scientist or user:

- approve or revise the Gate B figures after reviewing `review_report.html`;
- select journal-specific dimensions and formatting constraints;
- decide whether any excluded analysis is newly authorized;
- reconcile the scientist's parallel manuscript edits before production
  assembly;
- decide when the polished figure set is ready for legend, Results, Methods,
  and assembly workflows.

## Integration with parallel parent-branch work

The remote parent was still at the recorded base commit when this handoff was
created. Parallel scientist work may therefore be local or unpushed and must be
checkpointed before integration.

1. **Protect the scientist's work first.**
   Inspect `git status` and commit intended changes on a dedicated scientist
   branch. Do not use a broad `git add -A` in a dirty worktree, `git stash
   --all`, or `git clean`. Because the repository ignores `*.tex`, either copy
   the current manuscript outside the repository and record its SHA-256, or
   deliberately `git add -f` and commit it on the scientist checkpoint branch.
2. **Use a clean integration worktree.**
   Fetch the remote branches, then create a new integration branch from the
   scientist's latest committed tip. A separate worktree avoids collisions with
   unrelated untracked files.
3. **Verify ancestry before merging.**
   Confirm that both the scientist tip and `origin/manusctip_draft_ai` descend
   from `5a1a4aa689aeaf2f726baa98261d9f4bf2cb0a5a`. Stop and reassess if either
   ancestry check fails.
4. **Merge the final branch tip, not individual commits.**
   Merge `origin/manusctip_draft_ai` into the scientist-based integration
   branch. Do not cherry-pick `38529fc` and `c9f126d` in sequence: the
   intermediate history tracks a root `.tex` and can collide with the
   scientist's ignored manuscript, whereas the final handoff tip has no net
   root-TeX addition relative to the base.
5. **Review overlapping scientific sources.**
   Most handoff paths are new and package-local. If parallel work changed the
   model R/C++ sources used by Figure 2, inspect whether the schematic remains
   semantically accurate before accepting or updating manifest hashes.
   Preserve the scientist's scientific, code, and manuscript changes; preserve
   the frozen Gate B package unless deliberately revising it. Compare any
   overlapping vendored-skill changes manually.
6. **Restore and compare the authoritative manuscript.**
   If it was copied outside the repository, restore it after the merge. Compare
   it with the frozen snapshot using `git diff --no-index`; exit status 1 is
   expected when the scientist's manuscript is newer or otherwise different.
   Historical feedback may mention the former root `.tex`; do not rewrite
   archived feedback or provenance records for that reason.
7. **Validate after resolution.**
   Run the package validator above. If figure generators, frozen inputs, or
   semantic source code changed, regenerate the affected figures and provenance
   deliberately before accepting the integration.

One safe command outline is:

```bash
# Run these from the scientist's checkpointed branch.
git fetch origin
git rev-parse HEAD
git branch backup/20260725-pre-manuscript-integration HEAD
git merge-base --is-ancestor \
  5a1a4aa689aeaf2f726baa98261d9f4bf2cb0a5a HEAD
git merge-base --is-ancestor \
  5a1a4aa689aeaf2f726baa98261d9f4bf2cb0a5a \
  origin/manusctip_draft_ai
git worktree add ../miningcloneid-manuscript-integration \
  -b integrate/manuscript-draft HEAD
cd ../miningcloneid-manuscript-integration
git merge --no-ff origin/manusctip_draft_ai

# Restore the scientist's preserved .tex here if it was copied externally.
git diff --no-index -- \
  agent-dev/manuscript_work_packages/ltee_hypoxia_model/drafting/source_context/ltee_hypoxia_model.tex \
  ltee_hypoxia_model.tex

scripts/agentRrunner.sh \
  agent-dev/manuscript_work_packages/ltee_hypoxia_model/drafting/scripts/validate_drafting_package.R
```

The post-merge validator should again report six PNG/PDF pairs, all 31 drafted
PNGs represented, and all 231 provenance records verified. If Git reports that
an untracked file would be overwritten, do not delete it: preserve it, then
retry in a clean worktree. Do not perform the merge until the scientist's
parallel work, especially the current manuscript, has a durable checkpoint.
