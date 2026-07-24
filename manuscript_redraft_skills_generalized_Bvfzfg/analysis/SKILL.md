---
name: analysis
description: "Plan, prepare, execute, validate, and audit analysis work, including model checks, simulations, data or image reprocessing, statistic-definition checks, provenance audits, substantial compute-heavy jobs, and analysis-scoped responses to user feedback. Use when Codex needs to run or design analyses, decide whether analytical support changes are needed, decompose instructions or feedback into analysis work packages, coordinate subagents for parallel analysis work, or submit analysis-scoped feedback responses."
---

# Analysis

## Purpose

Use this skill for analysis work that needs scientific and computational
discipline: planning, preparation, execution, validation, provenance, and
analysis-scoped feedback response.

This skill owns analytical work. It does not own manuscript routing, figure
drafting, Results prose, Methods prose, claim/evidence integration, or global
feedback closure. When user feedback or a feedback packet is provided, interpret
it only within the analysis scope and sign off only on that scope.

## Operating Rule

Pursue analysis work immediately until the next necessary step is a major compute
job or irreversible write. Do not launch expensive jobs, long local runs, broad
image processing, new model fitting, large-memory jobs, GPU jobs, or irreversible
dataset writes until the user has approved a brief analysis plan.

Lightweight discovery, input inspection, feasibility checks, feedback triage
within analysis scope, small summaries, and smoke tests are allowed before
approval when they reduce ambiguity.

If the user already supplied explicit approval for a specific plan, proceed with
execution and record that approval in the analysis notes.

## Inputs

Accept either or both:

- Direct user instructions, such as a requested model check, sensitivity
  analysis, statistic definition audit, simulation, reprocessing run, or
  provenance investigation.
- A served feedback packet from the feedback-manager helper, including raw
  feedback spans, feedback item IDs, a response template, and relevant artifacts.

When feedback is provided, read the cited raw spans and relevant artifact
context before deciding what the feedback means for analysis. Do not submit a
response unless the item ID, source spans, evidence artifacts, and served
feedback directory are clear enough for the feedback-manager helper to record a
scoped response.

## Workflow

1. Intake the request.
   - Read user instructions and any served feedback packet.
   - Identify the analytical questions, requested checks, relevant artifacts,
     expected evidence outputs, and any constraints from the caller.
   - Prefer maintained scripts, shared helpers, and canonical inputs over new
     code or ad hoc derived objects.

2. Decompose the work into packages.
   - Create analysis work packages with a question, direct inputs, intended
     outputs, validation checks, expected evidence artifact, and stop condition.
   - Keep packages independent when possible so they can run in parallel.
   - Mark package dependencies explicitly when one package needs another
     package's output.
   - Do not assign fixed package types up front; let the next necessary action
     determine whether each package can continue now or needs an approval gate.

3. Start independent packages in parallel with subagents.
   - Subagents are mandatory for multi-package analysis work. Use the available
     subagent or multi-agent tool to start independent packages in parallel.
   - The main agent should keep orchestration, shared decisions, approval gates,
     and final synthesis local.
   - Give each subagent a bounded package, explicit input paths, expected output
     artifacts, validation expectations, and a write scope.
   - Tell subagents that other agents may be working in the codebase and that
     they must not revert or overwrite unrelated changes.
   - If subagent tooling is unavailable for multi-package work, stop and ask the
     user whether to enable subagents or approve a single-agent degraded pass.

4. Pursue each package until it completes or reaches a major compute gate.
   - Continue with local inspection, small scripts, smoke tests, provenance
     tracing, and validation when they can run quickly and safely.
   - Stop a package when the next necessary action is a long walltime job,
     high-CPU job, large-RAM job, GPU job, broad image/data reprocessing, new
     model fitting, or irreversible dataset write.
   - For stopped packages, write a brief analysis plan and ask for approval
     before launching the major job.
   - A package may also finish with a documented no-change or not-applicable
     analysis decision if the evidence supports that outcome.

5. Create a brief human-readable plan when a major compute gate is reached.
   - Write the plan before implementation under
     `agent-dev/major_analyses/<YYYYMMDD_slug>/analysis_plan.md` unless the user
     or caller specifies another root.
   - Keep it short enough for real review.
   - Include the question or instruction, direct inputs, feedback inputs when
     present, output contract, reusable-code plan, compute/SLURM plan,
     validation plan, expected runtime class, and approval status.
   - End with explicit user options: approve, revise scope, defer, or ask for a
     feasibility-only pass.

6. After approval, implement narrowly and deploy appropriately.
   - Reuse existing functions from `R/`, `src/`, maintained `scripts/`, and
     pipeline utilities wherever possible.
   - Put broadly reusable functions in an appropriate shared utility module,
     with a short comment only when the intent is not obvious.
   - Keep entrypoint scripts thin: parse arguments, call shared helpers, write
     declared outputs, and record provenance.
   - Keep package-specific glue inside the analysis root only when it is not
     reusable.
   - Run R scripts through `scripts/agentRrunner.sh`.
   - Run small jobs directly from the terminal.
   - Use SLURM for parallel or batch work, GPU jobs, large-memory jobs, or
     expected runtime over about 5 minutes.
   - Use arrays or independent shards for embarrassingly parallel work; combine
     results with a deterministic reducer.
   - Check current QOS limits before full submission with
     `sacctmgr show qos format=Name,Priority,MaxTRESPU,MaxJobsPU,MaxSubmitJobsPU,GrpTRES -P`.
   - Choose the smallest QOS/resource request that comfortably fits the job, and
     scale arrays to minimize wall time without exceeding per-user or group
     limits.
   - Smoke test one small unit before submitting a full batch.

7. Validate, collect, and synthesize.
   - Run the promised validation plan, or the smallest validation appropriate
     for lightweight packages.
   - Write `validation_report.md` or `validation_report.json` under the analysis
     root when a package root exists.
   - Record commands, SLURM job ids, logs, package versions when relevant, git
     commit, important input checksums, and output checksums or manifests.
   - Merge subagent outputs into a concise synthesis that states which packages
     are complete, which are blocked by major compute or missing inputs, and
     which analysis outputs are valid for reuse.

8. Submit analysis-scoped feedback responses when applicable.
   - Get feedback through
     `.agents/skills/feedback-manager/scripts/serve_feedback --consumer_name analysis`
     unless the caller already supplied a served feedback packet.
   - Fill the served response template only after analysis artifacts changed,
     validation was completed, a partial analysis change was made, or a
     documented analysis-scope no-change/not-applicable decision was reached.
   - Include evidence artifacts for every change, validation, or no-change
     decision.
   - Submit the response with
     `.agents/skills/feedback-manager/scripts/feedback_signoff <served_feedback_dir> <response.json>`.
   - Do not edit feedback-manager JSONL ledgers directly or create global
     closure fields such as `validated_complete`, `global_resolution_status`, or
     equivalent.

## Analysis Scope

Analysis scope includes:

- analysis validity, bug checks, sensitivity checks, model reruns, statistic
  definitions, provenance audits, and explicit no-change decisions about
  analytical support;
- canonical data and model input/output integrity;
- reproducible execution of analysis scripts, jobs, reducers, and validators;
- evidence that a feedback item does or does not require an analysis change.

Analysis scope excludes:

- deciding whether manuscript feedback is globally addressed;
- routing work to manuscript stages or maintaining broader planning state;
- figure design, legend wording, Results/Methods prose, or claim/evidence graph
  editing, except to describe concrete analysis artifacts that changed.

## Subagent Package Contract

For each delegated package, provide:

- package name and analysis question;
- direct input paths and relevant feedback spans/items, if any;
- allowed write scope and output root;
- expected evidence artifact and validation checks;
- stop condition for major compute, missing inputs, or ambiguous scientific
  interpretation;
- instruction to report commands, changed files, output paths, validation
  status, blockers, and whether feedback signoff evidence is ready.

Do not ask subagents to perform broad manuscript routing or global feedback
closure. Keep their assignment inside analysis scope.

## Feedback Response

When working from served feedback:

- Treat raw user wording as authoritative. Read raw spans, not only served
  metadata.
- Interpret feedback only as far as needed to decide analysis work.
- Prefer concrete outcomes: artifact changed, validation completed, partial
  change, or no-change decision with evidence.
- If the feedback is not about analysis after inspection, submit a no-change or
  not-applicable response only when a served feedback packet and evidence
  artifact are available.
- Do not use feedback responses to pass vague interpretation work elsewhere.
  Record only the analysis-scope response and the concrete artifacts checked,
  changed, or used as evidence.

## Data Lineage

Prefer a simple active dataset path. A good analysis package depends on a
canonical source such as a maintained input object, maintained count table,
versioned model output, or explicit reproducible pipeline output.

Avoid depending on a long chain of derived objects, especially when provenance is
unclear. If a derived object is necessary:

- Trace and record its direct builder, input files, parameters, commit, and
  checksum.
- Prefer rebuilding it from canonical inputs inside the package or adding a
  small maintained builder.
- Treat dubious provenance as a blocker or user decision, not as a silent
  dependency.
- Keep new derived objects near their workflow family, with a manifest
  explaining how to regenerate them.

## Output Layout

Use one analysis root for substantial package coordination artifacts:

```text
agent-dev/major_analyses/<YYYYMMDD_slug>/
  analysis_plan.md
  contract.json
  run_manifest.json
  notes.md
  validation_report.md
  logs/
  slurm/
```

Place durable scientific outputs in established project locations when they are
meant to be reused:

- Tables and audit exports: `data/report_exports/<analysis_slug>/`
- Figures or review plots: `figures/<analysis_slug>/`
- Model run outputs: the maintained `data/runs/...` family when the analysis
  belongs to an existing model pipeline
- Job records: `slurm/runs/<analysis_slug>/` or the relevant existing SLURM run
  family

Avoid scattering outputs across ad hoc folders. If a location deviates from
these conventions, explain why in `notes.md`.

## Plan Contents

Use this minimal structure for `analysis_plan.md`:

```markdown
# <Analysis Title>

Status: tentative pending user approval

## Question Or Instruction
<One short paragraph describing what this analysis will decide or produce.>

## Work Packages
- <Package> - <status: continuing locally | complete | blocked pending approval>

## Direct Inputs
- <Canonical input path> - <why it is appropriate>
- <Any derived input> - <provenance status and rebuild plan>

## Feedback Inputs
- <Feedback item/span/served packet path, or "none"> - <analysis-scope relevance>

## Implementation Plan
- Reused code:
- New reusable helpers:
- Thin entrypoints:
- Output root:

## Compute Plan
- Runtime class:
- Local vs SLURM:
- Parallelization:
- QOS/resource request:
- Smoke test:

## Validation Plan
- Input checks:
- Smoke/sanity checks:
- Scientific checks:
- Reproducibility/provenance checks:

## Output Contract
- Tables:
- Figures/review artifacts:
- Manifests/logs:
- Feedback response file/signoff submission:

## Approval
Approve, revise, defer, or request feasibility-only review before execution.
```

## Contract Fields

Write `contract.json` before or during implementation when the package will
produce reusable outputs. Keep it compact:

```json
{
  "contract_type": "analysis",
  "schema_version": "0.2",
  "analysis_root": "agent-dev/major_analyses/<YYYYMMDD_slug>",
  "approval_status": "approved|pending|revised|deferred|not_required",
  "question_or_instruction": "",
  "work_packages": [],
  "blocked_major_compute_packages": [],
  "direct_inputs": [],
  "derived_inputs": [],
  "feedback_inputs": [],
  "entrypoints": [],
  "reusable_helpers": [],
  "output_roots": [],
  "compute_plan": {
    "uses_slurm": false,
    "qos": "",
    "parallelization": "",
    "smoke_test": ""
  },
  "validation_outputs": [],
  "feedback_records": []
}
```

## Validation Expectations

Select validation checks that match the risk of the analysis. At minimum,
substantial analyses should include:

- Input integrity: expected rows, keys, dimensions, experimental factors,
  timepoints or conditions, missingness, and checksum/provenance where relevant.
- Smoke execution: one small unit or reduced dataset before full deployment.
- Output integrity: declared files exist, are nonempty, have expected schemas,
  and can be loaded by downstream code.
- Scientific sanity: compare against an existing baseline, null model, prior
  maintained output, simulated expectation, or independently recomputed summary.
- Reproducibility: record commands, script paths, parameters, git state, SLURM
  ids, logs, seeds, and package/session information.

If validation fails, do not promote outputs for reuse. Record the blocker and the
smallest next action.

## Finish Criteria

An analysis task is complete only when:

- All work packages completed, were signed off as no-change/not-applicable, or
  are explicitly blocked pending user approval for major compute.
- Declared analysis outputs exist or blockers are documented.
- Validation passed or unresolved failures are clearly marked as blockers.
- Commands, inputs, parameters, git state, logs or job ids, and output manifests
  are recorded at a level appropriate to the task.
- Reusable code is documented by path and does not duplicate an existing helper.
- Applicable feedback responses are submitted only for analysis scope and only
  with evidence artifacts.
