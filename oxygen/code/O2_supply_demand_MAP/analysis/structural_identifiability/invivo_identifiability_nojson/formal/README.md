# Formal in vivo workflow

This package owns only the separate tumor run `20261111_invivo_identifiability_formal_010554Z`. The approved plan, authorization, integer role caps, scientific inputs and task registry are immutable under that run's `manifests/`. The closed pilot contributes 30491 core-seconds; the additional ceiling is 648000 and the cumulative ceiling is 691200. Unused role capacity does not authorize extra attempts in another role.

All numerical R work runs in the exact SHA256-pinned SIF on HPC, against the read-only five-file production numerical snapshot. The adapter evaluates each point in a fresh child and exports its predictions in a separate clean child with the same canonical score/cache preparation. The shared model, runner, image, libraries and culture assets are outside this package's write scope.

`workflow.py` owns submission receipts, immutable dispatches, coherent TSV generations and accounting from Slurm parents. `advance.py` makes one bounded pass and has no background loop. It submits primary tasks whose dependencies are ready; it never activates a reserve. A lost receipt keeps its reservation and blocks new submissions. Unexpected execution/source/score failures halt dependent work. Search inadequacy is a distinct quality outcome and can use only explicitly selected finite reserve identities.

`formal01` is the original immutable numerical execution release. The optional `formal02` bridge only permits the single approved reserve global to select a previously checked archive center or uniform initialization; it changes no model, score, box or optimizer controls. Its presence alone does not mean it was dispatched. The actual release, initialization population, starting vector, native result and every improvement are recorded per task. Consult the coherent job ledger for actual attempts.

The two free baselines each have twenty active coordinates. Independent hybrid and full-box uniform families must recover within 0.002, coordinate probes must not improve the best by more than 0.002, and the result cannot be worse than a known feasible point by more than 0.002. The compatible pilot points retain their original provenance and are feasible bounds only. Warm propagation of a winning point does not supply new independent recovery. Native optimizer code 0, exact score replay and search recovery are separate facts.

No density task is authorized by a passing engineering gate alone. Both score policies must pass the current baseline gate before the fourteen anchor searches and their nuisance refits can start. Seven positive density anchors are conditional diagnostics, not calibrated likelihood-ratio confidence intervals. An unresolved baseline leaves density profiles unlaunched, even when the study still has unused core-hours.

## Metadata and evidence tools

- `collect_closed.py`: copy only closed owned payloads, minimal runtime identity files, logs, accounting and start snapshots; verify payload hashes after transfer.
- `audit_records.R`: verify every full-precision improvement, independent score arithmetic, parameter/fixed-target constraints and stock/volume identities using archived arrays, with no model calls.
- `check_starts.R`: compare immutable first local callback records with their actual starting points and source scores.
- `postprocess.R`: evaluate engineering, baseline and profile gates from closed evidence. Gate files bind the selected payload hashes.
- `report.R`: export fit/parameter/role-accounting tables and the source-backed execution report. Bound flags use 1e-6 of each transformed span.
- `export_candidates.R`: export selected arrays, observation predictions, terminal chromosome distributions and acceptance labels. An unresolved candidate remains explicitly labeled.
- `plot_results.R`: create static PNG/PDF diagnostics from those tables, with no browser or simulator.
- `close_audit.py`: require all owned jobs closed and verify receipts, coherent state, source/release/dispatch/payload/image proofs, improvement audits, role caps and parent allocation accounting. Collect the final remote `state/` and `dispatch/` evidence before using it locally.

Local R metadata work uses the unchanged `scripts/agentRrunner.sh`. Metadata tools must not be used to execute new model evaluations locally. All new operational files are typed TSV/CSV, RDS, plain text, English Markdown or static PNG/PDF; no JSON transport or Java/JavaScript computation is part of this workflow.
