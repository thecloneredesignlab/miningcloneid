# Chromosome state boundary audit: seed138 and seed101

Approved scope: 166 existing points (91 high-ploidy candidates + 75 ordinary
controls), N_MAX=154,220,308,440, N_MIN=22 and N_UNIT=22. Exact parameter vectors,
O2 and FAST N/M/phase identities are retained; no new FAST design or refitting.

Code: `oxygen/code/O2_supply_demand_MAP/analysis/eFAST/boundary_review/`.
The frozen case/task manifests contain all parameters and source hashes.
`selection.json` declares code, baseline, fitted-input and SIF provenance.

Execution: 24 smoke tasks (six cases at four ceilings), then 640 independent
main-array tasks, xxlarge/12h/1CPU/4GB, no specified node or array concurrency cap.
N_MAX=154 matrices are rebuilt exactly before reusing certified 50/100-digit
outputs. The other 498 tasks run both precisions and both vector sides.

Submission is not completion. See Slurm receipts/accounting, `completion.json`,
`collection_validation.json` and `independent_validation.json` for final status.
Results and original high-precision inputs are stored separately.

Scientific comparison uses all points including negative-growth states. This
selected-point audit diagnoses model-domain dependence and does not estimate
new ST values or certify the full eFAST design.
