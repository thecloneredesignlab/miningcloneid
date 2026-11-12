# Global eFAST refinement: submission status

Snapshot UTC: 2026-11-12T17:10:14.514580Z

Preparation job 21842774 completed successfully after its packaging dependency was restored and the same job requeued. The original failure is retained under validation/preparation_initial_failure/. No model array was resubmitted. The first full-oxygen production smoke block is running.

All three requested levels are submitted with sequential dependencies, five phases, 201 oxygen points, M4, xxlarge, 12h initially, no node constraint and no array throttle. All twelve scoped monitors were running at this snapshot; they add 6h when fewer than 6h remain and renew themselves as required.

| Stage | Job | Array tasks | Queue / terminal state | Extension monitor |
|---|---:|---:|---|---:|
| prepare | 21842774 | 0 | COMPLETED | 21842786 |
| smoke | 21842775 | 0 | RUNNING: 1 | 21842789 |
| model_32769 | 21842776 | 2310 | PENDING: 2310 | 21842790 |
| analyze_32769 | 21842777 | 70 | PENDING: 70 | 21842791 |
| report_32769 | 21842778 | 0 | PENDING: 1 | 21842792 |
| model_65537 | 21842779 | 4550 | PENDING: 4550 | 21842793 |
| analyze_65537 | 21842780 | 70 | PENDING: 70 | 21842794 |
| report_65537 | 21842781 | 0 | PENDING: 1 | 21842795 |
| model_131073 | 21842782 | 9030 | PENDING: 9030 | 21842796 |
| analyze_131073 | 21842783 | 70 | PENDING: 70 | 21842797 |
| report_131073 | 21842784 | 0 | PENDING: 1 | 21842798 |
| final | 21842785 | 0 | PENDING: 1 | 21842799 |

The exact cached reference point in the first block was validated and reused across all 201 oxygen levels. Completed blocks/checkpoint chunks are retained for recovery.

Sobol N131072 and archived neighborhoods are untouched. New FAST scientific indices have not yet completed; synthetic test fixtures are not scientific results.

- [Full launch receipt](launch_receipt.json)
- [Status validation](validation/submission_status.json)
- [Monitor ledger](time_extension/submission_jobs.tsv)
- [Preparation recovery](preparation_recovery_receipt.json)
