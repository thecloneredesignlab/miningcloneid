# Completed four-endpoint follow-up

All 280 main trajectories and both complete analysis chains finished on
2026-10-18. The 64 compact HPC source artifacts were collected and matched by
SHA256; an independent NPZ-to-table check reconciled 586,656 values/categories.
Six PDF pages were rendered and visually inspected; seven reporting tests passed.

The original README, submission_status.md and submission receipts are preserved
as historical submission snapshots. This file describes the completed review.

## Results

- [Detailed Chinese interpretation and recommendations](interpretation_zh.md)
- [Four-seed key-parameter diagnostics, all three O2 windows](ploidy_ST_key_parameter_details.tsv)
- [All 14 parameters](ploidy_ST_parameter_details.tsv)
- [Endpoint summary](seed_convergence_summary.tsv)
- [Collection receipt](results_collection_receipt.json)
- [Independent reconciliation](analysis_validation.json)
- [Delivery validation](deliverable_validation.json)

## Figures

| Experiment | ST curves | Numerical diagnostic overview | Growth ST curves |
|---|---|---|---|
| seed25/138, N4097 to N8193 at M4 | [PDF](N8193_M4/figures/ploidy_ST_comparison.pdf) | [PDF](N8193_M4/figures/ploidy_ST_diagnostic_overview.pdf) | [PDF](N8193_M4/figures/growth_ST_comparison.pdf) |
| seed101/174, M4 to M8 at N4097 | [PDF](N4097_M8/figures/ploidy_ST_comparison.pdf) | [PDF](N4097_M8/figures/ploidy_ST_diagnostic_overview.pdf) | [PDF](N4097_M8/figures/growth_ST_comparison.pdf) |

Each configuration also contains all-output diagnostics, phase and mechanism
summaries, measured runtime/projection tables, numerical recovery audits and
variance/gap diagnostics. The new plots are supplementary; original HPC PDFs
and PNGs retain their exact source hashes.

The N8193/M4 run did not establish uniform ploidy-ST convergence. M8 at fixed
N4097 was not uniformly more stable. The recommendations are diagnostic and
conditional; no new computation jobs were submitted during collection/review.

## Reproduction

See section 15 of the detailed interpretation. The frozen calculation sources
remain unchanged. The new collector and report are
`collect_resolution_nm_results.py` and `report_resolution_nm.py`, with
`test_resolution_nm_report.py`. Reporting requires NumPy>=2 and Matplotlib;
the standard-library collector can run with Python 3.6+.
