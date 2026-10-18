# Fixed-vector high-precision review

This is a point-level audit of existing separate in-vivo eFAST responses, using
the exact original 14-dimensional parameter vectors and original model config.
No FAST designs, N/M changes or changes to original sensitivity outputs occur.

## Frozen scope

401 unique matrix evaluations:

* 106 anomalous candidates: seed138 11, seed101 80 (ploidy >=6.5), and seed174 15
  (M8 phase5 influential responses above 1.01).
* 86 deduplicated ordinary controls from the same trajectory/phase/O2.
* 209 adjacent existing oxygen-grid points at the exact anomaly vector.

Candidates are all qualifying points from the previously disclosed top-three
influential-response table at representative O2 cases. This is not a census of
all anomalous responses at all O2 values and not a population estimate for 500
fit endpoints. `relationships.tsv` preserves anomaly/control/neighbor links.

Each point preserves the canonical 133x133 double matrix (chromosome states
22–154), its exact RDS representation, the entire double eigensystem, all
normalization stages, and 50/100-digit right and left vectors as decimal strings.
Multiprecision applies to eigensolves of the fixed double matrix; original model
rates are not recomputed at 100 digits. No epsilon edges or rate regularization
are introduced.

## Execution and provenance

All final jobs used qos=xxlarge, 12 hours, one CPU, 4 GB, no specified node and
no array-concurrency cap. SIF:
`o2_supply_demand_map_r442_hpc_exact_salib152_20260923.sif`, SHA-256
`0b60f6cdab8a91f6bbeea1ab6cf02dd79660a295f24bc6e8ad6c5fed04d982ad`.

* Initial code/selection commit: `a2c87d8`.
* `21139700`: known tiny-component fixture exposed premature convergence of an
  extremely small component. No scientific points ran. Componentwise convergence
  was added at `511cb1e7`.
* `21139701`: six-point scientific smoke passed; its complete outputs are kept
  in `history/smoke_21139701/`.
* `21139702`: initial 395-task array failed at cache loading before point
  evaluation. A temporary path embedded in the compiled cache was corrected at
  `9148617d`, using a stable container bind target and relocated-cache fixtures.
* Final smoke `21140097`: six points and cache-relocation fixtures passed.
* Final main array `21140098`: remaining 395 points completed with exit 0:0.

Historical selections, receipts and logs are retained. A compiled cache is a
rebuildable runtime artifact and is excluded from Git/download; its hash and
location are retained in `remote_collection_receipt.json`.

## Independent validation

The worker checks original sample/config/code identities. Its multiprecision
solver checks positivity, scaled residuals, a leading-root bound, L1 stability
and relative componentwise stability. The collector compares 50/100-digit
vectors and ploidy, and left/right eigenvalues using Python Decimal.

Local validation separately reconstructs the binary-double matrix values from
round-trip text via `Decimal.from_float`, then recomputes matrix-vector residuals
and leading-root upper bounds at 160 decimal digits. It verifies the full raw
file hashes and links every point to final Slurm accounting. Previously generated
phase/spectrum artifacts are checked against their immutable manifest.

## Reproduce local analysis

From the repository root, with NumPy, Matplotlib and pypdf:

```bash
CODE=oxygen/code/O2_supply_demand_MAP/analysis/eFAST/high_precision_review
ROOT=oxygen/results/eFAST/neighborhood10pct/pilot/N8193_M4_and_N4097_M8/phase_spectrum_diagnostics/high_precision_review
python3 "$CODE/finish.py" collection --root "$ROOT"
python3 "$CODE/review.py" collect --root "$ROOT"
python3 "$CODE/analyze.py" verify --root "$ROOT"
python3 "$CODE/analyze.py" figures --root "$ROOT"
python3 "$CODE/finish.py" final --root "$ROOT"
```

Render all four PDFs with Poppler and inspect them before passing
`--visual-render-reviewed` to the final command. The delivered snapshot records
that completed inspection. The large input cache used to select points remains
in the parent diagnostic's ignored `raw_cache/` and `processed_cache/`; all exact
point parameters are versioned here, and the original input paths and hashes
are retained.

## Reading and files

* `interpretation_zh.md`: detailed findings, scientific limits and recommendations.
* `point_evidence.tsv`, `group_summary.tsv`: independent point/group summaries.
* `independent_residuals.tsv`: recomputed residuals and root bounds.
* `full_vector_comparison.tsv`: all 133 right-vector components per point,
  including exact 100-digit strings and decimal log10 values.
* `points/point_*/`: complete matrices, contexts, eigensystems, vectors and hashes.
* `figures/`: original-versus-high-precision ploidy, complete vectors, oxygen
  neighbors and residual/normalization diagnostics, with PDF/PNG pairs.
* `slurm/`: actual submissions, resource settings, final accounting and job logs.
* `history/`: preserved prior smoke and selection revisions.
* Collection/independent/deliverable validation JSON and artifact manifest: QA.

`high_ploidy_retained` classifies all points with original ploidy >=6.5 and a
verified difference <=1e-6, including O2 neighbors. It is not the count of the
106 selected anomalous candidates. Numerical retention does not establish a
biological mechanism or convergence of the complete eFAST ST estimate.
