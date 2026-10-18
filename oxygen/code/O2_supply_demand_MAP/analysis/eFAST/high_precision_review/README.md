# Fixed-vector high-precision review

Point-level diagnostics for existing eFAST responses, with unchanged N/M and
parameter vectors. Original outputs are retained. No FAST samples are generated.

`review.py prepare` freezes all previously identified influential representative
points with seed138/101 ploidy >=6.5 and seed174 M8 phase5 response >1.01. Each
has a normal-response control from the same trajectory/phase/O2 and adjacent
existing O2-grid points at the original vector. Deduplication retains all roles
and relationships. This selection is not a census of all oxygen-grid anomalies.

`point.R` independently validates the original sample vector, config and code,
freezes the canonical double matrix, and saves the entire double eigensystem,
raw real/imaginary leading vector, sign-adjusted vector, clipped vector and
normalized vector. Scaled residuals are recorded at the relevant stages.

`perron_review.cpp` uses positive Noda inverse iteration on that same matrix at
50 and 100 decimal digits for both right and left vectors. No rates, graph edges
or matrix entries are perturbed. Internal working shifts are solver operations,
not modifications to the physical operator. Vectors/scalars are serialized as
decimal strings; tiny components are not first cast to double. Acceptance uses
scaled residuals, a leading-root upper bound and vector stability, with tighter
tolerances at 100 digits. Left/right agreement and vector precision differences
are independently checked using Python Decimal. The canonical matrix was built
in double precision, so this is not a 100-digit recomputation of all model rates.

Reduced-graph cases and possible nonunique dominant vectors are flagged. A
small residual alone is not a guarantee of a well-conditioned or unique vector.
The left/right overlap and eigenvalue condition estimate are recorded.

## Execution

1. Local: `python3 review.py prepare --root <result-root>` (NumPy required).
2. Commit/push current branch; fast-forward the HPC checkout.
3. HPC: `python3 submit.py smoke --root <result-root>`.
4. Wait for completed smoke; inspect fixtures and six anomaly/control results.
5. HPC: `python3 submit.py main --root <result-root>`.
6. Collect complete point directories and Slurm accounting; run
   `python3 review.py collect --root <result-root>`.

Slurm: qos=xxlarge, 12 hours, one CPU and 4 GB per point, no specified node and
no array-concurrency cap. The six smoke points are reused by the main collection.
Compiled code is shared through a smoke-built archive and copied to each task's
scratch. The submitted selection/code hashes must match the prepared snapshot.
The cache has a stable container bind target, and smoke explicitly reloads it
from a different scratch directory before main-array submission.

Each point has original-vs-high-precision results, full vectors, matrix/context
RDS, text matrix, residuals, graph information, timing, host, job and hashes.
Unconverged points are reported as unresolved; original-double disagreement
is a diagnostic finding rather than a workflow failure. `result.json` is the
final atomic point marker. Scientific classification is separate from Slurm
terminal state and collection completeness.

Python unit checks: `python3 -m unittest -v test_review.py` from this directory.
R fixtures in smoke test a known root/vector, preservation of a tiny component
and rejection of a non-Metzler matrix.
