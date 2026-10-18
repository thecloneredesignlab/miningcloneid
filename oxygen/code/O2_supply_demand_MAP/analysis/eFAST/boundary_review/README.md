# Chromosome state upper-bound review

The approved experiment retains all 91 high-ploidy candidates and 75 ordinary
controls for seed138/101 from the certified high-precision point review. It
changes only N_MAX to 154, 220, 308 or 440; N_MIN=22 and N_UNIT=22 remain fixed.
The 166 exact original natural parameter vectors/O2 point identities and their
FAST N/M, phase and source checksums are retained. No sampling or refitting.

There are 664 point/domain tasks. N_MAX=154 is independently reconstructed and
must match the saved matrix/config/parameters exactly; its certified 50/100-digit
outputs are reused byte-for-byte. The other 498 matrices receive full right/left
50/100-digit solves using the existing unmodified multiprecision solver.
Canonical rate/matrix construction remains double precision.

Boundary handling remains `drop`, with eps_tail=1e-8. Each matrix is rebuilt,
not padded. Shared-state rates and the common matrix block are compared with
the original. All original normalization stages and the double eigensystem are
retained. Complete decimal right/left vectors are saved, including tiny values.

The original canonical builder exposes boundary loss rates. A diagnostic-only
transition kernel with target/source ceiling 2*N_MAX captures all possible
daughters of included parents, separates upper/lower losses, and decomposes
upper losses into WGD/missegregation. This kernel is not used for eigenvalues
or ploidy. Weighted loss rates use the eigenvector of the actual analysis grid.
The WGD diagnostic uses the matching R proliferation-rate formula; its floating
roundoff is checked against the canonical total, without altering the operator.

## Execution

1. `python3 boundary.py prepare --root <result-root>` locally, after all code is saved.
2. Run unit checks, commit/push the current branch, and sync HPC with hashes.
3. `python3 submit.py smoke --root <result-root>` on RED.
4. Smoke runs six cases (boundary-high, interior-high and ordinary for each seed)
   at all four ceilings, known-answer fixtures, and compiled-cache relocation.
5. Only after smoke certification: `python3 submit.py main --root <result-root>`.
   The 24 smoke tasks are reused; the main array has 640 independent tasks.
6. Save Slurm terminal accounting and collection hashes, sync outputs locally,
   then `boundary.py collect`, `analyze.py verify`, `analyze.py figures`.
7. Run `mechanisms.py --root <result-root>` for diagonal/net-growth diagnostics
   and per-case domain classifications. This supplementary postprocessor is
   recorded separately from the frozen HPC evaluation code.
8. Run `finish.py collection --root <result-root>` to link every runtime job to
   completed Slurm accounting, then render/review all four PDFs and run
   `finish.py final --visual-render-reviewed --root <result-root>`.

Local postprocessing requires NumPy, Matplotlib and pypdf; the saved matrices and
decimal vectors suffice for repeating verification and plots. Repeating the
HPC evaluation additionally requires the original fitted-input and raw-sample
paths whose hashes are frozen in the manifests. Full parameter/config contexts
are also saved per task. Complete decimal source tables retain very small
components even when a plot uses a display floor (10^-120 in boundary-mass plots).

Every job uses qos=xxlarge, 12 hours, 1 CPU/4GB, no specified node and no array
concurrency cap. The original SALib1.5.2 SIF is hash-verified. Main tasks use
independent scratch caches restored from the smoke archive and bound to the
same stable container paths. Submission markers prevent accidental duplicates.

## Comparison and interpretation

Source tables report dominant mean ploidy, leading eigenvalue, full-distribution
TV distance (smaller grids are zero-padded without conditional renormalization),
peak state, mass at the ceiling, mass in the highest 22 states, mass above the
original N=154, and weighted boundary losses by mechanism. Negative growth is
retained and reported; it is not an exclusion rule.

Numerical domain-stability diagnostics compare N_MAX=308 to 440: absolute
ploidy change <0.01, growth change <1e-6 + 1e-3*max(abs(growth)), TV <1e-3,
and mass in the top 22 states <1e-4. All raw values are retained, and failure
is diagnostic rather than a reason to drop points. Strong connectivity,
left/right root consistency and precision agreement are checked separately.
Unresolved precision/reducibility is reported explicitly.

These selected points diagnose a possible model-domain source of unstable ST.
They do not constitute a new full FAST trajectory, do not estimate revised ST,
and cannot certify ST convergence or all 500 seeds. A peak that follows the
ceiling, or a material response change, requires reassessing the domain before
interpreting/increasing FAST sampling resolution.
