# Fixed-O2 finite-time validation

All 17 source-validation checks passed. The comparison contains 100500 seed-O2 pairs per time point (500 seeds x 201 O2 values), both 2N and 4N starts, and 13 horizons through day 1000.
At day 1000 the overall median absolute difference from the dominant-eigenvector ploidy is 1.3e-07 for 2N starts and 2.87e-07 for 4N starts.
This convergence is not uniform: in the spectral-gap <0.005 stratum, the corresponding medians are 0.331 and 0.601, whereas the >=0.01 stratum is essentially converged.
At day 1000, finite-time curve classes match the eigen class for 72.2% of 2N starts and 70.6% of 4N starts; the two finite-time initial conditions agree with each other for 96.2% of seeds.

Interpretation: dominant eigenvectors are valid steady-state attractor diagnostics for well-separated spectra, but they are not universal finite-time predictions. Small-gap results and curve-shape labels must be qualified or excluded from strong interpretation.
