# Independent in-vivo fixed-oxygen eFAST

The global analysis varies 14 parameters independently over the original in-vivo fit bounds, using the documented identity/log10 transforms. The 201 oxygen points span 0–5% by 0.025%.

Highest resolution: N=513; phase repetitions per cell: [5]. S1/ST have no sign. The existing Figure 4B Spearman panel supplies direction.

## Mechanism comparison

- dominant_mean_ploidy, S1: low-O2 death group mean 0.0229; low-O2 largest group buffering (0.0767). Missegregation 0.0677 -> 0.0298 and buffering 0.0767 -> 0.0447 from low to high O2. All three proposed conditions: not met in these point estimates.
- dominant_mean_ploidy, ST: low-O2 death group mean 0.3896; low-O2 largest group missegregation (0.5600). Missegregation 0.5600 -> 0.5149 and buffering 0.4126 -> 0.4225 from low to high O2. All three proposed conditions: not met in these point estimates.
- dominant_growth_rate, S1: low-O2 death group mean 0.0819; low-O2 largest group death (0.0819). Missegregation 0.0574 -> 0.0519 and buffering 0.0013 -> 0.0014 from low to high O2. All three proposed conditions: not met in these point estimates.
- dominant_growth_rate, ST: low-O2 death group mean 0.1319; low-O2 largest group missegregation (0.1561). Missegregation 0.1561 -> 0.1328 and buffering 0.0096 -> 0.0133 from low to high O2. All three proposed conditions: not met in these point estimates.

## Numerical diagnostics

- dominant_mean_ploidy, S1: p90 phase range 0.0704; p90 change from N=257 to N=513 0.0267.
- dominant_mean_ploidy, ST: p90 phase range 0.5209; p90 change from N=257 to N=513 0.2451.
- dominant_growth_rate, S1: p90 phase range 0.0824; p90 change from N=257 to N=513 0.0263.
- dominant_growth_rate, ST: p90 phase range 0.2054; p90 change from N=257 to N=513 0.0687.

Five repetitions measure phase variability; they do not establish resolution convergence. Fine ST rankings require caution wherever phase ranges or resolution changes remain large. Near-degenerate leading modes are separately recorded in spectral_gap_qc.tsv.

Group values are means of parameter indices, not additive group variance fractions. Total effects overlap through interactions. These indices are conditional on the chosen ranges and independent input distributions, not posterior uncertainty. The separate neighborhood10pct analysis examines robustness around each fitted endpoint.
