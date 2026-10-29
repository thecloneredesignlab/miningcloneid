# Supplementary Figure8_7: in-vivo finite-time oxygen–ploidy matrix

Ten subplots show population mean ploidy at 100, 200, 300, 400, 500 days
(first row) and 600, 700, 800, 900, 1000 days (second row). Each subplot
uses the same 0–5% imposed-oxygen grid in steps of 0.1%. Solid blue, pink
and orange curves represent pure 2N, pure 4N and an equal-cell-number
2N+4N starting population. The displayed value is the median across all
1000 optimizer seeds (500 C01 and 500 C02); shading is Q25–Q75, not a CI.
No seeds are deduplicated or selected by survival. No Q labels are drawn.

Each seed retains its complete 26-parameter native in-vivo fitted vector,
including its own fitted p_misseg. Only imposed oxygen changes. Dynamics
are continuous growth without culture passages. A spectrally shifted
100-day matrix exponential propagates the three initial states through
ten consecutive intervals, normalizing state and retaining the population
factor between intervals. This is finite-time propagation, not a stationary
eigenvector or an assumption of equilibrium. The 1:1 mixture is propagated
directly; it is not an arithmetic mean of the pure terminal ploidies.

Pilot validation compares every snapshot with direct exp(tM) calculations.
Full-ensemble validation compares every day-1000 seed/start/oxygen value
with the existing native-p v8 results and its pooled summary with Figure8C.
The source chromosome grid and oxygen-unit limitations still apply.

Run: `figure8_time_snapshots_20261007_v12`. Source: authenticated
`figure8_fitted_family_20261006_v8`. The new matrix does not change Figure8C
or any other main/supplementary figure. No supp_fig8_8 was created in
that 2026-10-07 extension. On 2026-10-08, the user separately requested
family-resolved supplements: supp_fig8_8 for in vitro and supp_fig8_9 for
in vivo. Those regroup the same seed snapshots without modifying this
pooled figure.
