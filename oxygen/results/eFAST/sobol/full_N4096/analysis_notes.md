# Complete independent in-vivo second-order Sobol analysis

One full-bound domain and all 500 clipped natural-span 10% neighborhoods are retained.
N4096, five independently scrambled repeats and the complete Figure 4 201-point grid.
Main indices use N4096; N2048 nested prefixes provide resolution diagnostics without extra model evaluations.
Repeat range <=0.05 and resolution delta <=0.02 are operational diagnostics, not proof of convergence.
No neighborhood is excluded for failing a diagnostic. Undefined indices remain undefined; endpoint medians
use defined values and report their counts explicitly. Negative estimates are sampling error, not antagonism.
S2 excludes main effects and interactions of order three or higher. S1/ST complement mechanistic interpretation;
ST cannot be summed into a process variance decomposition. Oxygen-band group summaries are group means
of parameter indices, not formal group Sobol indices. Endpoint fractions describe optimizer-solution robustness,
not biological confidence. Conditional response surfaces fix other parameters at the indicated optimum;
their nonadditive contrasts are not marginalized S2. Separate output/scope color scales must be read independently.

Source tables and diagnostic coverage: global/summaries and neighborhood10pct/summaries.
Individual 500-endpoint sources: neighborhood10pct/full/seed*/indices_summary.tsv.gz.
Main figures: global/figures and neighborhood10pct/figures. Per-endpoint figures remain separately available.
Raw float64 outputs, full repeated index arrays, checkpoints and exception matrices remain on HPC.
