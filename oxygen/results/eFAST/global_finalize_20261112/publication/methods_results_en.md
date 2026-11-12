# Global variance-based sensitivity analysis

We analyzed the separate in-vivo model at 201 externally fixed oxygen levels (0-5%, 0.025% increments). Fourteen active inputs were varied independently over the documented fitting bounds, uniformly in encoded space (log-uniform on the natural scale for log10-encoded inputs). Outputs were mean ploidy from the normalized dominant eigenvector and asymptotic net live-cell growth from the leading eigenvalue real part.

The completed Sobol/Saltelli design used N=131072 and five independent scrambles, with second-order estimation enabled. It supplied S1, ST and all 91 unordered S2 pairs. Extended FAST used SALib 1.5.2 fast_sampler and its compute_orders estimator, verified against the public fast.analyze interface, at N=2049, 4097, 8193, 16385, M=4 and five independent phase repeats. Identical model and numerical-recovery rules were used. FAST bootstrap intervals were not used as convergence evidence.

Operational diagnostics were a five-repeat range <=0.05 and an adjacent-resolution mean change <=0.02. They were reported per output, parameter/index and oxygen level without excluding entries. The first new FAST resolution has no adjacent new-design diagnostic. Numerical convergence, nonzero interaction evidence, conditional response-surface direction and biological interpretation were treated separately. Endpoint neighborhoods remained archived.

## Results

All 47,838 Sobol index-by-output-by-oxygen entries met the operational diagnostics. Ploidy sensitivities showed stronger low-oxygen effects of death parameters and stronger high-oxygen effects of baseline missegregation and buffering. This supported parts of the proposed shift, while the complete missegregation group did not uniformly increase at high oxygen and death was not uniquely dominant over other processes. Net growth showed a different buffering pattern. Among all 91 pairs, lam_max with p_misseg had the largest peak S2 for both outputs (ploidy 0.10495 at 0% oxygen; growth 0.23646 at 0.075% oxygen). Rankings based on integrated oxygen bands were reported separately from peak rankings.

For Ploidy S1, eFAST N16385 versus Sobol N131072 had mean absolute difference 0.00248 and maximum absolute difference 0.01558.

For Ploidy ST, eFAST N16385 versus Sobol N131072 had mean absolute difference 0.08118 and maximum absolute difference 0.20238.

For Net live-cell growth rate S1, eFAST N16385 versus Sobol N131072 had mean absolute difference 0.00219 and maximum absolute difference 0.02207.

For Net live-cell growth rate ST, eFAST N16385 versus Sobol N131072 had mean absolute difference 0.01001 and maximum absolute difference 0.04615.


## Additional robustness findings

All approved FAST resolutions completed without filtering entries by convergence. At N=16385, 78/2814 ploidy ST entries and 2066/2814 growth ST entries met the operational diagnostics. Mean absolute S1 differences from Sobol were 0.00248 for ploidy and 0.00219 for growth, whereas mean absolute ploidy ST difference was 0.08118. At zero oxygen, O2_crit and n_O are structurally inactive: Sobol gave zero, while FAST mean ploidy ST estimates were 0.04259 and 0.03935. Thus repetition/resolution checks alone do not exclude spectral estimation error. Cross-method agreement for ploidy total effects has not been established; gamma_mu also had opposite oxygen-preference patterns. Primary mechanistic interpretation therefore relies on the completed Sobol analysis, with FAST retained as a diagnostic comparison. The evidence supports relative shifts in selected death/stress, baseline missegregation and buffering parameters, but not exclusive death dominance or a uniform high-oxygen increase of all missegregation parameters. Sensitivity magnitudes do not identify effect direction or parameter identifiability.
