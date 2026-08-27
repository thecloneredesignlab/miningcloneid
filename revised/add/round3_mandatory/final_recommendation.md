# Round-3 scientific recommendation

## Overall decision

The integrated manuscript is suitable for internal scientific review and for preparing a reviewer response, but it is **not yet a GO for final high-impact-journal submission if the central claim remains a robust cross-context mechanism**. The current evidence supports a model-conditioned, hypothesis-generating comparison and several well-defined failure modes. It does not support a causal, anchor-invariant, regularization-invariant, or grid-independent oxygen--CIN--ploidy mechanism.

No revision can guarantee a journal impact factor. The appropriate target is a manuscript whose central statements are traceable, falsifiable, quantitatively stress-tested, and explicit about the limits of the data.

## Review requests judged scientifically necessary

1. **Causal and evidentiary language:** necessary and implemented. Direct observations, fitted states, numerical-search variability, and post-fit diagnostics are separated.
2. **O1/O2 endpoint mismatch:** necessary and implemented quantitatively. The shared fit assigns high-state probability 0.325; O2 has 16/20 high-state cells and an exact predictive model-check value of \(1.75\times10^{-5}\). The manuscript no longer treats access to both modes as reproduction of lineage-specific mixture weights.
3. **All-record adequacy:** necessary and implemented as an in-sample audit. The six-winner numerical range contains 14/92 nonzero burden observations and 4/8 terminal mean chromosome numbers. These ranges are not confidence intervals and are not held-out validation.
4. **Finite-time validation:** necessary and implemented. Overall day-1000 medians approach the dominant eigenvector, whereas the small-gap stratum retains median differences of 0.331 and 0.601 from 2N and 4N starts; asymptotic curve classes are therefore not universal finite-time predictions.
5. **Boundary, WGD, and necrosis audit:** necessary and implemented. Boundary-routed loss can contribute 26.7% of combined boundary plus within-grid CIN-associated loss for a subset of in-vivo response rows. Missing necrosis predictions were reconstructed and reproduce the reported objective. The WGD branch semantics are stated from the source implementation.
6. **Upper-grid sensitivity:** necessary and implemented after preserving the biological lower bound. Expanding only \(N_{\max}\) from 154 to 308 changes dominant mean ploidy by as much as 10.543 ploidy units and moves low-O2 modes near the new upper boundary. Low-O2 dominant modes cannot be called grid-independent biological equilibria.
7. **Single-anchor and saturated-penalty concern:** scientifically necessary, locally audited, but refit-dependent robustness remains incomplete. Eleven of 14 scalar directions agree across the current six winners, yet all share seed 10 and the median fraction of active coupling terms in the saturated Welsch region is 100%.

## Reasonable requests that cannot be completed from the present evidence alone

- Independent external-cohort validation and new wet-lab perturbations are scientifically valuable but require data that are not present in the repository. They should remain explicit future-validation requirements, not be simulated or replaced by an in-sample check.
- A guaranteed IF threshold is not scientifically decidable.
- A biologically stabilized high-chromosome mechanism may ultimately be needed to produce grid-independent low-O2 equilibria. That is a model-development question rather than a post hoc wording change.

## Required before a strong mechanistic final submission

1. Run the P0 alternative-anchor refits: six in-vivo warm starts crossed with culture anchors seed 132 and seed 157 under the current coupling specification.
2. Run the P1 less-saturating coupling refits for the six in-vivo warm starts with seed 10 under the specified \(c=1\) and \(c=10\) settings.
3. Implement and run the held-out folds specified in `hpc/heldout_validation_design_matrix.tsv`, keeping each held-out unit out of the fitting objective and evaluating it only after fitting.
4. Reassess whether the complete missegregation and survival-function directions persist across anchors, penalties, and held-out folds. If they do not, the central positive claim must be reduced to a descriptive conditional case study.
5. Repair the upstream necrosis exporter before final archival release, even though the reconstructed audit already verifies the current objective.

These optimization jobs have **not** been submitted. HPC submission requires a separate preview and user confirmation.

## Submission stance under two possible narratives

- **Strong cross-context mechanism narrative:** NO-GO until alternative-anchor, alternative-penalty, and held-out analyses are complete and concordant.
- **Transparent model-audit and hypothesis-generation narrative:** GO for internal circulation in the current integrated TeX, with Figure 6 explicitly retained as a truncated-grid failure-mode diagnostic rather than an equilibrium result.
