# Implementation plan: write the default Gemcitabine live/dead model as a LaTeX Methods section

## Goal

Create a manuscript-ready LaTeX Methods section that accurately describes **the default model now maintained in the Gemcitabine-model repository under `Code/in-vitro/pkpd_live_dead_model/src/invitro_fitting.py`**, and corrects/expands the user-provided draft.

The existing draft text is outdated because it describes:
- extracellular Gemcitabine exposure \(C(t)\);
- an active-drug ODE \(dC_{\mathrm{dna}}/dt=\eta C(t)-k_{\mathrm{decay}}C_{\mathrm{dna}}\);
- shared \((\eta,k_{\mathrm{decay}},n,k_{\mathrm{tr}})\).

The current implementation does **not** use that model. The default implementation is driven directly by a ploidy-specific intracellular dFdCTP signal surface derived from PKPD data. Parent Gemcitabine PK is not used as the live/dead driver.

Write the Methods section for the **default model only**. Do not document alternative disabled models, unless a one-sentence parenthetical is necessary to explain that they were not used.

---

## Input files to inspect

Read the current code in the Gemcitabine-model repository:

```text
Code/in-vitro/pkpd_live_dead_model/src/invitro_fitting.py
```

Also inspect the latest fit-output TSVs if available:
- `joint_fit_summary.tsv`
- `effective_dose_scaling_summary.tsv`
- `baseline_death_summary.tsv`
- `confluence_death_summary.tsv`
- `optimizer_attempts.tsv`
- `pk_tail_fit_diagnostics.tsv`

Use the code as the source of truth for equations and symbols. Use the TSVs only for examples of fitted values or default run configuration if needed.

---

## Output file

Create a LaTeX file or Markdown file containing LaTeX, preferably:

```text
methods_default_gemcitabine_model.tex
```

or, if this repo prefers Markdown:

```text
methods_default_gemcitabine_model.md
```

The section should be directly pasteable into a manuscript.

---

## Required structure

Write the following sections:

```latex
\subsection*{Delay-aware cytostasis and cytotoxicity model for gemcitabine}

\paragraph{Intracellular dFdCTP signal surface.}

\paragraph{Effective dose-response correction.}

\paragraph{Immediate cytostasis and distributed delay to death.}

\paragraph{Alive/dead cell-count dynamics.}

\paragraph{Ploidy-specific partial pooling and calibration.}

\subsection*{Optimization routine}
```

You may adjust titles slightly for clarity, but keep the structure close to this.

---

## Critical corrections to the existing draft

### 1. Remove the extracellular Gemcitabine / \(\eta,k_{\mathrm{decay}}\) active-drug ODE

Do **not** write:

```latex
\frac{d C_{\mathrm{dna}}}{dt} = \eta C(t)-k_{\mathrm{decay}}C_{\mathrm{dna}}.
```

Do **not** claim that \(C(t)\) is extracellular Gemcitabine exposure.

Do **not** claim that \(\eta\) or \(k_{\mathrm{decay}}\) are fitted or used in the default live/dead model.

Instead, write that the model is driven directly by a ploidy-specific intracellular dFdCTP signal surface:

\[
C_{\mathrm{dFdCTP},p}^{\mathrm{PK}}(t,d)
\]

where:
- \(p \in \{2N,4N\}\) is ploidy;
- \(d\) is the administered Gemcitabine dose in \(\mu\mathrm{M}\);
- \(t\) is time in days;
- the signal is derived from the `dFdCTP (ng/mL)` PKPD workbook column after conversion to \(\mu\mathrm{M}\)-equivalent units and baseline subtraction.

### 2. Describe the PK-derived dFdCTP signal surface

Use notation similar to:

\[
C_{\mathrm{dFdCTP},p}^{\mathrm{PK}}(t,d)
\]

State that:
- the PK sheets provide calibrated dFdCTP profiles at documented reference doses;
- for each calibrated dose, the measured dFdCTP profile is baseline-subtracted;
- terminal behavior is modeled using an exponential tail when identifiable, otherwise a fallback half-life is used;
- the surface interpolates across dose in log-dose space within the calibrated range;
- below the minimum calibrated dose, the minimum-dose temporal profile is scaled proportionally with dose;
- at zero dose, the dFdCTP signal is zero.

Keep the prose concise. This is a Methods section, not a code walkthrough.

### 3. Define the effective dFdCTP signal used in the live/dead ODE

The default model applies both:
- a ploidy-specific power-law dose correction, \(\beta_p\);
- a shared Hill dose gate, with \(EC_{50}\) and Hill coefficient \(h\).

Use:

\[
\widetilde{C}_{p}(t,d)
=
C_{\mathrm{dFdCTP},p}^{\mathrm{PK}}(t,d)
\left(\frac{d}{d_{\mathrm{ref},p}}\right)^{\beta_p-1}
G(d;EC_{50},h,d_{\mathrm{ref},p}),
\]

where:

\[
H(d;EC_{50},h)=\frac{d^h}{EC_{50}^h+d^h},
\]

\[
G(d;EC_{50},h,d_{\mathrm{ref},p})
=
\frac{H(d;EC_{50},h)}
     {H(d_{\mathrm{ref},p};EC_{50},h)}.
\]

Define:
- \(d_{\mathrm{ref},p}\): minimum calibrated PK reference dose for ploidy \(p\).
- \(\beta_p\): ploidy-specific dose-scaling exponent.
- \(EC_{50}\): shared Hill-gate dose.
- \(h\): shared Hill coefficient.

Mention that the normalization by \(H(d_{\mathrm{ref},p})\) ensures the Hill gate equals one at the PK reference dose.

At \(d=0\), define:

\[
\widetilde{C}_{p}(t,0)=0.
\]

### 4. Preserve the transit-chain notation \(Z_i\), but drive it by \(\widetilde{C}_p(t,d)\), not \(C_{\mathrm{dna}}\)

Use:

\[
\frac{dZ_{1,p}}{dt}
=
k_{\mathrm{tr},p}\left(\widetilde{C}_{p}(t,d)-Z_{1,p}\right),
\]

\[
\frac{dZ_{i,p}}{dt}
=
k_{\mathrm{tr},p}\left(Z_{i-1,p}-Z_{i,p}\right),
\qquad i=2,\ldots,n.
\]

Define:
- \(n\): number of transit compartments; selected by model comparison/grid search over candidate integer values, not estimated as a continuous parameter.
- \(k_{\mathrm{tr},p}\): ploidy-specific transit rate.
- \(Z_{n,p}\): delayed effective drug signal used for death.

State that this is a linear-chain approximation that produces a gamma-distributed delay from effective dFdCTP exposure to death commitment.

### 5. Add immediate cytostasis

The default model is `immediate_cytostasis_delayed_death`.

Define the cytostasis multiplier:

\[
m_{\mathrm{cyto},p}(t,d)
=
\frac{1}{1+k_{\mathrm{cyto},p}\widetilde{C}_{p}(t,d)}.
\]

Define:
- \(k_{\mathrm{cyto},p}\): ploidy-specific cytostatic potency.
- \(m_{\mathrm{cyto},p}\in(0,1]\): bounded multiplier on proliferation.

Explain that cytostasis uses the immediate effective signal \(\widetilde{C}_p(t,d)\), whereas death uses the delayed signal \(Z_{n,p}(t)\).

### 6. Define drug-induced death

Use:

\[
\kappa_{\mathrm{Gem},p}(t,d)
=
k_{\mathrm{kill},p} Z_{n,p}(t).
\]

Define:
- \(k_{\mathrm{kill},p}\): ploidy-specific cytotoxic potency per unit delayed effective signal.

### 7. Add drug-independent baseline and confluence-associated death

The current default model also includes:
- constant baseline death, \(\mu_{\mathrm{base},p}\);
- confluence-dependent death, \(\mu_{\mathrm{conf},p}(A_p/K_p)^q\).

Use:

\[
\mu_{\mathrm{ctrl},p}(A_p)
=
\mu_{\mathrm{base},p}
+
\mu_{\mathrm{conf},p}
\left(\frac{A_p}{K_p}\right)^q.
\]

Define:
- \(A_p(t)\): live-cell count for ploidy \(p\).
- \(K_p\): ploidy-specific carrying capacity.
- \(\mu_{\mathrm{base},p}\): ploidy-specific low-density/background death hazard.
- \(\mu_{\mathrm{conf},p}\): ploidy-specific confluence-associated death hazard scale.
- \(q\): fixed confluence-death exponent; the default implementation uses \(q=4\).

Do not mention the implementation cap \(A/K\le2\) unless you include it as a numerical guard in one sentence. If included, phrase it as a numerical safeguard rather than a biological assumption.

The total death hazard is:

\[
\lambda_p(t,d)
=
\mu_{\mathrm{ctrl},p}(A_p(t))
+
\kappa_{\mathrm{Gem},p}(t,d).
\]

### 8. Define the live/dead ODE system

Use:

\[
\frac{dA_p}{dt}
=
r_p\,m_{\mathrm{cyto},p}(t,d)\,
A_p\left(1-\frac{A_p}{K_p}\right)
-
\lambda_p(t,d)A_p,
\]

\[
\frac{dD_p}{dt}
=
\lambda_p(t,d)A_p
-
k_{\mathrm{clear},p}D_p.
\]

Define:
- \(A_p(t)\): observed/live-cell compartment.
- \(D_p(t)\): observed dead-cell compartment.
- \(r_p\): ploidy-specific maximum proliferation rate.
- \(K_p\): ploidy-specific carrying capacity.
- \(k_{\mathrm{clear},p}\): ploidy-specific disappearance/clearance rate of observed dead objects.

Mention that initial conditions are replicate-specific:
\[
A_p(0)=A_{0,\mathrm{rep}}, \qquad D_p(0)=D_{0,\mathrm{rep}}.
\]

Do not introduce cell-cycle states. The model is explicitly “without cell-cycle state awareness.”

### 9. Remove or replace the old \(\Phi(C)\) discussion

The original draft embeds the kill term into a generic operator \(\Phi(C)\). If the surrounding manuscript still uses \(\Phi(C)\), keep a short bridging sentence:

\[
\Phi_{\mathrm{Gem},p}(t,d)=\lambda_p(t,d)
\]

or, if \(\Phi\) should include only drug-induced killing:

\[
\Phi_{\mathrm{Gem},p}(t,d)=\kappa_{\mathrm{Gem},p}(t,d).
\]

Ask the code/manuscript author to decide which convention is used elsewhere. Do not invent a diagonal state matrix unless the current manuscript requires it.

If retaining the operator form, make clear that in this in vitro live/dead fit the explicit ODE equations above are the fitted model.

### 10. Calibration and optimization section

Write a clear optimization routine section matching the implementation.

Include:

- Data: time-resolved live/dead object counts from IncuCyte assays across doses, replicates, and ploidies.
- Time unit: days.
- Dose unit: \(\mu\mathrm{M}\).
- Objective: negative-binomial count likelihood for alive and dead channels.
- Overdispersion parameters: \(\theta_A\) and \(\theta_D\) for alive and dead counts, respectively.
- Fitting: primary MAP/penalized likelihood fit with partial pooling across ploidies in log-parameter space.
- Positive parameters optimized on log scale.
- Ploidy-specific parameters: \(r_p, K_p, k_{\mathrm{tr},p}, k_{\mathrm{kill},p}, k_{\mathrm{clear},p}, k_{\mathrm{cyto},p}, \beta_p, \mu_{\mathrm{base},p}, \mu_{\mathrm{conf},p}\).
- Shared fitted Hill parameters: \(EC_{50}\), \(h\).
- Fixed structural choices: confluence exponent \(q=4\); candidate transit-chain length \(n\) selected by comparing fits over a grid of integer values.
- Optimizer: L-BFGS-B on a normalized objective for numerical conditioning, while reporting raw negative log likelihood / posterior objective.
- ODE solver: LSODA through `solve_ivp`.
- Model selection over \(n\): choose the \(n\) with lowest raw posterior objective / data NLL plus prior penalty, as implemented.

### 11. Negative-binomial likelihood

Include a compact expression.

For observed count \(y\) and model mean \(\mu\), with dispersion \(\theta\):

\[
Y \sim \mathrm{NB}(\mu,\theta),
\qquad
\mathrm{Var}(Y)=\mu+\frac{\mu^2}{\theta}.
\]

Then define total NLL over alive and dead observations:

\[
\mathcal{L}
=
-\sum_{\mathrm{obs}}
\log p_{\mathrm{NB}}\left(A_{\mathrm{obs}} \mid A_{\mathrm{model}},\theta_A\right)
-
\sum_{\mathrm{obs}}
\log p_{\mathrm{NB}}\left(D_{\mathrm{obs}} \mid D_{\mathrm{model}},\theta_D\right)
+
\mathcal{P}_{\mathrm{prior}}.
\]

Use \(\mathcal{P}_{\mathrm{prior}}\) for log-scale partial-pooling penalties.

### 12. Partial pooling description

Use language like:

Each positive ploidy-specific parameter \(\alpha_p\) is represented as:

\[
\log \alpha_p = \mu_\alpha + \delta_{\alpha,p},
\]

where \(\mu_\alpha\) is a shared population-level log parameter and \(\delta_{\alpha,p}\) is a ploidy-specific deviation penalized by a zero-centered Gaussian prior.

Do not claim Bayesian posterior sampling; the implementation is MAP/penalized optimization.

### 13. Symbol table

Include a compact LaTeX table listing symbols and implementation names.

Required mappings, expanded with units, brief definitions, and values from
`Data/in-vitro/pkpd_live_dead_model/invitro_fitting_outputs/bestFitSoFar_20260513T164159` in the Gemcitabine-model repository:

| LaTeX symbol | Code parameter | Unit | Brief definition | Best-fit/default value |
|---|---|---|---|---|
| \(A_p(t)\) | alive state | cell/object count | Live-cell count state for ploidy \(p\). | Replicate-specific initial value \(A_{0,\mathrm{rep}}\); model prediction varies over time. |
| \(D_p(t)\) | dead state | cell/object count | Observed dead-cell count state for ploidy \(p\). | Replicate-specific initial value \(D_{0,\mathrm{rep}}\); model prediction varies over time. |
| \(C_{\mathrm{dFdCTP},p}^{\mathrm{PK}}(t,d)\) | `DfdctpSignalSurface` | \(\mu\mathrm{M}\)-equivalent dFdCTP | PK-derived intracellular dFdCTP signal surface before beta/Hill correction. | Not a fitted scalar; derived from PKPD workbook. |
| \(\widetilde{C}_p(t,d)\) | effective corrected dFdCTP signal | \(\mu\mathrm{M}\)-equivalent effective dFdCTP | Effective dFdCTP signal after ploidy-specific beta correction and shared Hill gate. | Not a fitted scalar; computed from \(C_{\mathrm{dFdCTP},p}^{\mathrm{PK}}\), \(\beta_p\), \(EC_{50}\), and \(h\). |
| \(d_{\mathrm{ref},p}\) | `min_calibration_dose_uM`, `reference_dose_uM` | \(\mu\mathrm{M}\) | Minimum calibrated PK reference dose used to normalize beta/Hill dose corrections. | 2N: 0.1; 4N: 0.1. |
| \(\beta_p\) | `beta_dose` | dimensionless | Ploidy-specific power-law dose-scaling exponent. | 2N: 0.2533; 4N: 0.4877. |
| \(EC_{50}\) | `dose_gate_ec50_uM` | \(\mu\mathrm{M}\) | Shared Hill-gate dose at half-maximal gate activation. | 0.018529 \(\mu\mathrm{M}\) (18.529 nM). |
| \(h\) | `dose_gate_hill` | dimensionless | Shared Hill coefficient for dose gate. | 2.6928. |
| \(Z_{i,p}\) | transit compartments | \(\mu\mathrm{M}\)-equivalent effective dFdCTP | Transit-chain state \(i\) for delayed drug-death signal. | Dynamic state; not separately fitted. |
| \(n\) | `n_tr` | compartment count | Number of transit compartments selected by model comparison. | 5. |
| \(k_{\mathrm{tr},p}\) | `k_tr` | day\(^{-1}\) | Ploidy-specific transit-chain rate. | 2N: 6.8427; 4N: 2.8202. |
| \(k_{\mathrm{cyto},p}\) | `k_cyto` | \((\mu\mathrm{M}\)-equivalent dFdCTP)\(^{-1}\) | Ploidy-specific cytostatic potency on proliferation. | 2N: 8040.9143; 4N: 1163.1135. |
| \(k_{\mathrm{kill},p}\) | `k_kill` | day\(^{-1}\,(\mu\mathrm{M}\)-equivalent dFdCTP)\(^{-1}\) | Ploidy-specific cytotoxic potency per delayed effective dFdCTP signal. | 2N: 217.2640; 4N: 129.8923. |
| \(r_p\) | `r` | day\(^{-1}\) | Ploidy-specific maximum proliferation rate. | 2N: 1.2464; 4N: 1.7809. |
| \(K_p\) | `K` | cell/object count | Ploidy-specific carrying capacity. | 2N: 25,358.2284; 4N: 20,512.0674. |
| \(k_{\mathrm{clear},p}\) | `k_clear` | day\(^{-1}\) | Ploidy-specific clearance/disappearance rate for observed dead objects. | 2N: 0.1518; 4N: 0.0616. |
| \(\mu_{\mathrm{base},p}\) | `mu_base_death` | day\(^{-1}\) | Ploidy-specific drug-independent baseline death hazard. | 2N: 0.1168; 4N: 0.0174. |
| \(\mu_{\mathrm{conf},p}\) | `mu_confluence_death` | day\(^{-1}\) | Ploidy-specific confluence-associated death hazard scale at \(A_p/K_p=1\). | 2N: 0.0514; 4N: 0.0541. |
| \(q\) | `confluence_death_exponent` | dimensionless | Fixed exponent for confluence-associated death. | 4.0. |
| \(\theta_A\) | `theta_alive` | count-dispersion parameter | Negative-binomial overdispersion parameter for alive counts. | 6.9519. |
| \(\theta_D\) | `theta_dead` | count-dispersion parameter | Negative-binomial overdispersion parameter for dead counts. | 7.7583. |

Do **not** include \(\eta\) or \(k_{\mathrm{decay}}\) in the default model table.

---

## Style requirements

- Write in manuscript Methods style, not code-comment style.
- Use precise mathematical notation.
- Avoid implementation noise unless it changes the mathematical model.
- Mention implementation names only in the symbol table or in brief parentheticals.
- Do not include speculative biological interpretation.
- Do not include results, fitted values, or claims about 2N vs 4N unless explicitly requested.
- Do not mention old failed model versions except to avoid outdated equations.
- Keep the section self-contained.

---

## Validation checklist before finalizing

Before saving the Methods section, verify against code:

1. The preferred live/dead driver is dFdCTP, not parent Gemcitabine.
2. There is no active-drug ODE using \(\eta\) and \(k_{\mathrm{decay}}\).
3. The Hill gate is normalized at \(d_{\mathrm{ref},p}\).
4. The beta correction is applied before both cytostasis and delayed death.
5. Cytostasis uses immediate effective dFdCTP signal.
6. Death uses delayed transit-chain signal.
7. Baseline and confluence death are drug-independent and add to the death hazard.
8. Alive/dead ODEs use total death hazard in both alive loss and dead accumulation.
9. The optimization routine describes negative-binomial alive/dead count likelihood.
10. Partial pooling is described as MAP/penalized optimization, not MCMC.
11. The symbols map to actual implementation parameters wherever possible.

---

## Existing draft text to correct

Use the user-provided draft as a starting point only. Keep symbols only when they map exactly to existing implementation parameters. Specifically:
- Keep \(Z_i\), \(n\), \(k_{\mathrm{tr}}\), \(k_{\mathrm{kill}}\), and \(\Phi\) only if needed.
- Replace \(C(t)\), \(C_{\mathrm{dna}}\), \(\eta\), and \(k_{\mathrm{decay}}\) with the dFdCTP signal-surface notation.
- Expand the section to include cytostasis, beta/Hill dose correction, baseline death, confluence death, and the optimization routine.

---

## Deliverable

Return:
1. the new LaTeX Methods section;
2. a short bullet list of what was corrected relative to the original draft;
3. any unresolved ambiguity, especially whether the surrounding manuscript wants \(\Phi_{\mathrm{Gem}}\) to represent drug-induced death only or total death hazard.
