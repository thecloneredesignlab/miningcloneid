# Figure 2 drafting note

## Purpose and directive coverage

- FD02 addressed with one integrated implementation schematic.
- FD18 respected: the panel translates existing model code and introduces no
  fit or new model behavior.

## Generation and provenance

Local generator:
`agent-dev/manuscript_work_packages/ltee_hypoxia_model/drafting/scripts/make_figure2.R`.

Semantic source:
`oxygen/code/O2_supply_demand_MAP/model/model_O2_supply_demand_MAP.R`,
and its implemented kernel
`oxygen/code/O2_supply_demand_MAP/model/model_O2_supply_demand_MAP.cpp`,
especially the effective-resource update from simulated live-state demand,
chromosome-number- and oxygen-dependent growth/death functions,
multi-chromosome missegregation shift distribution, WGD transition,
buffering/survival kernel, and chromosome-state generator.

## Visual and semantic decisions

The five callouts follow the approved mapping: A combines context-specific
external input, initial state, and latent effective resource; B is growth and
death; C is missegregation/WGD; D is post-missegregation survival; and E is the
integrated output and fit comparison. Measured tumor burden and terminal
karyotype appear only as fit targets. The in-vivo resource feedback arrow
comes from the simulated live state, matching the implemented supply-demand
update. Missegregation is labeled `N -> N +/- delta N` because the kernel
supports a modeled multi-chromosome shift distribution rather than only
neighboring-state moves.

Each element is explicitly labeled as an external input, initial state,
modeled state, fitted function, state transition, predicted output, or observed
fit target. Experimental context uses the repository in-vitro sky blue and
in-vivo bluish green. Chromosome-state output uses the canonical 2N blue and 4N
vermilion among additional state colors. Solid and dashed arrows plus direct
text labels preserve meaning without color alone.

The small output histogram is schematic and does not encode a particular fit or
dataset.

## Visual QC

Inspected `final_figures/recommended/figure2.png` at its 7.1-inch draft width
after the semantic and minimum-text-size corrections. All five callouts,
arrows, formulas, type labels, feedback loop, and output labels are visible
without clipping or overlap. All annotations are at least 2.50 mm
(approximately 7.1 pt).

## Caveats

- The diagram communicates implemented dependencies, not proof that any fitted
  mechanism is necessary.
- Effective resource state is latent in tumors.
- Measured tumor burden constrains the fit but does not drive the modeled
  resource feedback; simulated live-state demand does.
- CIN/missegregation and WGD are separate modeled transitions.
- The survival curve is fitted rather than directly measured.
