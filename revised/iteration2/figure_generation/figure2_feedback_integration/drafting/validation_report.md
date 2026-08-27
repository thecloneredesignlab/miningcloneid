# Figure 2 validation report

Status: **PASS with one documented wrapper limitation**

## Regeneration

- `scripts/agentRrunner.sh --check`: passed with R 4.5.1.
- Local revised source in worker mode: passed; generated PNG and PDF.
- Preserved prior source in worker mode: passed; regenerated the prior Figure 2
  for comparison.
- Canonical source in worker mode: passed; generated
  `revised/iteration2/Figures/assembled_fig2.{png,pdf}`.
- Canonical PNG and manuscript-staged PNG have the same SHA-256:
  `d34a97e34bbd36863dc9fa08cb412eac540c29c595acea5714fa6049f1688374`.

The convenience wrapper invocation without runtime arguments failed before
drawing because the shared runtime bootstrap requires an unrelated
`--invitro-result-root` option. Figure 2 has no fitting-data dependency, so the
documented worker mode was used for the applicable regeneration test. This
wrapper limitation is not reported as a successful test.

## Image and manuscript checks

- Review PNG: 2130 x 1905 pixels, RGB, non-interlaced.
- Direct visual inspection: no clipping, overlap, blank panel, or inconsistent
  panel label was observed. Lower- and higher-N curves differ by both color and
  line type; daughter rows are directly labeled.
- Canonical manuscript compilation: passed with `latexmk -pdf`.
- Iteration-3 review manuscript compilation: passed with `latexmk -pdf`.
- Rendered manuscript page 10 was inspected directly. The figure and complete
  caption are readable and remain on one page.

Document-wide float-size, table-width, and bibliography metadata warnings
remain outside this targeted revision. The Figure 2 float warning increased
because the more explicit caption is longer, but the rendered page was visually
checked and no content is clipped.
