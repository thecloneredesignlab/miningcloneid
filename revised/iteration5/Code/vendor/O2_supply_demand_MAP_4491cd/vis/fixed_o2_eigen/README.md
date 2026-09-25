# Fixed-O2 eigen-attractor visualization

## File registry

### `render_fixo2_eigen_attractor_figures.R`

- Recursively reads materialized coordinate/cluster CSV files under the result root.
- Draws embedding PDFs and PNGs without reading fit objects, best-parameter files,
  simulation code, or analysis code.
- Writes a figure manifest into the fixed-O2 eigen table contract.
- Called only by the fixed-O2 eigen runner.
