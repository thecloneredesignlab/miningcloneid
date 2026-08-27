# Stage 4 regression: functionally organized simulation

Stage 4 organized fitted-parameter numerical production by biological or
workflow function. Simulation owns concrete quantities required by figures and
downstream analysis; it does not own statistical comparison, visualization, or
report presentation.

## Standard in-vivo domains

The standard producer under `simulation/invivo/` delegates table production to
domain modules for:

- `o2/`
- `population/`
- `ploidy/`
- `cin/`
- `functional_response/`

For real completed seed366, all 31 legacy scientific TSV files matched exactly.
The largest checked table contained 1,065,064 rows and the maximum numerical
difference was zero. The downstream pure visualization produced all 28
expected figures; all four PNG files were byte-identical and 20 of 24 PDF first
pages were pixel-identical. The other four differed only in the corrected
manifest-derived fit-directory title.

## Standard in-vitro domains

The in-vitro producer under `simulation/invitro/` delegates to `o2/`,
`population/`, `ploidy/`, and `cin/` domain modules.

For the frozen seed10 fixture:

- all ten scientific tables matched at tolerance `1e-12`;
- all eleven fixed-seed PNG files were byte-identical;
- visualization succeeded without `fit_result.rds` or `best_params.tsv` in its
  input fixture, demonstrating that numerical reconstruction had moved out of
  the visualization layer.

## Specialized simulation families

- Mixed-ploidy perturbation reproduced all five legacy scientific files,
  including the compressed trajectory table, byte for byte.
- Factorial interaction reproduced all four untreated scientific files byte
  for byte.
- Process-fingerprint simulations now write numerical process trajectories and
  O2/ploidy event inputs for downstream analyses.
- Live-effective-`p_ms` materialization is organized under
  `simulation/invivo/cin/` and reproduced the real seed366 numerical fixture.
- Fit-result prediction materialization lives under `simulation/fit_results/`.
- Parameter-landscape numerical products live under
  `simulation/parameter_landscape/`.

## Fixed-O2 regression

The fixed-O2 workflow was separated into canonical simulation, analysis,
visualization, report, and runner layers. Against the real 500-seed fixture:

- all 15 attractor tables matched the legacy products;
- all eight expected PDFs were produced;
- the composite first-page raster SHA-256 matched the frozen reference;
- the consume-only HTML report was 9,542,096 bytes and embedded 28 images;
- all 200 functions exposed by historical fixed-O2 entrypoints remained
  available with unchanged formal arguments.

## Stage status

The component comparisons above passed. The final aggregate parse, unit,
architecture, immutable-core, and whitespace gate remains pending until all
parallel structural work is complete:

```bash
Rscript oxygen/tests/run_o2_reorganization_regression.R
```
