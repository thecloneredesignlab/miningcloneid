# Local global Sobol convergence comparison

This directory contains derivative analysis of completed global indices only.
No model evaluations or sensitivity re-estimation are performed. Neighborhoods
remain archived. The original N65536/N131072 exports are preserved unchanged.

## Reproduce from the repository root

Use Python 3.12 with the versions in `requirements.txt` (the pinned local plotting
environment, not a change to the HPC SIF):

```sh
python3 -m pip install -r oxygen/results/eFAST/sobol/global_convergence_20261021/comparison_N65536_N131072/requirements.txt
python3 oxygen/code/O2_supply_demand_MAP/analysis/eFAST/sobol/compare_global_convergence_exports.py --root oxygen/results/eFAST/sobol/global_convergence_20261021 --collection oxygen/results/eFAST/sobol/global_convergence_20261021/collection_20261029T194045Z
```

Add `--render-pdfs` to render every original/comparison PDF with Poppler and
regenerate the validation contact sheets; this requires `pdftoppm` on PATH.
Manual visual-review notes are separate from the automated rendering record.

For transfer, use `collect_global_convergence_exports.py` in the same code folder
with `--remote-root` and `--local-root`. It selects only completed exports and
small provenance files, verifies source/destination hashes, and preserves
different existing files by refusing to overwrite them.

```sh
python3 oxygen/code/O2_supply_demand_MAP/analysis/eFAST/sobol/collect_global_convergence_exports.py --remote-root /share/lab_crd/taoli/Project/soft_couping_eFAST/oxygen/results/eFAST/sobol/global_convergence_20261021 --local-root oxygen/results/eFAST/sobol/global_convergence_20261021
```

The exact collector executed for this download is archived as
`../collection_20261029T194045Z/provenance/collector_used.py`; its SHA-256 matches the collection
receipt. The current collector also supports filling missing exports in a
partially populated result directory while preserving existing identical files.

## Scientific definitions

- Two output-specific analyses; all 14 ST factors and 36 distinct S2 pairs that
  involve p_mis_base, mu_hp, or gamma_mu.
- Five independent scramble repeats, 201 oxygen points, unchanged global bounds.
- Repeat range <=0.05 and adjacent-resolution mean absolute change <=0.02.
- N65536 uses N32768 as its previous resolution; N131072 uses N65536.
- Repeat min-max bands are not confidence intervals; S2 sign is not a direction
  or compensation measure. Negative estimates are preserved as sampling noise.
- The original exported arrays and TSVs are checked independently before plotting.
- Oxygen keys are rounded to 12 decimals only for joins; index values and
  statistical calculations are not rounded. Original files stay byte-identical.
- Whole-curve ranking correlations include near-zero pairs; interpret them with
  absolute effect sizes. Top-six display selection uses maximum mean S2 across
  oxygen for each output, with all 36 pairs retained in the appendix and sources.

See `interpretation_zh.md`, `validation_data.json`, `validation_assets.json`,
`visual_review.md`, and the source/manifest files for evidence and coverage.
