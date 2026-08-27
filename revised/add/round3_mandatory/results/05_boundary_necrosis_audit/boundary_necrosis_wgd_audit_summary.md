# Boundary loss, necrosis, and WGD audit

Across the complete functional-response grids for all six frozen joint winners, the maximum boundary contribution is 26.75% of the combined boundary plus within-grid CIN-associated nonviable loss, and as many as 5.15% of rows in a context exceed the 1% materiality threshold.
The published seed-25 `necrosis_fit.tsv` contains six observed values but six missing predictions. Reconstructing predictions as terminal dead volume / terminal total volume yields six finite predictions and reproduces the reported normalized necrosis objective to maximum absolute error 1.82e-14 across the standalone fit and six joint winners.

WGD is implemented as a per-division competing branch: source state N maps to 2N with branch weight +1, mixed with the non-WGD branch by p_wgd, and out-of-grid transitions are dropped under the default boundary mode. The exact files and line ranges are recorded in `model_implementation_audit.tsv`.

Interpretation: boundary truncation is material for a subset of in-vivo functional-grid evaluations and therefore requires the separate expanded-grid sensitivity audit. The missing necrosis export was also a real audit defect; the reconstructed table repairs auditability without refitting and does not change the objective.
