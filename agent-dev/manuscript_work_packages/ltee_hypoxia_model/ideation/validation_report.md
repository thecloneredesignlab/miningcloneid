# Revised Ideation Package Validation

Date: 2026-07-24

Status: **PASS — Gate A approved and closed.** The ideation package is complete,
internally consistent, self-contained, and approved as the drafting contract for
Figures 1–5.

## Reproducible Build And Structural Checks

Commands:

```bash
Rscript agent-dev/manuscript_work_packages/ltee_hypoxia_model/ideation/build_ideation_report.R
Rscript agent-dev/manuscript_work_packages/ltee_hypoxia_model/ideation/validate_ideation_package.R
```

Results:

- 11 required ideation artifacts present;
- 16 unique historical-panel dispositions;
- 28 unique TeX panel-role records;
- 24 unique source-inventory records;
- 18 unique recorded feedback decisions;
- 3 unique fit-quality evidence records;
- 5 embedded historical PNGs;
- all five CSVs parse without missing fields;
- all six required report sections appear in the required order;
- report references are self-contained data URIs or internal anchors;
- no external stylesheet, script, or HTTP asset reference is present; and
- the current and earlier feedback-manager spans are both linked in the handoff
  context.

Final report SHA-256:

```text
585b76bf002970a1e52f7d771a7c4814851eb58b709dc01a0540a3594e203e6a
```

## Scientific And Source-Basis Checks

- Figure 3 controls were verified as the normoxic left-hand branches for both
  2N and 4N lineages. They were included in fitting, are not held-out
  validation, end earlier than the deprived branches, and are plotted by lineage
  passage rather than literal calendar time.
- Figure 3 seed 10 is rank 1/500, but 12 starts lie within 0.01 objective units;
  the local optimizer ended with code 1, five active parameters are at bounds,
  and the stored DEoptim convergence flags are 0/500.
- Historical Figure 4 panels A-E were verified as in-vivo analyses. Historical
  panel F deliberately uses a pooled in-vivo/in-vitro embedding. Following the
  user's correction, its point universe, saved coordinates, context markers,
  overlays, and geometry are preserved unchanged; no context-specific
  re-embedding is proposed.
- Formal in-vivo solution regions were verified as `vi_C01` (99/500), `vi_C02`
  (385/500), and `vi_C03` (16/500). The user confirmed that all three are
  already displayed in the pooled embedding and should remain; no two-regime
  mapping is required.
- The low-O2 feature was verified as `mu_hp` (AUC approximately 0.849), and the
  high-O2 feature as `p_mis_base` (AUC approximately 0.903). These are labeled as
  one-feature discrimination, not causal influence.
- Figure 4 seed 25 was verified as the best total weighted-MAP objective among
  500 starts, with 29 fits within 1% and 241 within 5%; its component ranks are
  52nd for ploidy and 189th for burden.
- All six approved Figure 5 pair winners and their objective arithmetic were
  verified. The stored checks are selection/provenance validation, not held-out
  predictive validation; no L-BFGS-B refinement was accepted.
- Separate and joint necrosis exports contain only `NA` predicted values. The
  revised fit-quality plans exclude necrosis pending resolution.

## Render Review

The rebuilt self-contained report was served locally and inspected with a
browser at 1200-pixel desktop width and 390-pixel mobile width.

- title, sticky navigation, metrics, section order, callouts, and embedded
  historical figures rendered;
- the desktop document remained within the viewport except for intentionally
  scrollable wide tables;
- the mobile layout stacked figure cards and figure plans, with tables contained
  by horizontal scroll wrappers; and
- the browser accessibility tree exposed all headings, navigation links, tables,
  lists, and decision text.

The corrected historical Figure 4 was also visually re-inspected at original
resolution. Panel F visibly retains both contexts in one embedding: in-vivo
points use circles and in-vitro points use triangles. The revised plan now
preserves that shared point universe and geometry.

## Feedback Archive Integrity

`shasum -a 256 -c SHA256SUMS` passed for:

- archive README;
- verbatim raw feedback;
- immutable reviewed HTML;
- archive manifest; and
- recorded Git context.

## Approved Drafting Boundary

- Gate A is closed and the five-figure architecture is approved for drafting.
- Compact direct observed–predicted blocks in Figures 3–5 are approved; objective
  distributions and optimizer diagnostics belong in supplementary material.
- All three formal in-vivo solution regions remain in the unchanged pooled
  embedding.
- Figures 3-5 can show fit adequacy and solution sensitivity, not formal
  statistical confidence.
- The exact result bundle is ignored/local and must be reduced to durable
  panel-ready tables plus hashes during drafting.
- The restored TeX still lacks referenced figure assets, two parameter-table
  TeX files, and the bibliography, so it is not yet a manuscript-ready compile.
- Figure 6, Figure 3E, optimizer refits, and new analytical grids remain
  explicitly outside the authorized scope.
- Production manuscript edits and resolution of the necrosis export remain
  outside the authorized scope.
