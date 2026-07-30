# miningcloneid manuscript and result review — joint-preprocessing revision

## Overall conclusion

The manuscript has a strong mechanistic premise and several directionally stable results, but it requires **major revision** before external submission or public release.

Incorporating `joint_pre.zip` materially improves the audit because the full-500 fixed-O2 curves, response classes, spectral-gap diagnostics, landscape clusters and warm-start selection can now be checked. It also exposes stronger limitations than the top-10 package alone:

- only 141/500 fixed-O2 curves are spectrally reliable;
- none of the six joint warm-start representatives is reliable;
- all six joint pairings use the same in-vitro seed;
- the Figure 6 objective overlay uses raw burden+ploidy likelihood rather than the fitted MAP objective;
- selected joint solutions remain warm-start dependent and all local refinements were rejected.

The most important scientific correction remains the altered-daughter survival function: **absolute fitted survival is higher in vivo, but the ploidy-dependent survival gradient is weaker in vivo; culture has the stronger ploidy dependence.**

## Main deliverables

- `01_core_results_evidence_index.md`: result-by-result evidence map, exact locations, supported conclusions, scientific significance and limitations.
- `02_expansion_and_correction_notes.md`: every substantive addition/correction made to complete the manuscript.
- `03_revised_manuscript.tex`: completed English LaTeX manuscript.
- `04_revision_diff.md`: section-by-section semantic difference record relative to Repo manuscript.
- `05_peer_review_report.docx`: formal copy/paste-ready review report; no tables.
- `06_major_comment_crosswalk.docx`: detailed evidence trace for each major comment; no tables.
- `07_analysis_validation_report.md`: integrated methodology, calculation, figure and reproducibility QA.

## Reanalysis material

- `analysis/`: primary v2 portable audit tables and data dictionary.
- `analysis_tables/`: additional detailed audit tables generated during independent checking.
- `scripts/reanalyse_results.py`: reproducible Python reanalysis script.
- `MANIFEST.sha256`: SHA-256 checksums for all files in the delivery.

## Highest-priority corrections before submission

1. Regenerate Figure 5 and all text with the corrected survival-function direction.
2. Describe all 14 parameters as soft-coupled with the actual Welsch penalty.
3. Expand joint starts across the in-vitro landscape and diagnose rejected local refinements.
4. Quantify practical identifiability and bound sensitivity of headline parameters and derived functions.
5. Display fixed-O2 spectral-gap reliability and avoid treating small-gap classes as stable mechanisms.
6. Recompute Figure 6 ranking with the MAP objective or explicitly relabel it as raw likelihood.
7. Repair `necrosis_fit.tsv`, align the loss equation to implementation, and add an objective-reconstruction test.
8. Release Figure 4 seed-level AUC inputs/results.
9. Separate biological CIN nonviability from finite-grid boundary loss and run grid sensitivity.
10. Rebuild the final figures from a clean tagged commit with resolved configuration and full provenance.

## LaTeX QA

`03_revised_manuscript.tex` passed:

- balanced-brace and environment checks;
- duplicate/missing label and reference checks;
- placeholder scan for `XX`, `TODO`, `TBD`, and unresolved editorial notes;
- two-pass `pdflatex` syntax compilation with placeholder copies of the Repo figures/tables.

The production build still requires the Repo's actual figures, parameter tables and bibliography. Figure 5 and Figure 6 should be regenerated or relabeled before that final build.

## Word-document QA

Both Word files were rendered with the canonical DOCX renderer and inspected page by page. They contain plain paragraphs and headings only, with no tables or advanced structures, so they can be copied directly into a journal review system.
