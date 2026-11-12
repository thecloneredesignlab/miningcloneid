# Completed global Sobol / eFAST comparison

**Primary result:** existing Sobol N131072. **FAST batch:** N2049,4097,8193,16385, M4, five phases, 201 oxygen points. All computations completed; FAST ploidy ST remains unstable. Convergence is diagnostic-only. Neighborhoods remain archived.

- [Short interpretation for Noemi](summary_for_noemi_en.md)
- [Chinese interpretation](interpretation_zh.md)
- [English Methods / Results with robustness findings](publication/methods_results_en.md)
- [Final figure package](publication/README.md)
- [Sobol S1/ST/focal-S2 summary](publication/sobol/sensitivity_summary.pdf)
- [Highest-N eFAST S1/ST summary](publication/efast/N16385/sensitivity_summary.pdf)
- [Matched ST heatmaps](publication/comparison/method_comparison_st.pdf)
- [Structural zero-effect check](publication/comparison/structural_null_O2_zero.pdf)
- [Conditional response surfaces](response_surfaces/README_zh.md)
- [Method alignment](method_alignment.tsv)
- [Independent numerical verification](validation/independent_exports.json)
- [Final nonoverlapping-band diagnostic counts](comparison/diagnostic_summary_disjoint_bands.tsv)
- [Original transfer manifest](collection_manifest.json)

Source tables and five-repeat indices are in sobol/, efast/N*/, comparison/ and inputs/. Original HPC exports remain unchanged. Publication curves correct title/oxygen-label and legend spacing without rerunning models or estimators. Runtime, source commits, scheduling and recovery receipts are retained.
