#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_path <- normalizePath(script_arg[[1]], mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(script_path), "..", "..", "..", ".."), mustWork = TRUE)

add_root <- file.path(repo_root, "revised", "add", "round3_mandatory")
zip_path <- file.path(
  repo_root,
  "revised", "SuppFiles", "round3", "LTEE_revision_review_package",
  "03_revised_manuscript_source.zip"
)
manuscript_dir <- file.path(repo_root, "revised", "iteration1", "manuscript")
output_tex <- file.path(manuscript_dir, "ltee_hypoxia_model_round3_integrated.tex")
corrected_auc_table <- file.path(
  add_root, "tables", "supp_invivo_continuous_o2_auc_round3_corrected.tex"
)
provenance_path <- file.path(add_root, "provenance", "integrated_manuscript_build.tsv")

dir.create(dirname(corrected_auc_table), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(provenance_path), recursive = TRUE, showWarnings = FALSE)

read_zip_entry <- function(zip_file, entry) {
  lines <- system2(
    "unzip",
    c("-p", zip_file, entry),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(lines, "status")
  if (!is.null(status) && status != 0L) {
    stop("Unable to read zip entry ", entry, "; unzip status=", status, call. = FALSE)
  }
  paste(lines, collapse = "\n")
}

match_count_fixed <- function(text, pattern) {
  hit <- gregexpr(pattern, text, fixed = TRUE)[[1]]
  if (identical(hit[[1]], -1L)) 0L else length(hit)
}

replace_once <- function(text, old, new, label) {
  n <- match_count_fixed(text, old)
  if (n != 1L) {
    stop(label, ": expected exactly one match, found ", n, call. = FALSE)
  }
  sub(old, new, text, fixed = TRUE)
}

insert_after <- function(text, anchor, addition, label) {
  replace_once(text, anchor, paste0(anchor, "\n\n", addition), label)
}

replace_between_markers <- function(text, start_marker, end_marker, replacement, label) {
  start_at <- regexpr(start_marker, text, fixed = TRUE)[[1]]
  end_at <- regexpr(end_marker, text, fixed = TRUE)[[1]]
  if (start_at < 1L || end_at < 1L || end_at <= start_at) {
    stop(label, ": markers not found in the expected order", call. = FALSE)
  }
  start_end <- start_at + nchar(start_marker) - 1L
  paste0(
    substr(text, 1L, start_end),
    "\n\n", replacement, "\n\n",
    substr(text, end_at, nchar(text))
  )
}

source_entry <- "manuscript/ltee_hypoxia_model.tex"
source_text <- read_zip_entry(zip_path, source_entry)
auc_text <- read_zip_entry(
  zip_path, "manuscript/tables/supp_invivo_continuous_o2_auc.tex"
)
writeLines(auc_text, corrected_auc_table, useBytes = TRUE)

source_text <- replace_once(
  source_text,
  "\\documentclass[11pt]{article}",
  paste(
    "\\documentclass[11pt]{article}",
    "% Round-3 integrated evidence-bound manuscript.",
    "% Generated deterministically by revised/add/round3_mandatory/scripts/06_build_integrated_manuscript.R.",
    "% Text baseline: 03_revised_manuscript_source.zip; original iteration1 manuscript is not modified.",
    sep = "\n"
  ),
  "provenance header"
)

new_abstract <- r"{Chromosome instability (CIN) generates karyotypic diversity, but whether altered progeny persist may depend on resource availability and ploidy. We combined lineage-matched near-diploid and near-tetraploid SUM159 culture and tumor datasets with a chromosome-number-structured mechanistic model of proliferation, death, whole-genome doubling, chromosome missegregation, and post-missegregation survival. Oxygen deprivation altered culture growth and low- versus high-chromosome components. Across the ten lowest-objective culture fits, direct hypoxia-origin death in the 4N lineage remained below 1.9\%, supporting redistribution through chromosome-number transitions within this model. Quantitative endpoint checking nevertheless showed that the shared terminal distribution assigned 0.325 probability to the high-chromosome state and was incompatible with the O2 lineage composition (16/20 high-state cells; exact predictive \(p=1.75\times10^{-5}\)), although O1 was not discrepant by this check. Six seed-10-anchored joint-fit winners yielded concordant directions for 11 of 14 scalar context contrasts and for the evaluated missegregation and survival functions, but these comparisons were in-sample, shared one culture anchor, and used a substantially saturated bounded-Welsch penalty. Across all tumor records, observed burden lay inside the six-winner numerical range at only 14 of 92 time points, and terminal mean chromosome number did so for four of eight tumors. Fixed-O$_2$ dominant-eigenvector summaries had near-zero overall median day-1000 differences, but small-spectral-gap cases retained median differences of 0.331 and 0.601 from 2N and 4N starts. Expanding only the chromosome-grid upper bound changed dominant mean ploidy by as much as 10.543 ploidy units and drove low-O$_2$ attractors toward the expanded boundary. The supported contribution is therefore a conditional, hypothesis-generating comparison of fitted response functions together with explicit predictive, identifiability, finite-time, and boundary failure modes; it is not a causal or grid-independent resolution of the oxygen--CIN--ploidy landscape.}"
source_text <- replace_between_markers(
  source_text,
  "\\section*{Abstract}",
  "\\section*{Introduction}",
  new_abstract,
  "abstract"
)

methods_anchor <- r"{We distinguished three evidentiary levels throughout the analysis. Experimental burden, chromosome-number, flow-density, and necrosis records are direct observations. Fitted trajectories, latent \textit{in vivo} O$_2$, and context-specific response functions are conditional model outputs constrained by those observations. Fixed-O$_2$ dominant-eigenvector curves are post-fit asymptotic diagnostics under standardized conditions and are therefore more extrapolative than fitted finite-time trajectories. Numerical optimizer seeds and warm-start pairs assess search dependence; they do not create additional biological replication or sampling-based confidence intervals.}"
methods_addition <- r"{\paragraph{Round-3 endpoint and full-sample adequacy audits.}
For the two parallel oxygen-deprived 2N lineages, we compared each 20-cell terminal empirical distribution with the same fitted passage-34 chromosome-number distribution. Prespecified diagnostics were total-variation distance, chromosome-number Wasserstein-1 distance, Jensen--Shannon distance, and a two-sided exact binomial predictive check for the number of cells with \(N\geq80\). This predictive check tests compatibility with the fitted high-state probability and is not a test of biological replication. Full-sample in-sample checks used all 114 finite passage-growth records, 12 karyotyped passages (220 cells), 20 flow-density curves, 92 nonzero tumor-burden time points, and eight terminal tumor chromosome-number profiles. For joint fits, min--max ranges across the six frozen winners quantify numerical solution multiplicity only; they are neither confidence nor posterior intervals.

\paragraph{Finite-time and chromosome-grid sensitivity.}
We validated the separate-fit fixed-O$_2$ dominant-eigenvector summaries against matrix-exponential trajectories for 500 parameter sets, 201 oxygen values, 2N and 4N initial states, and 13 time horizons through day 1000. Differences were stratified by the prespecified spectral-gap regions \(<0.005\) and \(\geq0.01\). To isolate upper-grid truncation without redefining biological viability, we then recomputed the six joint-winner operators on the reported \(N=22\)--154 grid and on an upper-expanded \(N=22\)--308 grid, holding parameters fixed and performing no refitting. The lower boundary remained at \(N=22\); transitions below it therefore remained classified as loss from the modeled viable state space.

\paragraph{Necrosis, boundary, and WGD implementation audit.}
Boundary-routed mass was separated from within-grid nonviable daughter production over all exported oxygen and reference-ploidy rows for the six joint winners. Missing row-level necrosis predictions were reconstructed as terminal dead volume divided by terminal total volume and checked against the reported logit-scale objective with standard deviation 0.75. Source-code inspection verified that WGD is implemented as a constant per-division branch from \(N\) to \(2N\), with one offspring unit on the WGD branch, mixed with the non-WGD division branch before subtraction of the source division operator.

\paragraph{Anchor and regularization sensitivity scope.}
We audited the selected seed-10 culture anchor against the same-cluster near-optimal seed 132 and the best retained anchor from the alternative culture cluster, seed 157. We also quantified directional agreement and Welsch-penalty saturation among the six existing joint winners. Alternative-anchor and less-saturating-penalty refits require new optimization runs; they were specified as reproducible task matrices but were not represented as completed evidence in this manuscript.}"
source_text <- insert_after(
  source_text, methods_anchor, methods_addition, "round-3 methods"
)

o1_anchor <- r"{The supported result is therefore a shared process capable of generating and expanding high-chromosome states, rather than deterministic prediction of each lineage-specific endpoint.}"
o1_addition <- r"{The quantitative endpoint audit sharpened this limitation (Supplementary Table~\ref{tab:round3_o1_o2_adequacy}; Supplementary Figure~\ref{fig:round3_o1_o2_endpoint_adequacy}). The fitted high-state probability was 0.325. O1 contained eight of 20 high-state cells (observed fraction 0.400; exact predictive \(p=0.480\)), whereas O2 contained 16 of 20 (0.800; \(p=1.75\times10^{-5}\)). O1 and O2 total-variation distances from the shared fitted distribution were 0.844 and 0.894, and their Wasserstein-1 distances were 8.16 and 24.26 chromosomes, respectively. Thus, access to both modes is supported within the fitted process, but quantitative recovery of the lineage-specific endpoint distributions is not.}"
source_text <- insert_after(
  source_text, o1_anchor, o1_addition, "O1/O2 result"
)

invivo_anchor <- r"{Thus, the representative fit captured the direction and cohort-level endpoint of ploidy remodeling together with tumor growth, while its broader predicted terminal distributions and model-derived intermediate states preclude interpretation as precise mouse-level trajectory recovery.}"
invivo_addition <- r"{An all-record audit of the six frozen joint winners further constrained the strength of this fit claim (Supplementary Table~\ref{tab:round3_full_sample_adequacy}; Supplementary Figure~\ref{fig:round3_full_sample_adequacy}). Across 92 nonzero burden observations from eight tumors, the observed value lay inside the min--max range of the six numerical winners at 14 time points (15.2\%). Terminal observed mean chromosome number lay inside the corresponding six-winner range for four of eight tumors. These ranges are not uncertainty intervals, and the audit is not held-out validation; it shows that the fitted family can reproduce broad scales while leaving substantial record-level mismatch.}"
source_text <- insert_after(
  source_text, invivo_anchor, invivo_addition, "full-sample result"
)

joint_anchor <- r"{The observed directional concordance is therefore limited to the six sampled warm-start-conditioned search ensembles; it does not resolve the global optimum or behavior in unsampled regions of the objective landscape and does not identify unique context effect sizes.}"
joint_addition <- r"{The explicit anchor and penalty audit reinforced that boundary (Supplementary Table~\ref{tab:round3_anchor_joint_sensitivity}; Supplementary Figure~\ref{fig:round3_anchor_sensitivity}). Eleven of 14 scalar parameter contrasts had the same direction in all six current winners, but all six winners inherited seed 10 as their culture anchor. Seed 132 occupied the same culture-fit cluster and differed in objective by only 0.000722, whereas seed 157 was the best retained fit in the alternative cluster and had an objective 0.419 higher than seed 10. Among active coupling terms, the median fraction in the saturated Welsch region was 100\%. Without the prespecified alternative-anchor and less-saturating-penalty refits, the current analysis supports conditional directional concordance, not anchor- or regularization-invariant context effects.}"
source_text <- insert_after(
  source_text, joint_anchor, joint_addition, "anchor sensitivity result"
)

fixed_anchor <- r"{Because these curves are asymptotic operator diagnostics rather than fitted finite-time observations, small spectral gaps weaken their biological interpretation over experimental time horizons.}"
fixed_addition <- r"{Direct finite-time validation confirmed that this limitation is material (Supplementary Table~\ref{tab:round3_fixed_o2_finite_time}; Supplementary Figure~\ref{fig:round3_fixed_o2_finite_time}). By day 1000, the overall median absolute differences from the dominant-eigenvector mean were \(1.30\times10^{-7}\) and \(2.87\times10^{-7}\) from 2N and 4N starts, but the corresponding medians remained 0.331 and 0.601 in the spectral-gap-\(<0.005\) stratum. Curve-class agreement with the asymptotic classification was only 72.2\% and 70.6\%, although the two finite-time initial states agreed with one another for 96.2\% of seed--oxygen combinations.

Upper-grid sensitivity exposed a separate failure mode (Supplementary Tables~\ref{tab:round3_boundary_necrosis_wgd} and~\ref{tab:round3_expanded_grid}; Supplementary Figures~\ref{fig:round3_boundary_necrosis} and~\ref{fig:round3_expanded_grid}). Expanding only the upper chromosome bound from 154 to 308 changed a joint winner's dominant mean ploidy by as much as 10.543 ploidy units. At 0--0.5\% O$_2$, the median dominant mean shifted from approximately 6.4--6.8 on the reported grid to approximately 14 on the expanded grid, close to the new upper boundary. Consequently, these low-O$_2$ dominant modes cannot be interpreted as grid-independent equilibrium ploidies; they are retained as truncated-grid diagnostics of upward propensity and model instability.}"
source_text <- replace_once(
  source_text,
  fixed_anchor,
  paste0(fixed_addition, "\n\n", fixed_anchor),
  "finite-time and grid result"
)

source_text <- replace_once(
  source_text,
  r"{\caption{\textbf{Fixed-O$_2$ diagnostics retain competing asymptotic oxygen--CIN--ploidy responses.} Panel-by-panel legend continues on the following page.}}",
  r"{\caption{\textbf{Fixed-O$_2$ truncated-grid diagnostics expose competing and boundary-dependent oxygen--CIN--ploidy responses.} Panel-by-panel legend continues on the following page.}}",
  "Figure 6 title"
)
source_text <- replace_once(
  source_text,
  r"{These pooled class counts include reliable, caution, and unreliable spectral-gap strata.}",
  r"{These pooled class counts include reliable, caution, and unreliable spectral-gap strata. Finite-time and upper-grid audits show that low-gap curve classes and low-O$_2$ dominant modes cannot be read as universal finite-time trajectories or grid-independent biological equilibria.}",
  "Figure 6 panel A limitation"
)

source_text <- replace_once(
  source_text,
  r"{The main result is not a uniquely identified parameter vector, but a small set of function-level contrasts that remained directional across six warm-start-conditioned search ensembles.}",
  r"{The primary positive result is a set of function-level contrasts that remained directional within six seed-10-anchored, bounded-Welsch search ensembles, rather than an anchor-invariant effect or uniquely identified parameter vector.}",
  "Discussion primary result"
)
source_text <- replace_once(
  source_text,
  r"{This separation sharpens the biological significance of the model.}",
  r"{Within that conditional specification, this separation sharpens the biological hypothesis generated by the model.}",
  "Discussion hypothesis framing"
)

discussion_o1_anchor <- r"{The supported claim is therefore a common process capable of generating and filtering a mixed distribution, rather than complete reconstruction of lineage-specific evolution.}"
discussion_o1_addition <- r"{The broader adequacy audit leads to the same distinction at the tumor level. Growth RMSE and endpoint averages can coexist with extensive timepoint-level and distributional residuals: the six-winner numerical range contained only 14 of 92 nonzero burden observations and four of eight terminal mean chromosome numbers. Because these are in-sample fits and numerical ranges rather than predictive intervals, the present data support conditional reconstruction of selected aggregate features, not general predictive validation.}"
source_text <- insert_after(
  source_text, discussion_o1_anchor, discussion_o1_addition, "Discussion adequacy"
)

discussion_ident_anchor <- r"{The directional concordance reported here is limited to these six warm-start-conditioned search ensembles; a unique effect size, global optimum, and behavior in unsampled objective-landscape regions remain unresolved.}"
discussion_ident_addition <- r"{The new checks further show that asymptotic response topology is not fully robust to either time horizon or state-space truncation. Median day-1000 convergence was excellent overall, but remained poor for small spectral gaps, and the upper-expanded chromosome grid produced low-O$_2$ dominant modes near its new boundary. Figure~\ref{fig:iteration1-o2-linked-model-selection} should therefore be read as a map of competing truncated-grid model responses and failure modes, not as evidence for biological equilibrium ploidy classes.}"
source_text <- insert_after(
  source_text,
  discussion_ident_anchor,
  discussion_ident_addition,
  "Discussion finite-time and grid"
)

source_text <- replace_once(
  source_text,
  r"{Two implementation limitations require separate caution. First, out-of-range chromosome-transition mass is routed to the CIN-associated dead compartment, so that compartment combines biological nonviability with finite-grid boundary loss. Second, the necrosis likelihood is active in the objective, but the current per-sample \texttt{necrosis\_fit.tsv} exports contain missing predicted values; the terminal necrosis-fraction audit was therefore reconstructed from the stored model states and objective components. The export should be repaired before claiming fully automated row-level reproduction of the necrosis fit.}",
  r"{Two implementation audits require separate caution. First, out-of-range chromosome-transition mass is routed to the CIN-associated dead compartment, so that compartment combines within-grid biological nonviability with boundary-routed loss. Across the six winners, the maximum boundary contribution was 26.7\% of combined boundary plus within-grid CIN-associated loss, although upper-grid expansion removed at most 5.81\% of that current-grid quantity because much of the retained boundary term was below the biologically defined lower viable state. Second, the necrosis likelihood is active in the objective, but six predicted rows are missing from the seed-25 \texttt{necrosis\_fit.tsv} export. Reconstruction from terminal dead and total volumes supplied all six predictions and reproduced the reported objective exactly; this repairs the present audit trail but not the upstream exporter. WGD source-code inspection additionally confirmed a constant per-division \(N\rightarrow2N\) branch with one offspring unit, rather than a context-dependent WGD rate or two-offspring WGD branch.}",
  "Discussion implementation audit"
)

source_text <- replace_once(
  source_text,
  "\\input{tables/supp_invivo_continuous_o2_auc.tex}",
  "\\input{../../add/round3_mandatory/tables/supp_invivo_continuous_o2_auc_round3_corrected.tex}",
  "corrected AUC table"
)

source_text <- insert_after(
  source_text,
  "\\input{tables/supp_invivo_top10_fit_quality.tex}",
  paste(
    "\\input{../../add/round3_mandatory/tables/round3_full_sample_predictive_adequacy.tex}",
    "\\input{../../add/round3_mandatory/tables/round3_fixed_o2_finite_time_validation.tex}",
    "\\input{../../add/round3_mandatory/tables/round3_boundary_necrosis_wgd_audit.tex}",
    "\\input{../../add/round3_mandatory/tables/round3_expanded_grid_boundary_sensitivity.tex}",
    sep = "\n\n"
  ),
  "in-vivo audit tables"
)

source_text <- insert_after(
  source_text,
  "\\input{tables/supp_invitro_o1_o2_final_karyotype.tex}",
  "\\input{../../add/round3_mandatory/tables/round3_o1_o2_endpoint_adequacy.tex}",
  "O1/O2 audit table"
)

source_text <- insert_after(
  source_text,
  "\\input{tables/supp_invitro_top10_4n_trajectory.tex}",
  "\\input{../../add/round3_mandatory/tables/round3_anchor_joint_sensitivity.tex}",
  "anchor audit table"
)

invivo_figure_anchor <- paste(
  "\\label{fig:supp_all18_cluster_prior_violins}",
  "\\end{figure}",
  sep = "\n"
)
invivo_audit_figures <- r"{\clearpage
\begin{figure}[p]
\centering
\includegraphics[width=\textwidth,height=0.80\textheight,keepaspectratio]{../../add/round3_mandatory/figures/round3_full_sample_predictive_adequacy.png}
\caption{\textbf{All-record tumor in-sample adequacy audit.} The left panel compares every nonzero tumor-burden observation with predictions from the six frozen joint winners; the right panel compares all eight terminal mean chromosome-number observations. Winner ranges represent numerical solution multiplicity, not confidence or posterior intervals, and the figure is not held-out validation.}
\label{fig:round3_full_sample_adequacy}
\end{figure}

\clearpage
\begin{figure}[p]
\centering
\includegraphics[width=\textwidth,height=0.80\textheight,keepaspectratio]{../../add/round3_mandatory/figures/round3_fixed_o2_finite_time_validation.png}
\caption{\textbf{Finite-time validation of the fixed-O$_2$ dominant-eigenvector summaries.} Matrix-exponential trajectories from 2N and 4N initial states were compared with the dominant eigenvector through day 1000. Small spectral gaps retain material deviations and lower curve-class agreement even when the overall median error is near zero.}
\label{fig:round3_fixed_o2_finite_time}
\end{figure}

\clearpage
\begin{figure}[p]
\centering
\includegraphics[width=\textwidth,height=0.80\textheight,keepaspectratio]{../../add/round3_mandatory/figures/round3_boundary_necrosis_audit.png}
\caption{\textbf{Boundary-loss and necrosis-objective audit.} Boundary-routed and within-grid nonviable components are separated for the six frozen joint winners. Missing row-level necrosis predictions are reconstructed from terminal state volumes and reproduce the fitted objective.}
\label{fig:round3_boundary_necrosis}
\end{figure}

\clearpage
\begin{figure}[p]
\centering
\includegraphics[width=\textwidth,height=0.80\textheight,keepaspectratio]{../../add/round3_mandatory/figures/round3_expanded_grid_boundary_sensitivity.png}
\caption{\textbf{Upper-grid sensitivity without refitting.} The biological lower bound remains \(N=22\), while the upper bound is expanded from 154 to 308. Low-O$_2$ dominant modes move toward the expanded upper boundary, demonstrating that they are not grid-independent equilibrium ploidies.}
\label{fig:round3_expanded_grid}
\end{figure}}"
source_text <- insert_after(
  source_text,
  invivo_figure_anchor,
  invivo_audit_figures,
  "in-vivo audit figures"
)

o1_figure_anchor <- "\\input{../../add/round3_mandatory/tables/round3_anchor_joint_sensitivity.tex}"
o1_anchor_figures <- r"{\clearpage
\begin{figure}[p]
\centering
\includegraphics[width=\textwidth,height=0.80\textheight,keepaspectratio]{../../add/round3_mandatory/figures/round3_o1_o2_endpoint_adequacy.png}
\caption{\textbf{Quantitative shared-distribution check for O1 and O2.} Each lineage's 20-cell endpoint is compared with the same fitted passage-34 distribution. Exact binomial values are predictive model checks; the lineages are not asserted to be biological replicates.}
\label{fig:round3_o1_o2_endpoint_adequacy}
\end{figure}

\clearpage
\begin{figure}[p]
\centering
\includegraphics[width=\textwidth,height=0.80\textheight,keepaspectratio]{../../add/round3_mandatory/figures/round3_invitro_anchor_sensitivity.png}
\caption{\textbf{Culture-anchor and current joint-winner sensitivity audit.} Seed 132 is a same-cluster near-optimal culture fit; seed 157 is the best retained alternative-cluster fit. Existing joint-winner directionality and Welsch saturation are descriptive of the seed-10-anchored analysis and do not replace alternative-anchor or alternative-penalty refits.}
\label{fig:round3_anchor_sensitivity}
\end{figure}
\clearpage}"
source_text <- insert_after(
  source_text,
  o1_figure_anchor,
  o1_anchor_figures,
  "O1/O2 and anchor figures"
)

writeLines(source_text, output_tex, useBytes = TRUE)

hash_rows <- data.frame(
  artifact = c(
    "round3_source_zip",
    "corrected_auc_table",
    "integrated_manuscript"
  ),
  relative_path = c(
    file.path(
      "revised", "SuppFiles", "round3", "LTEE_revision_review_package",
      "03_revised_manuscript_source.zip"
    ),
    file.path(
      "revised", "add", "round3_mandatory", "tables",
      "supp_invivo_continuous_o2_auc_round3_corrected.tex"
    ),
    file.path(
      "revised", "iteration1", "manuscript",
      "ltee_hypoxia_model_round3_integrated.tex"
    )
  ),
  md5 = unname(tools::md5sum(c(zip_path, corrected_auc_table, output_tex))),
  bytes = unname(file.info(c(zip_path, corrected_auc_table, output_tex))$size),
  stringsAsFactors = FALSE
)
write.table(
  hash_rows,
  provenance_path,
  sep = "\t", row.names = FALSE, quote = FALSE
)

message("Wrote integrated manuscript: ", output_tex)
message("Wrote corrected AUC table: ", corrected_auc_table)
