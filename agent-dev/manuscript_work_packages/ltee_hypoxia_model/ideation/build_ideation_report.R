#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
output_path <- if (length(args) >= 1L) {
  args[[1L]]
} else {
  "agent-dev/manuscript_work_packages/ltee_hypoxia_model/ideation/ideation_report.html"
}

repo_root <- Sys.getenv("MININGCLONEID_REPO_ROOT", unset = getwd())
package_root <- dirname(output_path)

input_path <- function(name) file.path(package_root, name)
required <- c(
  input_path("existing_panel_disposition.csv"),
  input_path("tex_figure_requirements.csv"),
  input_path("source_inventory.csv"),
  input_path("feedback_decisions.csv"),
  input_path("fit_quality_evidence.csv"),
  input_path("existing_context_manifest.tsv"),
  input_path("baseline_build.md"),
  input_path("feedback_manager_context.md"),
  input_path("ideas.md")
)

missing <- required[!file.exists(required)]
if (length(missing) > 0L) {
  stop("Missing report input(s): ", paste(missing, collapse = ", "))
}
if (!requireNamespace("base64enc", quietly = TRUE)) {
  stop("Package base64enc is required to build the self-contained report")
}

disposition <- read.csv(
  input_path("existing_panel_disposition.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
requirements <- read.csv(
  input_path("tex_figure_requirements.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
inventory <- read.csv(
  input_path("source_inventory.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
decisions <- read.csv(
  input_path("feedback_decisions.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
fit_quality <- read.csv(
  input_path("fit_quality_evidence.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
context_manifest <- read.delim(
  input_path("existing_context_manifest.tsv"),
  sep = "\t",
  quote = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

html_escape <- function(x) {
  x <- ifelse(is.na(x), "", as.character(x))
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub('"', "&quot;", x, fixed = TRUE)
  x
}

root_path <- function(path) {
  if (grepl("^/", path)) path else file.path(repo_root, path)
}

data_uri <- function(path) {
  resolved <- root_path(path)
  if (!file.exists(resolved)) return(NA_character_)
  paste0("data:image/png;base64,", base64enc::base64encode(resolved))
}

table_html <- function(data, columns, class_name = "") {
  columns <- intersect(columns, names(data))
  data <- data[, columns, drop = FALSE]
  header <- paste0("<th>", html_escape(columns), "</th>", collapse = "")
  if (nrow(data) == 0L) {
    body <- paste0(
      "<tr><td colspan='", length(columns),
      "' class='muted'>No records.</td></tr>"
    )
  } else {
    rows <- lapply(seq_len(nrow(data)), function(index) {
      values <- unlist(data[index, , drop = FALSE], use.names = FALSE)
      paste0(
        "<tr>",
        paste0("<td>", html_escape(values), "</td>", collapse = ""),
        "</tr>"
      )
    })
    body <- paste0(rows, collapse = "\n")
  }
  paste0(
    "<div class='table-wrap ", html_escape(class_name), "'><table><thead><tr>",
    header, "</tr></thead><tbody>", body, "</tbody></table></div>"
  )
}

badge <- function(text, class_name = "neutral") {
  paste0(
    "<span class='badge ", html_escape(class_name), "'>",
    html_escape(text), "</span>"
  )
}

context_record <- function(filename) {
  hits <- context_manifest[
    basename(context_manifest$local_review_copy) == filename,
    ,
    drop = FALSE
  ]
  if (nrow(hits) == 0L) return(NULL)
  hits[1L, , drop = FALSE]
}

figure_card <- function(filename, title, role, concern) {
  record <- context_record(filename)
  relative_path <- if (is.null(record)) {
    file.path(package_root, "existing_context", filename)
  } else {
    record$local_review_copy[[1L]]
  }
  uri <- data_uri(relative_path)
  image_html <- if (is.na(uri)) {
    paste0(
      "<div class='missing'>Historical image unavailable<br><code>",
      html_escape(relative_path), "</code></div>"
    )
  } else {
    paste0(
      "<img src='", uri, "' alt='", html_escape(title),
      "' loading='lazy'>"
    )
  }
  provenance <- if (is.null(record)) {
    ""
  } else {
    paste0(
      "<p class='micro'><strong>Historical source:</strong> commit <code>",
      html_escape(substr(record$source_commit[[1L]], 1L, 12L)),
      "</code>; SHA-256 <code>",
      html_escape(substr(record$sha256[[1L]], 1L, 16L)),
      "…</code>; ", html_escape(record$width_px[[1L]]), " × ",
      html_escape(record$height_px[[1L]]), " px.</p>"
    )
  }
  paste0(
    "<article class='figure-card'><div class='card-head'><h3>",
    html_escape(title), "</h3>", badge("historical context", "warn"),
    "</div>", image_html,
    "<p><strong>Useful role:</strong> ", html_escape(role), "</p>",
    "<p><strong>Main concern:</strong> ", html_escape(concern), "</p>",
    provenance,
    "<p class='micro'><code>", html_escape(relative_path), "</code></p></article>"
  )
}

status_counts <- as.data.frame(
  table(requirements$evidence_status),
  stringsAsFactors = FALSE
)
names(status_counts) <- c("Evidence status", "Panel roles")

disposition_counts <- as.data.frame(
  table(disposition$disposition),
  stringsAsFactors = FALSE
)
names(disposition_counts) <- c("Disposition", "Historical panels")

figure_counts <- aggregate(
  panel_id ~ sub("[A-Z]$", "", requirements$panel_id),
  data = requirements,
  FUN = length
)
names(figure_counts) <- c("Figure", "Proposed panel roles")

existing_cards <- paste0(
  figure_card(
    "assembled_fig1.png",
    "Historical Figure 1",
    "A source-table-driven matched-lineage timeline.",
    "It shows design and sampling, but not the observed chromosome trajectories used to motivate the Results claim."
  ),
  figure_card(
    "assembled_fig3.png",
    "Historical Figure 3",
    "Fit, predicted distribution, survival function, and mechanism diagnostics.",
    "Panel A is too composite. Its left-side control comparators are valid and should be retained, but the revised A should show growth rate over passage only. Figure 3E is rejected."
  ),
  figure_card(
    "assembled_fig4.png",
    "Historical Figure 4",
    "Latent effective-O2 trajectories, AUC screens, and solution-landscape context.",
    "Panels A-E are in-vivo analyses, while panel F intentionally embeds both in-vivo and in-vitro fits. The pooled embedding must remain unchanged; AUCs should still be described as separators rather than causal effects."
  ),
  figure_card(
    "assembled_fig5.png",
    "Historical Figure 5",
    "Context-specific fitted functions and parameter-ratio context.",
    "Panels B-D and E came from incompatible joint runs. The revised figure must show fit adequacy first, retain a full all-six ratio view second, and regenerate the three functions from the same six winners."
  ),
  figure_card(
    "assembled_fig6.png",
    "Historical Figure 6",
    "Fixed-O2 response-curve classes.",
    "The user has deferred Figure 6, its framing, and its missing analytical map. It remains outside the approved five-figure drafting scope."
  )
)

css <- "
:root {
  --ink:#1f2a32; --muted:#5c6972; --paper:#ffffff; --canvas:#eaf0f2;
  --navy:#163747; --blue:#26708f; --line:#ccd7dc; --soft:#f5f8f9;
  --green:#16785b; --green-soft:#e8f6f0; --amber:#ad6a00;
  --amber-soft:#fff5dc; --red:#a23b27; --red-soft:#fff0eb;
  --purple:#67528d; --purple-soft:#f3effa;
}
* { box-sizing:border-box; }
html { scroll-behavior:smooth; }
body {
  margin:0; color:var(--ink); background:var(--canvas);
  font-family:Inter,Arial,Helvetica,sans-serif; line-height:1.48;
}
header {
  padding:46px max(28px,calc((100vw - 1240px)/2));
  color:white; background:linear-gradient(125deg,#102c39,#1c5269);
  border-bottom:7px solid #70b7c9;
}
header h1 { margin:0 0 10px; max-width:1000px; font-size:36px; line-height:1.12; }
header p { margin:0; max-width:940px; color:#dcebf0; font-size:17px; }
main { max-width:1240px; margin:0 auto; padding:26px 38px 70px; background:var(--paper); }
nav {
  position:sticky; top:0; z-index:5; display:flex; flex-wrap:wrap; gap:6px 18px;
  margin:-26px -38px 28px; padding:14px 38px; background:rgba(255,255,255,.97);
  border-bottom:1px solid var(--line);
}
nav a { color:#22566d; text-decoration:none; font-size:13px; font-weight:700; }
section { scroll-margin-top:70px; }
h2 {
  margin:46px 0 16px; padding-bottom:8px; border-bottom:2px solid #b8c8cf;
  color:var(--navy); font-size:25px;
}
h3 { margin:24px 0 9px; color:#263b46; font-size:18px; }
h4 { margin:16px 0 6px; font-size:15px; }
p, li { font-size:15px; }
li { margin:6px 0; }
code { overflow-wrap:anywhere; font-size:.9em; }
.lede {
  padding:18px 20px; background:#edf7fa; border-left:5px solid var(--blue);
  font-size:18px; font-weight:700; color:#163f52;
}
.notice, .warning, .critical, .boundary {
  margin:16px 0; padding:14px 17px; border-left:4px solid;
}
.notice { border-color:var(--green); background:var(--green-soft); }
.warning { border-color:var(--amber); background:var(--amber-soft); }
.critical { border-color:var(--red); background:var(--red-soft); }
.boundary { border-color:var(--purple); background:var(--purple-soft); }
.summary-grid {
  display:grid; grid-template-columns:repeat(auto-fit,minmax(205px,1fr));
  gap:12px; margin:18px 0;
}
.metric { padding:14px 16px; border:1px solid var(--line); background:var(--soft); }
.metric strong { display:block; font-size:23px; color:var(--navy); }
.metric span { display:block; color:var(--muted); font-size:12px; }
.cards {
  display:grid; grid-template-columns:repeat(auto-fit,minmax(360px,1fr));
  gap:18px; align-items:start;
}
.figure-card { padding:15px; border:1px solid var(--line); background:#fff; }
.card-head { display:flex; justify-content:space-between; gap:8px; align-items:start; }
.card-head h3 { margin-top:0; }
.figure-card img {
  display:block; width:100%; height:310px; object-fit:contain;
  border:1px solid #e5ecef; background:#fff;
}
.figure-card p { margin:9px 0; font-size:13px; }
.micro { color:var(--muted); font-size:11px !important; }
.missing {
  min-height:220px; display:flex; align-items:center; justify-content:center;
  text-align:center; padding:20px; color:var(--muted); background:var(--soft);
  border:1px dashed #aab9c0;
}
.badge {
  display:inline-block; padding:3px 8px; border-radius:999px; white-space:nowrap;
  font-size:11px; font-weight:800; letter-spacing:.01em;
}
.badge.good { color:#0f6449; background:#dff2e9; }
.badge.warn { color:#754700; background:#ffedbc; }
.badge.stop { color:#7d2c1c; background:#f8dcd5; }
.badge.neutral { color:#3f515a; background:#e5ecef; }
.figure-plan {
  display:grid; grid-template-columns:84px 1fr; gap:0; margin:10px 0;
  border:1px solid var(--line);
}
.figure-plan .label {
  padding:14px 12px; color:white; background:#2d596d; font-weight:800;
}
.figure-plan .content { padding:12px 16px; background:#fbfcfc; }
.figure-plan .content p { margin:4px 0; }
.option {
  margin:17px 0; padding:18px 20px; border:1px solid var(--line);
  border-left:6px solid var(--purple); background:#fff;
}
.option.recommended { border-left-color:var(--green); background:#f3fbf7; }
.option h3 { margin-top:0; }
.table-wrap {
  width:100%; overflow:auto; margin:12px 0 22px; border:1px solid var(--line);
}
table { width:100%; min-width:820px; border-collapse:collapse; font-size:12px; }
th {
  padding:9px 10px; text-align:left; vertical-align:bottom;
  color:white; background:#344f5c; position:sticky; top:44px;
}
td { padding:9px 10px; vertical-align:top; border-bottom:1px solid #e2e8eb; }
tr:nth-child(even) td { background:#f7f9fa; }
.compact table { min-width:520px; }
details { margin:12px 0; padding:11px 14px; border:1px solid var(--line); }
summary { cursor:pointer; font-weight:800; color:#294d5e; }
.checklist { padding-left:26px; }
.checklist li { margin:10px 0; }
.muted { color:var(--muted); }
footer {
  margin-top:52px; padding-top:16px; color:var(--muted);
  border-top:1px solid var(--line); font-size:12px;
}
@media (max-width:720px) {
  header { padding:30px 20px; } header h1 { font-size:29px; }
  main { padding:20px; } nav { position:static; margin:-20px -20px 24px; padding:12px 20px; }
  .cards { grid-template-columns:1fr; }
  .figure-plan { grid-template-columns:1fr; }
  th { position:static; }
}
"

html <- paste0(
  "<!doctype html><html lang='en'><head><meta charset='utf-8'>",
  "<meta name='viewport' content='width=device-width,initial-scale=1'>",
  "<link rel='icon' href='data:,'>",
  "<title>Approved LTEE hypoxia manuscript figure plan</title>",
  "<style>", css, "</style></head><body>",
  "<header><h1>LTEE hypoxia manuscript: approved five-figure plan</h1>",
  "<p>Gate A closed on July 24, 2026. This surface records the approved architecture, resolved corrections, evidence boundaries, and drafting scope. It contains no draft figures and authorizes no work beyond the explicit drafting boundary.</p></header>",
  "<main><nav>",
  "<a href='#existing'>1 · Existing context</a>",
  "<a href='#feedback'>2 · Review decisions</a>",
  "<a href='#disposition'>3 · Disposition</a>",
  "<a href='#ideas'>4 · Working architecture</a>",
  "<a href='#decisions'>5 · Gate A resolution</a>",
  "<a href='#appendix'>6 · Evidence appendix</a>",
  "</nav>",

  "<p class='lede'>Gate A approved · Matched observations → model assumptions → in-vitro fit and mechanism → in-vivo fit and solution landscape → joint-fit context differences</p>",
  "<div class='summary-grid'>",
  "<div class='metric'><strong>5</strong><span>active main figures</span></div>",
  "<div class='metric'><strong>3</strong><span>fit-quality openings in Figures 3–5</span></div>",
  "<div class='metric'><strong>1</strong><span>figure explicitly deferred</span></div>",
  "<div class='metric'><strong>0</strong><span>refits authorized</span></div>",
  "</div>",

  "<section id='existing'><h2>1. Existing Figures And Panels</h2>",
  "<p>The restored TeX defines six historical main-figure concepts. Five historical composites are embedded below only to make the requested changes concrete; Figure 2 is caption-only. These images are review context, not draft assets.</p>",
  "<div class='warning'><strong>Historical-use boundary.</strong> Every accepted panel must later be regenerated from declared observations, exact fit tables or objects, and repository-local generator code. No historical raster pixels may enter a manuscript figure.</div>",
  "<div class='cards'>", existing_cards,
  "<article class='figure-card'><div class='card-head'><h3>Figure 2</h3>",
  badge("greenfield", "neutral"), "</div>",
  "<div class='missing'>No historical panel exists.<br>The TeX contains only a five-part model-overview caption.</div>",
  "<p><strong>Useful role:</strong> One integrated schematic connecting resource context, growth/death cost, CIN/WGD generation, post-missegregation survival, and chromosome-state outcomes.</p>",
  "<p><strong>Main concern:</strong> It must show implemented model assumptions and mechanisms, not imply independent experimental proof.</p></article>",
  "</div>",
  "<h3>Restored TeX scaffold mapped to the revised proposal</h3>",
  table_html(
    requirements,
    c("panel_id", "current_caption_role", "evidence_status", "proposed_role", "claim_guardrail")
  ),
  "</section>",

  "<section id='feedback'><h2>2. Relevant Feedback And Recorded Decisions</h2>",
  "<p>The review and subsequent corrections were archived and served through feedback-manager as <code>SPAN-20260724-001</code> through <code>SPAN-20260724-005</code>. The table below records accepted choices, rejected proposals, source corrections, and final Gate A approval; no drafting decision remains open.</p>",
  table_html(
    decisions,
    c("decision_id", "figure_or_topic", "review_state", "recorded_decision", "implementation_in_revised_ideation", "remaining_question")
  ),
  "<h3>Figure-by-figure interpretation</h3>",
  "<div class='figure-plan'><div class='label'>F1</div><div class='content'><p><strong>Approved.</strong> Keep the matched design and add the observed chromosome-number changes that motivate the Results claim.</p></div></div>",
  "<div class='figure-plan'><div class='label'>F3</div><div class='content'><p><strong>Simplify A.</strong> Show observed and predicted growth rate over passage only. The historical left-side control comparators in A and B are valid and remain.</p><p><strong>Remove E.</strong> No negative-control/ablation panel and no associated refit are active.</p></div></div>",
  "<div class='figure-plan'><div class='label'>F4</div><div class='content'><p><strong>Keep the science in main.</strong> Fit adequacy comes first, followed by latent O2, the low/intermediate/high feature result, and the solution landscape.</p><p><strong>Preserve the pooled embedding.</strong> Historical A–E are in-vivo analyses, while historical F deliberately includes both in-vivo and in-vitro fits in one reference landscape. Keep its point universe, saved geometry, context markers, and overlays unchanged.</p></div></div>",
  "<div class='figure-plan'><div class='label'>F5</div><div class='content'><p><strong>Reorder.</strong> Fit adequacy → full all-six parameter ratios → proliferation → missegregation → survival.</p><p><strong>Use one universe.</strong> All panels come from the same six approved July pair winners.</p></div></div>",
  "<div class='figure-plan'><div class='label'>F6</div><div class='content'><p><strong>Deferred.</strong> Park the figure, its framing, and its missing analytical grid. It is not part of this five-figure drafting plan.</p></div></div>",
  "<div class='boundary'><strong>Meaning of “confidence.”</strong> The available evidence supports in-sample observed–predicted fit adequacy, optimizer-start rankings, and sensitivity across competitive solutions. It does not support formal confidence intervals, held-out predictive validation, posterior uncertainty, or unique parameter identification.</div>",
  "</section>",

  "<section id='disposition'><h2>3. Existing-Panel Disposition</h2>",
  "<p>", badge("targeted fix", "good"), " ",
  badge("replace", "stop"), " ",
  badge("defer", "warn"),
  " Useful evidence roles are retained, incompatible source universes are separated, and rejected/deferred work is removed from the active path.</p>",
  table_html(
    disposition,
    c("artifact_id", "current_role", "disposition", "target_role", "rationale", "required_action")
  ),
  "<details><summary>Disposition counts</summary>",
  table_html(disposition_counts, names(disposition_counts), "compact"),
  "</details></section>",

  "<section id='ideas'><h2>4. Revised Working Figure Architecture</h2>",
  "<article class='option recommended'><h3>Approved five-figure core</h3>",
  "<p>This is no longer a choice among the earlier three broad alternatives. It incorporates the approved Figure 1, removes rejected Figure 3E, preserves the requested Figure 4 science, reorders Figure 5, and parks Figure 6.</p>",
  "<p><strong>Design rule:</strong> Figures 3–5 each establish how well the model reproduces the fitted data before interpreting latent variables, fitted functions, parameter differences, or solution regions.</p></article>",

  "<div class='figure-plan'><div class='label'>F1</div><div class='content'>",
  "<p><strong>Matched systems and motivating observations.</strong></p>",
  "<p><strong>A</strong> simplified matched design and sampling timeline; <strong>B</strong> observed in-vitro chromosome-number trajectory across informative passages; <strong>C</strong> observed in-vivo starting versus terminal chromosome distributions with compact burden/harvest context.</p>",
  "<p class='micro'>The in-vivo evidence is a start/end comparison, not longitudinal tumor karyotyping.</p></div></div>",

  "<div class='figure-plan'><div class='label'>F2</div><div class='content'>",
  "<p><strong>Mechanistic model overview.</strong></p>",
  "<p>One integrated schematic linking context/resource input, chromosome-number-dependent proliferation and death, stress-linked missegregation and WGD, post-missegregation survival, and the predicted chromosome-state distribution.</p>",
  "<p class='micro'>Encode observed inputs, latent states, fitted functions, transitions, and outputs differently.</p></div></div>",

  "<div class='figure-plan'><div class='label'>F3</div><div class='content'>",
  "<p><strong>In-vitro fit adequacy and inferred reshaping mechanism.</strong></p>",
  "<p><strong>A</strong> observed/predicted growth rate over lineage passage only, with 2N/4N and control/deprived facets; <strong>B</strong> predicted state distributions with observed karyotypes overlaid and control/deprived paired; <strong>C</strong> fitted survival function; <strong>D</strong> exact-fit nonviable-daughter fraction versus missegregation rate across reference ploidies.</p>",
  "<p class='micro'>No E. Controls are fitted baselines, not held-out validation, and end before the deprived branches. Passage is not literal calendar time. C–D show how the selected fit operates; they do not prove mechanistic necessity.</p></div></div>",

  "<div class='figure-plan'><div class='label'>F4</div><div class='content'>",
  "<p><strong>In-vivo fit adequacy, resource regimes, and solution landscape.</strong></p>",
  "<p><strong>A</strong> seed25 observed/predicted tumor burden plus terminal chromosome distributions; <strong>B</strong> target and latent effective-O2 trajectories; <strong>C</strong> compact O2=0/1/5 discrimination triptych; <strong>D</strong> unchanged pooled in-vivo/in-vitro embedding; <strong>E</strong> focused <code>p_mis_base</code> and <code>n_O</code> distributions.</p>",
  "<p class='micro'>Panel D keeps both contexts, the existing point universe, saved coordinates, context symbols, solution overlays, and geometry. It must not be filtered or recomputed as an in-vivo-only embedding.</p>",
  "<p class='micro'>Omit necrosis from A because exported predictions are all NA. Keep full feature screens and objective diagnostics supplementary, while retaining the triptych and landscape in main.</p></div></div>",

  "<div class='figure-plan'><div class='label'>F5</div><div class='content'>",
  "<p><strong>Joint-fit adequacy and context-specific functions.</strong></p>",
  "<p><strong>A</strong> direct in-sample fit evidence for both contexts; <strong>B</strong> full all-six log-ratio view for every soft-coupled parameter; <strong>C</strong> proliferation; <strong>D</strong> stress-linked missegregation; <strong>E</strong> post-missegregation survival.</p>",
  "<p class='micro'>B–E share the same six-winner source universe. Across-winner bands are solution-sensitivity envelopes, not confidence intervals.</p></div></div>",

  "<h3>Approved fit-quality design</h3>",
  "<div class='notice'><strong>Approved encoding.</strong> Plot observations directly against fitted values, use a small seed/rank annotation, and move the full objective distributions and optimizer diagnostics to supplementary panels. An objective bar by itself is not fit quality; packing every diagnostic into the first panel would recreate the historical Figure 3A density problem.</div>",
  table_html(
    fit_quality,
    c("figure", "recommended_main_evidence", "objective_or_rank_evidence", "robustness_evidence", "known_data_issue", "interpretive_limit")
  ),

  "<h3>Evidence-level corrections for Figure 4</h3>",
  "<div class='cards'>",
  "<article class='option'><h3>Low O2</h3><p><code>mu_hp</code>, the hypoxia-associated high-ploidy death-strength term, is the strongest reproduced one-feature separator at O2=0 (AUC ≈ 0.849).</p></article>",
  "<article class='option'><h3>High O2</h3><p><code>p_mis_base</code>, baseline per-chromosome missegregation probability, is strongest at O2=5 (AUC ≈ 0.903). The buffering feature <code>buffer_beta</code> is lower-ranked (AUC ≈ 0.672).</p></article>",
  "</div>",
  "<div class='warning'><strong>Interpretation limit.</strong> These AUCs discriminate model-defined lower- and higher-ploidy fixed-O2 attractor modes. They are not causal effect sizes for endpoint ploidy.</div>",
  "</section>",

  "<section id='decisions'><h2>5. Gate A Resolution</h2>",
  "<p>Gate A is closed. The following decisions form the approved drafting contract.</p>",
  "<ol class='checklist'>",
  "<li><strong>Figure 4 solution regions:</strong> retain all three already displayed in-vivo regions: <code>vi_C01</code>, <code>vi_C02</code>, and <code>vi_C03</code>. No two-regime remapping is required.</li>",
  "<li><strong>Fit-quality density:</strong> use compact direct observed–predicted blocks in Figures 3–5, with objective distributions and optimizer diagnostics supplementary.</li>",
  "<li><strong>Architecture:</strong> draft the approved five-figure set described above.</li>",
  "</ol>",
  "<div class='notice'><strong>Drafting authorized within bounds.</strong> Gate A approval does not authorize optimizer refits, production manuscript edits, resolution of the necrosis export, revival of Figure 6, or a new analytical grid.</div>",
  "</section>",

  "<section id='appendix'><h2>6. Appendix: Sources, Caveats, And Deferred Work</h2>",
  "<h3>Source inventory</h3>",
  "<p>The July fit bundle contains the needed scientific inputs without refitting, but it is an ignored/local export. Drafting must promote the minimal panel-ready tables, hashes, and generators instead of depending indefinitely on the full local bundle.</p>",
  table_html(
    inventory,
    c("source_id", "figure_scope", "path", "status", "reuse_class", "requires_refit", "required_action", "caveat")
  ),
  "<h3>Key fit-quality caveats</h3>",
  "<ul>",
  "<li><strong>Figure 3:</strong> seed10 is the top in-vitro objective, but 12 starts lie within 0.01 objective units, the local optimizer ended with code 1, five active parameters are at bounds, the stored table records 0/500 DEoptim convergence flags, and no Hessian-based uncertainty artifact exists.</li>",
  "<li><strong>Figure 4:</strong> seed25 is the best total weighted-MAP objective among 500 starts, yet ranks 52nd for ploidy and 189th for burden. Twenty-nine fits are within 1% and 241 within 5% of the best objective.</li>",
  "<li><strong>Figure 5:</strong> the six winners span total objectives 18.852–19.978. Stored checks validate selection provenance and arithmetic, not held-out prediction; no L-BFGS-B refinement was accepted.</li>",
  "<li><strong>Necrosis:</strong> separate and joint exports contain observed endpoints but all saved predicted necrosis values are NA. Necrosis is excluded from proposed fit-quality panels pending source resolution.</li>",
  "</ul>",
  "<h3>Evidence and disposition summaries</h3>",
  "<div class='cards'><div>",
  table_html(status_counts, names(status_counts), "compact"),
  "</div><div>",
  table_html(figure_counts, names(figure_counts), "compact"),
  "</div></div>",
  "<h3>Restored manuscript build context</h3>",
  "<div class='warning'><strong>The TeX source is restored and identified.</strong> Its current SHA-256 is <code>2fcc918e6d095a5b2962eac9bc8d2b581ba0d15a724fb6f975c61bd21e3ff3b0</code>. The baseline build still lacks referenced figure assets, two parameter-table TeX files, and the bibliography; these do not block ideation but do block a manuscript-ready compile.</div>",
  "<details><summary>Historical context manifest</summary>",
  table_html(
    context_manifest,
    c("context_id", "source_commit", "source_path", "sha256", "width_px", "height_px", "use_boundary")
  ),
  "</details>",
  "<details><summary>Explicitly deferred work</summary><ul>",
  "<li>Figure 6 response-curve framing, analytical grid, and selection overlay.</li>",
  "<li>Journal-specific dimensions, main-figure limit, and color constraints until polishing.</li>",
  "<li>Any restricted-model refit or negative-control Figure 3E.</li>",
  "<li>Formal uncertainty claims not supported by the present artifacts.</li>",
  "</ul></details>",
  "<p><strong>Feedback-manager context:</strong> the July 24 report review is archived under <code>feedback_archive/20260724_ideation_report_review</code> and served as <code>SPAN-20260724-001</code>. The earlier July review remains preserved under <code>SPAN-20260703-001</code>.</p>",
  "<p><strong>Companion files:</strong> <code>ideas.md</code>, <code>feedback_decisions.csv</code>, <code>fit_quality_evidence.csv</code>, <code>existing_panel_disposition.csv</code>, <code>tex_figure_requirements.csv</code>, <code>source_inventory.csv</code>, <code>baseline_build.md</code>, and <code>feedback_manager_context.md</code>.</p>",
  "<footer>Self-contained approved ideation report generated from repository-local manifests and embedded historical context. Gate A is closed; bounded figure drafting is authorized.</footer>",
  "</section></main></body></html>"
)

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
writeLines(html, con = output_path, useBytes = TRUE)
cat("Wrote", output_path, "\n")
