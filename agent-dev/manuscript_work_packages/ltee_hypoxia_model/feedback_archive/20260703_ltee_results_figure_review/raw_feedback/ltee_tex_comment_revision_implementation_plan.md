# LTEE Hypoxia Results Revision Plan

## Source Used For Line Mapping

Target manuscript source: `ltee_hypoxia_model.tex` in the repository root.

Line references below are based on the current working-tree file as inspected on 2026-07-03. If the manuscript is edited before implementation, re-run:

```bash
nl -ba ltee_hypoxia_model.tex | sed -n '83,217p'
```

The relevant current regions are:

- Results opening: lines 83-113.
- Figure 1 caption: lines 121-128.
- Duplicate model overview figures: lines 130-154.
- In vitro figure caption: lines 156-169.
- In vivo figure caption: lines 171-187.
- Joint-fit Results section: lines 109-113.
- Joint figure caption: lines 189-204.
- Fixed-O2 response-class figure caption: lines 206-217.

## Transcript Sentence To Line Map

| ID | Transcript sentence or comment | Current line(s) | Required interpretation |
|---|---:|---:|---|
| T01 | "Matched SUM-159 lineages motivate a resource stress model of ploidy evolution." | 87 | Current subsection title. |
| T02 | "I don't like that SUM-159 is in the title of this subsection." | 87 | Revise title to remove SUM159/SUM-159. |
| T03 | "We first used the matched SUM-159 near-diploid and near-tetraploid lineages to place the culture and tumor experiments on a common interpretive axis." | 89 | Existing opening sentence; keep concept but do not lead with it. |
| T04 | "The same isogenic lineage pair was observed in vitro and in vivo with measured culture passages, tumor burden trajectories, tumor chromosome number profiles, and histologic necrosis measurements aligned in a data stream overview (Figure 1A)." | 89 | Existing data-stream sentence; keep after the biological observation. |
| T05 | "This design allowed differences in ploidy evolution to be interpreted as environment-dependent selection, rather than as differences between unrelated model systems." | 89 | Existing interpretive sentence; keep after the observed ploidy behavior. |
| T06 | "This data motivated a model in which resource limitation creates opposing pressures on chromosome number evolution." | 91 | Current model-motivation sentence; revise because the motivation should be the observed ploidy reduction and in vitro WGD sequence. |
| T07 | "Okay." | 89-91 | Transition only; no manuscript action. |
| T08 | "I think that what motivated a model, such a model, is a specific observation in both of these datasets, namely the observation of a ploidy reduction." | 89-91 | Add ploidy reduction as the lead observation. |
| T09 | "Surprisingly, in the in vitro setting, the ploidy reduction was preceded by whole genome doubling." | 89, 95-97 | Add WGD-before-reduction to first Results paragraph and in vitro interpretation. |
| T10 | "So that's interesting for their similarities but also differences between the systems." | 89-91 | Frame shared ploidy reduction plus in vitro-specific WGD as the biological rationale. |
| T11 | "So that needs to be included in the first paragraph." | 89 | Rewrite paragraph 1 around ploidy reduction and WGD. |
| T12 | "That's the most interesting part of the experiment, these observations." | 89 | Make the observed trajectories the narrative lead. |
| T13 | "So, then going to the second paragraph, this data motivated a model in which resource limitation creates opposing pressures on chromosome number evolution." | 91 | Rewrite paragraph 2 as model motivation after the observation. |
| T14 | "In the conceptual model, resource stress can oppose haploidy by imposing growth and death costs on chromosome-rich states (Figure 3B), while also increasing sen and whole genome doubling pressure and generating deleterious selection (X, Figure 3C)." | 91 | Current text already uses high ploidy/CIN, not haploidy/sen. Keep concept but correct figure citation after deduplication. |
| T15 | "Haploidy can partly buffer chromosome missegregation events that would be lethal or strongly deleterious near diploidy (Figure 3D), thus the predicted ploidy outcome depends on the balance among growth suppression, stress-associated death, sen generation, whole genome doubling generation, and post-missing survival filtering." | 91 | Current text already uses high ploidy/CIN/post-missegregation. Keep concept; verify no speech-to-text artifacts remain. |
| T16 | "Okay, so Figure 3 should not be cited at all in this first section." | 91, 130-154 | Resolve figure-number duplication before final citation updates. |
| T17 | "Figure 2 should be cited here instead." | 91, 130-154 | Cite the single retained model-overview figure. |
| T18 | "Okay. Let's read the captions for Figures 1 and 2 first, next I mean." | 121-154 | Transition to caption cleanup. |
| T19 | "Figure 1, matched SUM-159 lineages reveal environment-dependent ploidy selection." | 125 | Figure 1 title; not specifically rejected. |
| T20 | "Existing optimization, in A we have our existing optimization data stream overview showing the matched SUM-159 near-diploid and near-tetraploid lineages in vitro and in vivo on a comparable experimental time axis." | 126 | Remove "Existing optimization" / workflow phrasing; keep biological description. |
| T21 | "The in vitro portion shows control and two hypoxia/anoxia passage histories, biological replicates, with known target oxygen encoded by blue intensity and measurement events marked by flow cytometry and karyotyping." | 126 | Keep biological content, phrase cleanly. |
| T22 | "The in vivo portion shows tumor burden trajectories for the fitted untreated harvest mice, with tumor ploidy and histologic necrosis measurements indicated at harvest." | 126 | Keep biological content, phrase cleanly. |
| T23 | "This panel was copied from the automatically generated data stream figures." | 126 | Remove entirely. |
| T24 | "No, this should not be in here." | 126 | Remove provenance text. |
| T25 | "We shouldn't state where this was copied from." | 126, 163, 181-185, 198-202, 212 | Apply globally to all captions: no report names, file names, extraction provenance, or workflow notes in manuscript captions. |
| T26 | "This is not, this should be removed from the entire draft." | 126 and all captions | Global provenance cleanup. |
| T27 | "And then the last sentence of Figure 1 caption is kind of, the second part of the section might fit into the results subsection 1, but otherwise I definitely would not keep it." | 126 | Remove final interpretive sentence from caption; optionally reuse core idea in Results paragraph only if needed. |
| T28 | "The same isogenic system underlies both settings, the two settings span different time frames and oxygen is known by passage in vitro." | 126 | Current final caption sentence; remove from caption. |
| T29 | "No, we can remove this entire sentence entirely from the draft." | 126 | Delete final interpretive sentence rather than moving it. |
| T30 | "Okay, on Figure 2, resource limitation creates opposing pressures on ploidy evolution." | 134, 147 | Figure 2 title; revise as model overview, not a proven result. |
| T31 | "Panel A, central model schematic, resource limitation shapes ploidy evolution through competing forces." | 135, 148 | Keep if model-overview figure retained. |
| T32 | "Panel B, cost arm, haploidy is costly under resource stress through growth depression." | 136, 149 | Keep concept using "high ploidy", not haploidy. |
| T33 | "Panel C, instability arm, stress can increase whole genome doubling pressure." | 137, 150 | Keep concept; use CIN/WGD. |
| T34 | "Panel D, buffering arm, haploidy buffers otherwise lethal missegregations." | 138, 151 | Keep concept using high ploidy. |
| T35 | "Panel E, integrated outcome. The outcome depends on growth suppression, death, sen, whole genome doubling generation, and post-missegregation survival filtering." | 139-140, 152 | Keep concept; use CIN, WGD, post-missegregation survival. |
| T36 | "This is good, except the title should convey that this is an overview of the model." | 134, 147 | Rename figure title to "Model overview..." |
| T37 | "The title currently reads as if the figure proves that resource limitation creates opposing pressures on ploidy evolution." | 134, 147 | Avoid causal/proof framing. |
| T38 | "But that's not the case. These are the assumptions of the model, essentially." | 134, 147 | Caption should state model mechanisms/assumptions. |
| T39 | "All right. In vitro ploidy, next section, next subsection of the results." | 93 | Transition; no direct edit except section review. |
| T40 | "In vitro ploidy evolution during hypoxia occurs without strong hypoploid elimination." | 93 | Current title already says high-ploidy; keep or lightly polish. |
| T41 | "The protein vitrofit supported a culture regime in which oxygen stress reshapes ploidy without requiring strong elimination of hypoploidy cells." | 95 | Current text already says "Separate in vitro fits" and "high-ploidy"; no artifact remains. |
| T42 | "The in vitro fit jointly tracked proliferation, chromosome count observation, and live-dead burden across the 204 control and oxygen-deprived branches, as in Figure 4A." | 95 | Current text says 2N and 4N; keep. Ensure figure reference remains correct after numbering. |
| T43 | "Oxygen slowed proliferation and altered ploidy trajectories, but it did not require a large net death or strong direct removal of hypoploidy states." | 95 | Current text already matches with high-ploidy. |
| T44 | "The resulting mechanism is therefore not that more a cells are killed off or removed simply, but that hypoploidy cells provide more opportunities for chromosome missegregation." | 97 | Revise mechanism language to high-ploidy cells missegregate more often because they carry more chromosomes. |
| T45 | "And the model-implied survival filter allows hypoploidy states to tolerate some chromosome loss products, as in Figure 4C." | 97, 165 | Current concept present; keep as high-ploidy survival buffering. |
| T46 | "This buffering is double-edged." | 97 | Keep, but revise explanation. |
| T47 | "It preserves hypoploidy cells after missegregation, but it also makes some chromosome loss daughters viable." | 97 | Replace "preserves high-ploidy cells after missegregation" framing. |
| T48 | "And I'm not sure if this is framed right." | 97 | Signals required rewrite. |
| T49 | "Ploidy's double-edged should not be that it preserves hypoploidy cells after missegregation." | 97 | Remove preservation framing. |
| T50 | "It should be that hypoploidy cells missegregate more chromosomes, but hypoploidy also makes some chromosome loss daughters viable, gearing the population route back toward lower chromosome number." | 97, 164-166 | Main rewrite target. Use high-ploidy, not hypoploidy. |
| T51 | "That's in Figure 4BD." | 97, 164, 166 | Current manuscript calls this Fig. 3B/D within label `fig:iteration1-invitro-ploidy-reshaping`; restore panels and cite B-D. |
| T52 | "Direct chromosome burden effects may still bias which daughters expand, but the in vitro fit supports making the missegregation kernel plus ploidy-dependent survival the central explanation, with strong hypoploidy elimination treated as a rejected or secondary interpretation, in Figure 4E." | 97, 167 | Keep the idea and keep panel E present as an explicit panel-to-be-generated placeholder if direct evidence is not available yet. |
| T53 | "Let's read Figure 4." | 156-168 | Transition to in vitro figure caption. |
| T54 | "In vitro ploidy reshaping occurs without strong hypoploidy elimination." | 162 | Figure title; current high-ploidy wording is acceptable. |
| T55 | "I just realized Figure 3 is actually Figure 2." | 130-154 | There is a duplicated model-overview figure block. |
| T56 | "So the citations in the first result sections were right after all." | 91 | After deduplication, ensure the model paragraph cites the retained model overview. |
| T57 | "There's just a duplication." | 130-154 | Remove one duplicate model-overview block. |
| T58 | "Figure, what is currently Figure 2, should be removed entirely." | 130-154 | Consolidate to a single model-overview figure. Decide which block to retain by label/citations before editing. |
| T59 | "Okay, so back to Figure 4." | 156-168 | Transition. |
| T60 | "In vitro ploidy reshaping occurs without strong hypoploidy elimination." | 162 | Keep high-ploidy title. |
| T61 | "In vitro fit summary panel." | 163 | Panel A caption. |
| T62 | "This panel was extracted from an HTML report and aligns growth rate observations, chromosome count observations, and live-dead burden." | 163 | Remove extraction/report provenance; keep what the panel shows. |
| T63 | "It supports the fit summary panel by showing that the culture fit jointly tracks proliferation and ploidy across passages rather than fitting either stream alone." | 163 | Keep as manuscript content if concise. |
| T64 | "Okay, so we need to remove information as to where this figure was extracted from." | 163 and all captions | Global caption cleanup. |
| T65 | "We don't need reference to HTML reports here or some PNG files." | 124, 159-161, 174-179, 192-196, 209-212 | Remove manuscript-facing report/file references. Comments may remain only if useful for build provenance, but not inside captions. |
| T66 | "Panel 3B, missegregation source." | 164 | Restore panel B; do not leave missing. |
| T67 | "And the internal panel should show that hypoploidy cells missegregate more often." | 164 | Panel B should show high-ploidy cells have more missegregation opportunities/events. |
| T68 | "And survivor buffering shows model viability after one chromosome missegregation event as a function of chromosome number, decreasing curve visualizes the fitted buffer, with which higher ploidy states better tolerate missegregation products." | 165 | Caption panel C; verify curve orientation against actual plot, but concept is higher ploidy better tolerates missegregation. |
| T69 | "Figure 3D missing downward reshaping." | 166 | Restore panel D; do not leave missing. |
| T70 | "The internal panel should show surviving chromosome loss daughters generating viable lower ploidy descendants." | 166 | Panel D should represent downward route. |
| T71 | "No distribution heatmap showed ploidy reshaping over time, but did not directly show the cellular route." | 166 | Earlier assessment to avoid overclaiming, later superseded by user disagreement; if using a heatmap, caption as trajectory-level evidence. |
| T72 | "Figure 3E missing negative control rejected interpretation." | 167 | Panel E still lacks direct evidence, so keep it present as a panel-to-be-generated placeholder rather than deleting it. |
| T73 | "The internal panel should make clear that the in vitro experience is not simply that diploid cells are killed off." | 167 | If retained, panel E should visualize rejected high-ploidy elimination interpretation. |
| T74 | "And no current iteration 1 panel directly tests or visualizes the rejected interpretation." | 167 | If still true, keep panel E in the manuscript caption as a panel-to-be-generated placeholder so the needed figure remains visible. |
| T75 | "Okay, so I disagree." | 164, 166 | Supersedes missing-panel assessment for B and D. |
| T76 | "Figure 3B, Figure 3C, sorry, I disagree Figure 3B and Figure 3D are not missing." | 164, 166 | Restore B and D panels. |
| T77 | "I think I liked the content you had there earlier in iteration 1." | 164, 166, `oxygen/figures/iteration1` | Recover prior iteration-1 content. |
| T78 | "Put it back." | 164, 166, assembled Fig. 3 | Restore files and rerun assembler. |
| T79 | "So Figure 3B, sorry, Figure 4B and Figure 4E, missegregation source and downward reshaping, they should be put back." | 164, 166 | Restore the missegregation-source and downward-reshaping panels, regardless of numbering confusion. |
| T80 | "There was also a figure under iteration 1 at some point showing the models ability to recapitulate in vitro WGD as well." | 95-97, 162-168, iteration-1 files | Restore or regenerate WGD-recapitulation panel. |
| T81 | "That should be mentioned in the main text and the figure should be put back into the assembled version if not already mentioned in the plan." | 95-97, 158, 162-168 | Main text and assembled in vitro figure must explicitly include WGD recapitulation. |
| T82 | "All right, let's move on to..." | 99 | Transition to in vivo section. |
| T83 | "Let's go to the in vivo fits, identify resource-dependent ploidy regimes and unstructured parameter landscapes." | 99 | Current title says structured parameter landscapes; review title wording, but no explicit requested change besides section review. |
| T84 | "For the in vivo fits, oxygen should be interpreted as the formal model variable and as a first-order proxy for proliferative effective resources." | 101 | Current text matches; keep. |
| T85 | "Limited effective oxygen trajectories provide a latent model-implied resource axis rather than direct oxygen measurements." | 101 | Current text matches with fitted effective O2; keep. |
| T86 | "Figure 5a, this distinction is important because in vivo hypoxia is not equivalent to in vitro low oxygen in rich medium." | 101, 181 | Keep distinction; update figure numbering only via labels. |
| T87 | "Tumor stress can also reflect nutrient limitation, perfusion heterogeneity, waste buildup, cell density, stroma constraints, and tissue structure." | 101 | Current text matches; keep. |
| T88 | "Fixed oxygen analysis showed that the dominant ploidy state depends on resource regime." | 103 | Keep result but define analysis first. |
| T89 | "At low oxygen, feature AUC analysis supported a death-growth-dominated interpretation in which stress-associated death terms provided the most compact explanation of lower versus higher ploidy modes." | 103 | Keep after defining feature-AUC and mode. |
| T90 | "So the reader will not understand what feature AUC analysis is." | 103 | Add explanatory sentence before result. |
| T91 | "We need to explain this in a sentence or two what we did, how many parameters, and why we did this, rather than using some sort of undefined terms like fixed O2 analysis." | 103 | Add concise methods-in-results setup; verify exact parameter count from code/report before writing. |
| T92 | "At intermediate oxygen, mode identity was less reducible to a single parameter because of a trade-off among death selectivity, misversivity, buffering strength, oxygen thresholding, and baseline missegregation." | 103 | Keep but define mode; correct "misversivity" to intended term, likely missegregation sensitivity/selectivity. |
| T93 | "Again, this, it didn't define what mode is here." | 103 | Define mode before low/intermediate/high O2 statements. |
| T94 | "At high oxygen, reduced death pressure made baseline missegregation and buffering-related parameters more important because they determine how often buffering capacity is tested during growth." | 103 | Keep after definitions. |
| T95 | "The in vivo parameter landscape analysis further suggested that multiple fitted solution regions may represent distinct biological regimes rather than purely optimizer noise." | 105 | Remove "rather than purely optimizer noise." |
| T96 | "Remove then rather than purely optimizer noise." | 105 | Delete phrase. |
| T97 | "Cluster-associated differences in interpretable fitted axes, including basal missegregation and oxygen response shape, connect the solution regions to biological mechanisms." | 105 | Keep concept. |
| T98 | "Remove the last part, rather than only to embedding coordinates." | 105 | Delete phrase. |
| T99 | "This interpretation should remain primary." | 105 | Keep biological interpretation as primary. |
| T100 | "Clustering and dimensionality reduction need to account for parameter transformation scale and prior bound structure." | 105 | Move earlier into explanation of how clustering was done. |
| T101 | "Okay, all this last part here should be moved upwards, because again, the reader will not understand how this analysis was done." | 103-105 | Add setup before interpreting clusters. |
| T102 | "What cluster, you're talking about clusters associated differences, but there is no mentioning how the clustering was done." | 105 | Explain what was clustered and how. |
| T103 | "So there needs to be a short couple of sentences describing what was clustered." | 105 | Add cluster-method setup sentence. |
| T104 | "And then the whole dimensionality reduction, parameter transformation scale, prior bound structure, that should all be part of that explanation." | 105 | Include prior-aware preprocessing in setup. |
| T105 | "As a complementary model selection view, fixed oxygen dominant ploidy response curves provided a richer phenotype than binary mode labels." | 107 | Remove from the in vivo section, but reintroduce later as the Results-ending model-selection logic after the joint-fit axes. |
| T106 | "The ploid response class analysis example showed complex non-monotonic, increasing U-shaped, and so on, behaviors across candidate fits, Figure 7a-c." | 107, 211-215 | Move/reframe into the final Results section using the current Figure 6 label/block if retained. |
| T107 | "Oh, no, Figure 6a-c, I guess." | 107, 211-215 | Current manuscript uses Figure 6-style response-class block lines 206-217; final numbering should be label-based. |
| T108 | "These classes can be used to prioritize parameter regimes whose oxygen ploidy response resembles empirically plausible behavior." | 107 | Preserve this idea for the final Results ending, not in the in vivo section. |
| T109 | "But the final selection should combine response class fit quality and the parameter landscape context, rather than relying on curve shape alone." | 107 | Superseded by later comment: use fixed-O2 response curves to select among remaining explanations; still avoid relying on curve shape alone if fit/landscape context is available. |
| T110 | "I think I want to leave out this entire paragraph, the last paragraph of the in vivo section can stay out." | 107 | Delete line 107 from the in vivo section; do not discard the response-curve concept entirely. |
| T111 | "Moving on to, um, well, let me do the caption first for in vivo figures." | 180-185 | Signals in vivo caption will need review next; current transcript stops before detailed caption comments. |
| T112 | "Joint fits distinguish in vivo and in vitro stress response mechanisms." | 109 | Current final Results subsection title; keep or lightly align with caption wording. |
| T113 | "Joint soft coupled fits were used to ask which mechanisms were shared across culture and tumor settings, and which remained context-specific despite regularization." | 111 | Current opening sentence; keep with minor punctuation/style cleanup. |
| T114 | "Figure 6A, the strongest, 6A, figure 5A now." | 111, 197-203 | Numbering confusion; maintain label-based references to the joint-fit figure. |
| T115 | "The strongest interpretation is that in vivo differs from in vitro less in the existence of oxygen stress itself, and more in how that stress is translated into proliferation, chromosome instability, and post-mitosis survival filtering." | 111 | Current sentence says post-missegregation survival filtering; keep that correction. |
| T116 | "We could add here presumably because it's not just oxygen stress, because it's oxygen stress plus all the other stressors that are present in vivo." | 111 | Add caveat that in vivo stress reflects oxygen plus additional tumor-environment stressors. |
| T117 | "Three context-specific axes emerged from the joint fit summaries." | 113 | Current paragraph opener; keep. |
| T118 | "First, in vivo context showed a lower effective proliferation ceiling, consistent with a broader growth constraint in tumors than in culture, figure 6B." | 113, 199 | Current text matches; keep label-based citation. |
| T119 | "Second, stress-linked chromosome mis-segregation was stronger in vivo, suggesting tighter coupling between resource stress and syn generation in tumor settings, figure 6C." | 113, 200 | Current text uses CIN generation; keep CIN, not speech-to-text `syn`. |
| T120 | "Third, the post-mitosis survival filter was more ploidy-dependent in vivo, indicating that the tumor environment changes not only how many altered daughters are generated, but also which altered progeny persists, figure 6D." | 113, 201 | Current text uses post-missegregation survival; keep and polish grammar to `persist`. |
| T121 | "The integrated progeny ratio summary therefore supports a context comparison branch centered on proliferative ceiling, stress-to-syn coupling, and survival filter differences, figure 6E." | 113, 202 | Current final sentence should be removed; if Fig. E remains, caption can describe it, but Results should not end on this sentence. |
| T122 | "The last sentence here can be removed." | 113 | Delete final sentence of line 113 beginning `The integrated parameter-ratio summary...`. |
| T123 | "Instead, we want to bring back here the fixed O2 dominant ploidy response curves." | 107, 211-215, after 113 | Move/rewrite response-curve material after joint-fit axes. |
| T124 | "Essentially, how I envision this manuscript results section ending is that we will still have multiple potential explanations for, you know, the regimes that could explain the in vivo ploidy evolution." | after 113 | Add transition: joint/in vivo analyses leave multiple plausible regimes/explanations. |
| T125 | "And to select among them, we will then use the eigenvector solutions for fixed oxygen." | after 113 | Add fixed-O2 eigenvector/steady-state or dominant-eigenvector analysis as model-selection discriminator; verify exact term used in code before writing. |
| T126 | "And the essential expectation is that higher oxygen correlates with higher ploidy, and only a subset of solutions have this particular behavior." | after 113, 211-214 | Final Results ending should emphasize monotone increasing O2-ploidy behavior as the expected/biologically plausible response. |
| T127 | "So they're more likely to be the ones that are closer to reality." | after 113 | Phrase cautiously: solutions with higher-ploidy-at-higher-O2 behavior are prioritized as more biologically plausible, not definitively true. |
| T128 | "And here we can also essentially explain under which regimes, like which parameters predict a monotonically increasing relationship between ploidy, between fixed oxygen and ploidy." | after 113, 211-215 | Add explanation of which parameter regimes/features predict monotone increasing fixed-O2 ploidy response. |
| T129 | "A general comment is also that there are several ad hoc terms used that are not defined, fixed-O2 analysis is just one example." | 101-113, 180-215 | Add a global Results/caption pass: replace shorthand analysis labels with one or two explanatory sentences saying what was computed and why. |
| T130 | "Find all such instances and specify that they should be replaced with a sentence or two of what was done and why." | 101-113, 180-215 | Add the undefined/ad hoc term inventory below and require each term to be defined or rewritten before interpretation. |
| T131 | "For example, to understand whether oxygen influences steady state ploidy in a predictable way we..." | 103, after 113, 211-215 | Use this style for fixed-O2 response-curve language: state purpose, computational operation, and biological reason before using the term. |

## Undefined Or Ad Hoc Term Inventory

Use this inventory during implementation. For each term, either remove the shorthand or replace it with a concise purpose-first explanation. Keep Results prose short: define only enough for the reader to understand the claim, and refer technical details to Methods or Supplement, which can be written later. Do not rely on the label alone.

| Current term/phrase | Current line(s) | Problem | Required replacement strategy |
|---|---:|---|---|
| `culture regime` | 95 | Reads like an unexplained category. | Replace with what the fitted culture simulations show, e.g. oxygen stress slowed growth and shifted chromosome-number distributions without requiring strong direct high-ploidy killing. |
| `model-implied survival filter` | 97, 165, 201 | Technical model term introduced as if self-evident. | Add a short definition: the fitted function that maps a missegregated daughter/mother chromosome number to survival probability after chromosome loss/gain. |
| `missegregation kernel` | 97 | Mathematical shorthand not defined in Results. | Replace with "the model rule that determines how chromosome gains/losses are produced during division," or define before use. |
| `effective resource stress` | 101 | Could sound like a measured quantity. | State that in vivo O2 is a latent model variable used as a proxy for broader tumor resource stress, not a direct oxygen measurement. |
| `latent, model-implied resource axis` | 101, 181 | Dense wording. | Explain that the model infers an effective O2 trajectory from tumor growth and fitted supply-demand parameters, then uses it to order simulated tumor states by resource limitation. |
| `Fixed-O2 analyses` | 103 | Shorthand; reader will not know what was held fixed or why. | Replace with purpose-first wording, e.g. "To ask whether oxygen alone would predict the long-term chromosome-number composition for a fitted parameter set, we simulated/evaluated the model at fixed O2 values and compared the resulting dominant ploidy." Verify whether this used simulation, steady state, attractor, or eigenvector. |
| `dominant ploidy state` | 103, 107, 212 | Undefined outcome. | Define as the chromosome-number/ploidy class with the largest predicted fraction in the fixed-O2 steady-state/long-time distribution, or use the exact code definition. |
| `resource regime` / `low-O2 regime` / `intermediate-O2 regime` / `high-O2 regime` | 99, 103, 180-184 | Regime labels are not meaningful until the O2 grid/reference values are stated. | Introduce as the low, intermediate, and high fixed-O2 reference points used for the analysis; give the actual values if available. |
| `feature-AUC analyses` / `feature-AUC panel` | 103, 182-184 | Statistical shorthand; AUC target and features are undefined. | Replace with a sentence saying that fitted parameters or model-implied features were scored by how well they separated lower- versus higher-ploidy modes across candidate fits; AUC was used as a one-feature discrimination score. Put exact parameter/feature counts in Methods/Supplement unless a short parenthetical is essential. |
| `death/growth-dominated interpretation` | 103, 182 | Interpretive label without setup. | State which death/growth parameters or derived features had the highest discrimination scores and what biological process they represent. |
| `lower- versus higher-ploidy modes` / `mode identity` / `binary mode labels` / `ploidy mode` | 103, 107, 182 | Mode is undefined. | Define how modes were assigned, e.g. threshold or class labels based on the dominant chromosome-number solution. Use exact thresholds/classes from the code/report. |
| `parameter-landscape analyses` | 105, 185 | Vague; could mean embedding, clustering, or parameter comparison. | Replace with a concise purpose-first description: fitted parameter sets were compared after scale/prior-aware preprocessing to see whether good fits occupy distinct biological regions. Put input matrices, algorithms, and cluster settings in Methods/Supplement. |
| `fitted solution regions` | 105 | Ambiguous. | Define as groups/clusters of fitted parameter sets in the preprocessed parameter/phenotype space. |
| `cluster-associated differences` | 105, 185 | Cluster construction is not defined. | Precede with how clusters were generated, then describe the original-parameter or model-phenotype differences that distinguish clusters. |
| `interpretable fitted axes` | 105 | Too abstract. | Name the axes directly, e.g. baseline missegregation probability and oxygen-response steepness/threshold. |
| `baseline missegregation and oxygen-response shape` | 105, 185 | Needs biological translation. | Translate to parameter names and meanings, e.g. baseline chromosome missegregation rate and how sharply/at what O2 stress turns on. |
| `prior-aware` / `prior-aware z positions` | 105, 185 | Undefined preprocessing method. | Briefly explain that parameters were transformed and scaled relative to their prior/bound ranges so large numeric units do not dominate. Put the exact preprocessing recipe in Methods/Supplement. |
| `dimensionality reduction` / `embedding coordinates` / `t-SNE clusters` | 105, 185 | Methods labels without what/why. | State only the purpose in Results: visualization/clustering of fitted solutions after scale/prior-aware preprocessing. Put input feature matrices and algorithm settings in Methods/Supplement. Remove "rather than only to embedding coordinates" from manuscript prose. |
| `model-implied features such as net growth, turnover, mean O2, missegregation-associated death, and fixed-O2 attractor behavior` | 105 | List of derived quantities without definitions. | Either define each derived quantity briefly or move the list to Methods/Supplement. In Results, use only the derived features actually shown in the figure. |
| `fixed-O2 attractor behavior` | 105 | Undefined and possibly computationally inaccurate. | Verify exact computation; replace with fixed-O2 steady state, dominant eigenvector, long-time composition, or attractor only if that is the actual object. |
| `model-selection view` | 107, 211-215 | Abstract framing. | Replace with purpose-first sentence: after several fitted mechanisms remain plausible, fixed-O2 response curves are used to prioritize solutions whose predicted ploidy increases with oxygen. |
| `fixed-O2 dominant-ploidy response curves` | 107, after 113, 211-214 | Needs what/why. | Define as curves obtained by evaluating each fitted parameter set at a range of fixed O2 values and plotting the predicted dominant/steady-state ploidy versus O2. State the goal: test whether oxygen influences ploidy predictably. |
| `pooled response-class examples` / `response classes` | 107, 212-214 | Undefined classification. | Explain that curves were grouped by shape, such as monotone increasing, monotone decreasing, flat, U-shaped, or nonmonotone, and why the monotone increasing class is prioritized. |
| `monotone-increasing priority class` | 214 | Jargon. | Replace with "solutions in which predicted ploidy increases as fixed O2 increases," then explain why this matches the biological expectation. |
| `selection overlay` | 215 | Figure-planning jargon. | Replace or remove from manuscript-facing caption. If retained in plan, define it as overlaying response-curve class, fit quality, and parameter-landscape position. |
| `joint soft-coupled fits` / `regularization` | 111, 197-202 | Methods shorthand; probably acceptable after Methods but should still be readable. | Add a compact reminder: paired in vivo/in vitro parameters were penalized for unnecessary divergence, so divergence indicates context-specific support from the data. |
| `context-specific axes` | 113 | Abstract. | Immediately name the three axes: proliferative ceiling, stress-to-CIN coupling, and post-missegregation survival. |
| `effective proliferative ceiling` | 113, 199 | Technical parameter concept. | Translate as maximum fitted division/proliferation capacity under favorable oxygen, or name the model parameter if appropriate. |
| `stress-linked chromosome missegregation` / `CIN axis` | 113, 200 | Needs biological definition. | Explain as the fitted increase in chromosome missegregation probability under oxygen/resource stress. |
| `post-missegregation survival filter` / `survival-filter axis` | 113, 201 | Needs definition. | Define as the ploidy-dependent probability that daughters survive after chromosome missegregation. |
| `integrated parameter-ratio summary` / `integrated contrast` / `context-comparison branch` | 113, 202 | Abstract and the sentence is marked for removal. | Remove from Results final sentence. If the figure panel remains, caption it as a plot of fitted in vivo/in vitro parameter ratios and explain what ratios above/below one mean. |

## Implementation Plan

### Cross-Cutting Rule: Keep Analysis Definitions Succinct

Before revising each Results paragraph or caption, check the `Undefined Or Ad Hoc Term Inventory` above. Any shorthand analysis label should be replaced by one concise sentence or clause that gives the reader just enough context to follow the claim. The Results should answer:

1. What was computed?
2. Why was it useful for the biological question?

Technical details should be deferred to Methods/Supplement, including exact input matrices, preprocessing steps, parameter counts, thresholds, algorithms, and implementation details. If those Methods/Supplement sections are not written yet, the Results can still contain a forward reference such as `Methods; Supplementary Methods`.

Example style:

```text
To ask whether oxygen alone would push fitted models toward predictable long-term ploidy states, we evaluated each fit across fixed O2 levels and grouped the resulting O2-ploidy curves by shape (Methods; Supplementary Methods). Solutions in which predicted ploidy rises with oxygen match the expected resource-permissive behavior.
```

Use exact computational terms after verifying the code/report. Do not write `fixed-O2 eigenvector solution` unless the code actually uses a dominant eigenvector; otherwise use `steady state`, `attractor`, `long-time composition`, or the exact implemented term.

### Step 1: Establish The Editing Target And Preserve Line Context

Actions:

1. Treat root `ltee_hypoxia_model.tex` as the target source, not the removed `docs/ltee_hypoxia_model.tex`.
2. Before editing, capture current relevant lines:

```bash
nl -ba ltee_hypoxia_model.tex | sed -n '83,217p'
```

Validation:

```bash
test -f ltee_hypoxia_model.tex
rg -n "\\section\\*\\{Results\\}|Matched SUM159|In vitro ploidy|In vivo fits" ltee_hypoxia_model.tex
```

### Step 2: Rewrite The First Results Subsection Around The Observation

Targets:

- Title at line 87.
- Paragraph 1 at line 89.
- Paragraph 2 at line 91.

Actions:

1. Change the subsection title so it does not contain SUM159/SUM-159.
   - Candidate: `Matched culture and tumor trajectories motivate a resource-stress model of ploidy evolution`.
2. Rewrite the first paragraph to lead with:
   - ploidy reduction in both in vivo and in vitro datasets;
   - in vitro ploidy reduction being preceded by WGD/high-ploidy expansion;
   - why the matched near-diploid and near-tetraploid SUM159 lineages make this comparison interpretable.
3. Rewrite the second paragraph to say these observations motivated the model.
4. Keep model concepts:
   - resource stress can oppose high ploidy through growth/death costs;
   - stress can increase CIN/WGD;
   - high ploidy can buffer otherwise lethal missegregations;
   - outcome depends on the balance among growth suppression, death, CIN/WGD generation, and post-missegregation survival.
5. After resolving model-figure duplication, update all citations in this paragraph to the single retained model-overview label.

Validation:

```bash
rg -n "SUM-159|SUM159 lineages motivate|whole-genome doubling|WGD|ploidy reduction|resource-stress model" ltee_hypoxia_model.tex
rg -n "haploidy|hypoploidy|sen generation|post-missing|Figure 3B|Figure 3C|Figure 3D" ltee_hypoxia_model.tex
```

Expected:

- SUM159 may appear in body text, but not in the subsection title.
- Opening Results paragraph contains ploidy reduction and WGD/high-ploidy expansion.
- No speech-to-text artifacts remain.

### Step 3: Clean Figure 1 Caption

Targets:

- Figure 1 caption lines 125-126.

Actions:

1. Keep the title unless later comments reject it.
2. Remove workflow/provenance language:
   - `Existing optimization-data-stream overview`;
   - `This panel was copied...`;
   - any mention of automatically generated figures.
3. Remove the final interpretive sentence:
   - `the same isogenic 2N/4N system underlies both settings...`.
4. Keep concise content:
   - matched near-diploid and near-tetraploid SUM159 lineages;
   - in vitro control and hypoxia/anoxia passage histories;
   - target oxygen encoded by blue intensity;
   - flow cytometry and karyotyping events;
   - in vivo untreated harvest tumor-burden trajectories;
   - terminal chromosome-number and necrosis measurements.

Validation:

```bash
rg -n "Existing optimization|copied|automatically generated|used here because|Panel file" ltee_hypoxia_model.tex
```

### Step 4: Resolve The Duplicate Model-Overview Figure

Targets:

- First model figure block lines 130-142, label `fig:iteration1-mechanismDiagram`.
- Second model figure block lines 145-154, label `fig:iteration1-resource-limitation-mechanism-missing`.
- Results citation line 91.

Actions:

1. Decide which model-overview block to retain by checking which label should be used in the Results text.
2. Remove the duplicate block entirely.
3. Retained model-overview caption should be framed as assumptions/mechanism overview, not proof.
   - Candidate title: `Model overview: resource limitation can impose opposing pressures on ploidy evolution.`
4. Keep intended manuscript panels visible even if their artwork is not available yet. Use concise "panel to be generated" language for missing artwork rather than hiding the panel from the compiled manuscript.
5. Update line 91 citations to the retained label.

Validation:

```bash
rg -n "iteration1-mechanismDiagram|iteration1-resource-limitation-mechanism-missing" ltee_hypoxia_model.tex
rg -n "Resource limitation creates opposing pressures|Model overview" ltee_hypoxia_model.tex
```

Expected:

- Only one model-overview figure block remains.
- Title includes `Model overview` or equivalent.
- Opening Results cites the retained model-overview label.

### Step 5: Rewrite In Vitro Results Mechanism

Targets:

- Subsection title line 93.
- Paragraphs lines 95 and 97.

Actions:

1. Keep or lightly polish the title. Current high-ploidy wording is already better than the transcript's hypoploid wording.
2. In line 95, explicitly mention that the separate in vitro fit recapitulates transient WGD/high-ploidy expansion before later lower-ploidy reshaping.
3. In line 97, replace the current double-edged sentence.
   - Avoid: "it preserves high-ploidy cells after missegregation".
   - Use: high-ploidy cells missegregate more often because they carry more chromosomes, and high ploidy also makes some chromosome-loss daughters viable; together this gives a route from high ploidy back toward lower chromosome number.
4. Keep the rejected/secondary interpretation that the fit does not require strong direct high-ploidy elimination, and keep the intended negative-control panel present as a panel-to-be-generated placeholder if the real evidence panel does not exist yet.

Suggested replacement concept:

```text
The double-edged feature is that high-ploidy cells generate more chromosome-loss opportunities, while the same ploidy-dependent buffering makes a subset of those chromosome-loss daughters viable. This gives the model a route from high-ploidy expansion back toward lower chromosome number without requiring strong direct killing of high-ploidy cells.
```

Validation:

```bash
rg -n "preserves high-ploidy cells after missegregation|hypoploidy|hypoploid|protein vitrofit|WGD|whole-genome" ltee_hypoxia_model.tex
```

### Step 6: Restore In Vitro Panels B, D, And WGD Recapitulation

Targets:

- Assembled figure include line 158.
- In vitro caption lines 162-168.
- Iteration-1 panel files under `oxygen/figures/iteration1`.
- Assembled output `oxygen/figures/assembled_fig3.png`.

Actions:

1. Restore or regenerate the prior iteration-1 missegregation-source panel if the artwork exists or can be rebuilt now; otherwise keep Fig. 3B present as a panel-to-be-generated placeholder.
   - Current line 164 must not say Fig. 3B is missing; use "panel to be generated" if the artwork is unavailable.
   - Panel should show high-ploidy cells missegregate more often or have more missegregation opportunities.
2. Restore or regenerate the prior downward-reshaping panel if the artwork exists or can be rebuilt now; otherwise keep Fig. 3D present as a panel-to-be-generated placeholder.
   - Current line 166 must not say Fig. 3D is missing; use "panel to be generated" if the artwork is unavailable.
   - If the available panel is a distribution/trajectory heatmap, caption it as trajectory-level evidence for WGD/high-ploidy expansion followed by lower-ploidy reshaping, not direct single-cell lineage tracing.
3. Restore or regenerate the prior WGD-recapitulation panel.
   - The prior selection log named `fig3b_invitro_predicted_ploidy_distribution.png` as a removed in vitro distribution panel; check whether that is the WGD-recapitulation panel.
   - If yes, restore it and assign it to the most appropriate in vitro panel letter.
4. Keep panel E in the main manuscript as an intended negative-control/rejected-interpretation panel.
   - If no direct evidence panel exists yet, label it as a panel to be generated rather than removing it from the caption.
5. Rerun the assembler:

```bash
python3 oxygen/figures/assemble_iteration_panels.py oxygen/figures/iteration1
```

Validation:

```bash
find oxygen/figures/iteration1 -maxdepth 1 -type f | sort
file oxygen/figures/assembled_fig3.png
rg -n "Fig\\. 3B: missing|Fig\\. 3D: missing|HTML report|Panel file|extracted|PNG|WGD|whole-genome|recapitulat" ltee_hypoxia_model.tex oxygen/figures/iteration1/figureCaptions.txt
```

Expected:

- In vitro assembled figure includes the restored mechanism/trajectory panels.
- Main text and caption explicitly mention WGD recapitulation.
- Captions do not contain report/file provenance.

### Step 7: Clean In Vitro Caption

Targets:

- Lines 162-168.

Actions:

1. Remove:
   - `Panel file`;
   - `extracted from`;
   - HTML report names;
   - "current iteration-1";
   - "removed panel" wording.
2. Keep:
   - panel A: fit tracks growth, chromosome counts, and live/dead burden across 2N/4N control and oxygen-deprived branches;
   - panel B: missegregation source/opportunity;
   - panel C: post-missegregation viability as a function of chromosome number;
   - panel D: WGD/downward reshaping trajectory evidence;
   - panel E as the rejected strong-elimination interpretation; if no real panel exists yet, keep it as a panel-to-be-generated placeholder.
3. Verify panel C curve orientation against the actual plot before saying increasing/decreasing.

Validation:

```bash
rg -n "Panel file|extracted|HTML report|removed|iteration-1|current iteration" ltee_hypoxia_model.tex
```

### Step 8: Revise In Vivo Fixed-O2 And Landscape Text

Targets:

- In vivo section title line 99.
- Oxygen/resource paragraph line 101.
- Fixed-O2 paragraph line 103.
- Parameter-landscape paragraph line 105.
- Response-class paragraph line 107.

Actions:

1. Keep line 101 concept:
   - oxygen is the formal model variable and a first-order proxy for broader effective resource stress;
   - fitted effective O2 trajectories are latent model-implied variables, not direct measurements;
   - in vivo hypoxia is not identical to low in vitro oxygen in rich medium.
2. Before interpreting fixed-O2 feature-AUC results, add a concise explanation:
   - what was evaluated at fixed O2;
   - what "dominant ploidy mode" means;
   - that fitted parameters or derived features were scored for how well they separated lower- versus higher-ploidy outcomes;
   - why AUC was used to identify compact discriminating axes.
   - Defer exact feature counts, thresholds, and scoring details to Methods/Supplement unless a short parenthetical is essential.
3. Keep low/intermediate/high O2 interpretation after this explanation.
4. Before interpreting landscape clusters, add a concise explanation:
   - what objects were clustered, e.g. fitted in vivo parameter sets/seeds;
   - that scale/prior-aware preprocessing was used so parameter units and bounds did not drive the comparison;
   - why clustering was useful for the biological question.
   - Defer exact preprocessing, dimensionality-reduction, and clustering settings to Methods/Supplement.
5. Clean line 105:
   - remove `rather than purely optimizer noise`;
   - remove `rather than only to embedding coordinates`;
   - keep biological interpretation through baseline missegregation and oxygen-response shape.
6. Move the prior-aware preprocessing statement upward into the analysis-description sentence, not as an afterthought after the biological interpretation.
7. Delete line 107 paragraph from the main Results section.
   - The transcript explicitly says to leave this paragraph out of the in vivo section.
   - Do not discard the concept entirely; reintroduce fixed-O2 dominant-ploidy response curves in the final Results subsection as the model-selection discriminator after the joint-fit axes.

Validation:

```bash
rg -n "feature-AUC|feature AUC|dominant ploidy mode|mode identity|cluster|prior|bound|optimizer noise|embedding coordinates|model-selection view" ltee_hypoxia_model.tex
```

Expected:

- Fixed-O2 and landscape analyses are introduced before interpretation.
- "Mode" and "feature-AUC" are defined.
- No optimizer-noise or embedding-coordinate contrast language remains.
- The response-class paragraph is absent from the in vivo subsection.
- Fixed-O2 response-curve material is reserved for the final joint/model-selection subsection.

### Step 9: Revise Final Joint-Fit Results And Reintroduce Fixed-O2 Response Curves

Targets:

- Joint-fit subsection title line 109.
- Joint-fit interpretation paragraph line 111.
- Joint-fit axes paragraph line 113.
- Prior in vivo response-curve paragraph line 107, to be moved/reframed.
- Fixed-O2 response-class figure lines 206-217, if retained.

Actions:

1. Keep the joint-fit subsection as the final Results subsection unless later manuscript structure changes.
2. In line 111, keep the main interpretation:
   - in vivo and in vitro differ less in whether oxygen-linked stress exists;
   - they differ more in how stress maps to proliferation, CIN, and post-missegregation survival.
3. Add the caveat requested in the transcript:
   - in vivo stress should be framed as oxygen-linked stress plus additional tumor-environment stressors, not oxygen alone.
4. Keep the three context-specific joint-fit axes:
   - lower effective proliferative ceiling in vivo;
   - stronger stress-linked missegregation/CIN generation in vivo;
   - more ploidy-dependent post-missegregation survival filtering in vivo.
5. Preserve corrected terminology:
   - use `CIN`, not transcript artifact `syn`;
   - use `post-missegregation survival`, not `post-mitosis survival`;
   - use `parameter-ratio summary`, not `progeny ratio summary`, if the ratio panel is mentioned.
6. Remove the final sentence of line 113:
   - `The integrated parameter-ratio summary therefore supports...`
7. Replace that sentence with the fixed-O2 response-curve ending:
   - After the joint/in vivo analyses, multiple plausible parameter regimes remain that could explain the in vivo ploidy evolution.
   - Use fixed-O2 dominant-ploidy response curves, computed from the fixed-O2 eigenvector or steady-state solution, to discriminate among these regimes.
   - The biologically expected behavior is that higher oxygen should correlate with higher dominant ploidy.
   - Only a subset of fitted solutions shows this monotone increasing O2-ploidy relationship, so those solutions should be prioritized as more biologically plausible.
   - Explain which parameter regimes/features predict the monotone increasing fixed-O2 ploidy response.
8. Before writing `eigenvector`, verify the exact computational object in the code/report.
   - If the code computes a dominant eigenvector, use that term.
   - If it computes a steady state, attractor, or long-time fixed-O2 composition, use the exact term.
9. Reuse the Figure 6 response-class material if it supports this final argument.
   - Remove provenance language from the Figure 6 caption.
   - Consider retitling the figure to emphasize fixed-O2 response curves as model-selection evidence, not generic "joint-fit behaviors".

Validation:

```bash
rg -n "oxygen stress itself|other stressors|CIN generation|syn|post-mitosis|post-missegregation|integrated parameter-ratio summary|fixed-O\\$_2\\$|dominant-ploidy response|monotone" ltee_hypoxia_model.tex
rg -n "eigenvector|steady state|attractor|fixed-O2|fixed_O2|dominant ploidy" oxygen/code oxygen/figures docs -g '*.{R,py,md,txt,tex}'
```

Expected:

- The Results ending explains why fixed-O2 response curves are used after joint fitting.
- Monotone increasing O2-ploidy behavior is framed as a biologically plausible selection criterion.
- The text identifies or promises to identify the parameter regimes that predict that behavior.
- Transcript artifacts `syn` and `post-mitosis survival` do not appear.

### Step 10: Clean In Vivo, Joint, Response-Curve, And Remaining Captions For Provenance

Targets:

- In vivo caption lines 180-185.
- Joint caption lines 197-202.
- Fixed-O2 response-class caption lines 211-215.
- Comments lines 117-119 and panel-file comments throughout, if they leak into manuscript-facing text.

Actions:

1. Remove manuscript-facing file/report provenance from all captions.
2. Keep scientific descriptions of what each panel shows.
3. If Figure 6 is retained as the response-curve/model-selection figure, revise its caption to support the final Results ending:
   - fixed-O2 response classes across candidate fits;
   - monotone increasing O2-ploidy behavior as the expected/prioritized class;
   - parameter regimes/features associated with that behavior if the evidence panel exists.

Validation:

```bash
rg -n "Panel file|HTML report|report panel|extracted|copied|automatically generated|current iteration|first-pass" ltee_hypoxia_model.tex
```

### Step 11: Compile/Static Validate

Run after edits:

```bash
rg -n "TODO|XX|Panel file|HTML report|extracted|copied|automatically generated|optimizer noise|embedding coordinates|hypoploid|hypoploidy|haploidy|sen generation|post-missing|syn|post-mitosis" ltee_hypoxia_model.tex
file oxygen/figures/assembled_fig*.png
```

If LaTeX tooling is available:

```bash
latexmk -pdf ltee_hypoxia_model.tex
```

If `latexmk` is unavailable, record that explicitly and rely on static checks.

## Open Decisions Before Editing

1. Which duplicate model-overview block should be retained?
   - Preferred: keep one model-overview figure label and update all references to it.
2. Which prior iteration-1 files correspond to:
   - missegregation source;
   - downward reshaping;
   - WGD recapitulation?
3. Should the rejected-interpretation panel E stay in the in vitro figure?
   - Yes. Keep it as a panel-to-be-generated placeholder until direct evidence/artwork exists.
4. What is the exact computational term for the fixed-O2 solution used in the final Results ending?
   - Use `dominant eigenvector` only if that is what the code/report computes.
   - Otherwise use `steady state`, `attractor`, or `long-time fixed-O2 composition`.
5. Which panel or table identifies the parameters/regimes that predict monotone increasing fixed-O2 ploidy?
   - If no panel exists yet, the main text should state the intended analysis cautiously and the manuscript caption should keep the intended panel visible as a panel-to-be-generated placeholder.
