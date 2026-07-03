# Figure Panel Backlog V0

This backlog is generated from `oxygen/figures/claim_graph_manifest.csv`. It is panel-first and evidence-neutral: every entry states the claim, the intended visual form, and the one-sentence reader takeaway the panel should deliver.

## Figure 1: Experimental System And Data Streams

| Panel | Node | Claim | Intended visual form | Reader should learn |
|---|---|---|---|---|
| Fig. 1A | C1 | Matched SUM159 near-2N and near-4N lineages enable environment-dependent in vitro/in vivo comparison. | Timeline/data-stream panel showing the same 2N/4N starting system split into in vitro and in vivo measurements. | The culture and tumor comparisons come from the same lineage system, so differences can be framed as environment dependent. |

## Figure 2: Model Mechanism

| Panel | Node | Claim | Intended visual form | Reader should learn |
|---|---|---|---|---|
| Fig. 2A | R0 | Resource limitation rewires ploidy evolution through opposing chromosome-content costs and buffering advantages. | Summary schematic linking resource stress to chromosome-content cost, CIN/WGD pressure, and high-ploidy missegregation buffering. | The manuscript thesis is a balance between resource-linked costs and high-ploidy buffering advantages. |
| Fig. 2B | C2 | Resource limitation can both oppose and promote high ploidy. | Mechanism schematic showing chromosome-content cost and stress death opposing high ploidy while CIN/WGD and survival buffering preserve or create high-ploidy states. | The same resource-stress axis can push ploidy upward or downward depending on the balance of model mechanisms. |

## Figure 3: Separate In Vitro And In Vivo Fits

| Panel | Node | Claim | Intended visual form | Reader should learn |
|---|---|---|---|---|
| Fig. 3A | C7 | Separate in vitro fits show oxygen reshapes ploidy without requiring strong high-ploidy elimination. | In vitro fit panel showing growth and ploidy trajectories with weak lower-bound high-ploidy growth-penalty parameters annotated. | Culture can remain permissive for high ploidy while oxygen stress still reshapes ploidy trajectories. |
| Fig. 3B | C8 | In vitro ploidy reduction can arise from CIN-generated variation plus ploidy-dependent missegregation survival buffering. | Mechanistic trajectory panel showing high-ploidy cells missegregating more often, viable chromosome-loss daughters, and downward ploidy movement. | Downward ploidy reshaping can come from the CIN kernel and survival filtering without making direct high-ploidy elimination the central explanation. |
| Fig. 3C | S10 | Separate in vivo fits provide the in vivo side of the context comparison by testing tumor burden and terminal ploidy/necrosis constraints. | Separate in vivo fit panel comparing observed and fitted tumor burden trajectories with terminal ploidy and necrosis constraint summaries, keyed to the in_vivo_top_10 report set. | The reader sees the separate in vivo fit evidence plan before the joint model is used to compare context-specific axes. |

## Figure 4: Joint In Vivo/In Vitro Context Differences

| Panel | Node | Claim | Intended visual form | Reader should learn |
|---|---|---|---|---|
| Fig. 4A | C3 | The tumor setting is broader than lower oxygen and should be interpreted as an effective in vivo resource-stress regime. | In vivo versus in vitro stress-context schematic contrasting fixed known O2 in culture with effective tumor resource stress. | Oxygen is the formal model variable, but the in vivo interpretation is a broader tumor resource-stress axis. |
| Fig. 4B | C4 | Joint soft-coupled fits infer context-specific proliferation, CIN, and survival-filter axes between in vivo and in vitro. | Paired parameter or ratio forest plot comparing in vivo and in vitro soft-coupled axes. | The joint model should show which proliferation, missegregation, and survival-filter axes differ across contexts. |

## Figure 5: Fixed-O2 Resource-Regime Behavior

| Panel | Node | Claim | Intended visual form | Reader should learn |
|---|---|---|---|---|
| Fig. 5A | C5 | Dominant ploidy state is resource-regime dependent. | Fixed-O2 regime panel with low/intermediate/high O2 bands and dominant drivers annotated. | Ploidy outcome depends on the resource regime rather than one universal driver. |
| Fig. 5B | S1 | Steady-state ploidy mode is predictable from fitted parameters, but predictability depends on O2/resource level. | AUC-by-O2 panel with top single parameters or small feature sets annotated. | Mode predictability should be read as discrimination across O2 levels, not as variance explained. |
| Fig. 5C | S2 | Low-resource ploidy behavior is mainly death dominated. | Low-O2 feature-importance or coefficient panel emphasizing stress-associated death scale and death-reduction benefit. | At very low O2, death terms provide the compact explanation for lower versus higher ploidy mode behavior. |
| Fig. 5D | S3 | Intermediate-resource behavior is the least reducible and most tradeoff dependent. | Intermediate-O2 interaction or multi-feature panel covering death selectivity, missegregation sensitivity, buffering, thresholding, and baseline missegregation. | Intermediate O2 is where outcome boundaries depend most on multi-parameter tradeoffs. |
| Fig. 5E | S4 | High-resource ploidy behavior is mainly governed by baseline missegregation, with buffering modifying the outcome. | High-O2 feature panel emphasizing baseline missegregation and buffering strength. | When death pressure is reduced, ongoing missegregation and buffering control which ploidy state persists. |
| Fig. 5F | S5 | Continuous dominant-ploidy analyses support the discrete mode analyses. | Side-by-side scatter or calibration plot comparing continuous dominant ploidy with Mode 1/Mode 2 classification. | The discrete mode analysis should agree with continuous dominant-ploidy behavior before it is treated as robust. |

## Figure 6: Parameter Landscape And Candidate Biological Regimes

| Panel | Node | Claim | Intended visual form | Reader should learn |
|---|---|---|---|---|
| Fig. 6A | C6 | Multiple fitted solution regions may represent distinct biological regimes rather than purely optimizer noise. | Prior-aware parameter-landscape plot with clusters or regions, representative fixed-O2 curve classes, and candidate regime labels. | Structured fitted-parameter regions can be interpreted as candidate biological regimes if preprocessing and provenance are controlled. |
| Fig. 6B | S6 | Full fixed-O2 ploidy-response curve classes provide a richer phenotype than binary mode labels. | Curve-class panel across the fixed-O2 grid showing nonmonotone, monotone, U-shaped, inverted-U-shaped, flat, and transition-plateau classes. | Curve classes preserve response-shape information that binary mode labels collapse. |

## Supplemental And Methods Panels

| Panel | Node | Claim | Intended visual form | Reader should learn |
|---|---|---|---|---|
| Fig. S2A | S8 | CIN-inducing therapies may have ploidy-dependent effects under the survival-filter framework. | Prediction schematic contrasting near-diploid sensitivity to missegregation dosage imbalance with high-ploidy buffering of chromosome gains or losses. | The therapy point is an implication of the mechanism, not a primary result unless direct drug-response analyses are added. |
| Fig. S2B | P3 | The manuscript description of dead compartments should match the implementation. | Model-routing schematic showing boundary-dropped ordinary missegregation daughters and out-of-grid WGD events entering the CIN-associated or boundary/nonviability-associated dead compartment. | The mechanism text should name the nonviable compartment consistently with the implemented routing. |
| Fig. S3A | S9 | The current in vitro fit does not prove that division rate is ploidy independent. | Ablation-design panel comparing the current weak-penalty model with a ploidy-neutral division-rate refit or fully neutral alpha_o2/gamma_growth variant. | A separate ablation is needed before claiming true ploidy-independent division. |
| Fig. S4A | S7 | Oxygen is an informative but incomplete proxy for resource limitation. | Conceptual limitation schematic showing oxygen as the formal variable within broader tumor resource constraints. | The model formalizes oxygen-linked stress, while the biology likely includes nutrient, perfusion, waste, density, stromal, and tissue-structure constraints. |
| Fig. S4B | P1 | Default soft-coupling settings and parameter lists affect comparability across joint-fit runs. | Provenance table or schematic listing soft-coupling settings, parameter lists, and branch/config compatibility. | Joint-fit comparisons need config and parameter-list provenance before biological differences are interpreted. |
| Fig. S4C | P4 | The manuscript-described analysis should match the default implementation, including warm-start initialization when used. | Run-provenance panel indicating warm-start initialization, analyzed config, and whether the manuscript workflow is default or explicitly specified. | The manuscript should either use the default workflow or explicitly state how analyzed runs were initialized. |
| Fig. S5A | M1 | Fixed-O2 attractor analyses and finite-time experimental trajectories are related but not identical. | Methods schematic contrasting dominant-eigenvector fixed-O2 attractor behavior with finite-time experimental trajectory constraints. | Long-time fixed-O2 attractors should qualify, not replace, finite-time fit interpretation. |
| Fig. S6A | M2 | Parameter-landscape analyses need prior-aware preprocessing before PCA, UMAP, t-SNE, or clustering are interpreted biologically. | Workflow panel showing transformed/scaled prior-aware vectors feeding embeddings and clustering, then characterization in original or derived feature space. | Landscape interpretation depends on preprocessing that respects transforms, scales, priors, and bounds. |
| Fig. S6B | M3 | Model-implied biological features help translate fitted parameters into interpretable regimes. | Derived-feature overlay panel for net growth, turnover, death/proliferation balance, O2 thresholds, missegregation burden, post-missegregation viability, and effective oxygen supply regime. | Biological regime labels should be supported by interpretable derived features, not embedding position alone. |
| Fig. S6C | P2 | Start tables and bound expansion can silently change the optimizer region. | Traceability panel showing start-table source, bound expansion status, and resulting feasible parameter region. | Parameter-landscape comparisons need optimizer-region provenance because starts and bounds can change the search space. |
