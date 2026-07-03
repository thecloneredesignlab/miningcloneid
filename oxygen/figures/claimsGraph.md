 panel-first claim graph: every node is allowed into the graph only if we can state what visual panel would represent it. Evidence sufficiency comes later; for now each node gets a proposed panel contract.

  Graph Structure
  Use four layers:

  1. System Layer
     What was measured and why the comparison is interpretable.

  2. Mechanism Layer
     How the model says oxygen/resource stress can change ploidy evolution.

  3. Inference Layer
     What the separate and joint fits infer about in vitro, in vivo, and context differences.

  4. Landscape/Regime Layer
     How fitted solutions behave across oxygen regimes and parameter-space regions.

  Edges should mean one of three things only:

  - supports: child provides evidence for parent.
  - mechanism: child explains how parent could occur.
  - qualifies: child limits or contextualizes parent.

  Main Figure Logic
  I would map the main claims into figures like this:

   Node    Claim                                                                                Panel concept
  ━━━━━━  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
   R0      Resource limitation rewires ploidy evolution through opposing costs and buffering    Summary schematic: resource stress increases costs but also increases CIN/WGD, with high ploidy buffering missegregation
  ──────  ───────────────────────────────────────────────────────────────────────────────────  ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
   C1      Matched SUM159 2N/4N lineages enable in vitro/in vivo comparison                     Timeline panel showing same 2N/4N starting system, split into in vitro and in vivo measurements
  ──────  ───────────────────────────────────────────────────────────────────────────────────  ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
   C2      Resource limitation can oppose and promote high ploidy                               Mechanism panel: high chromosome content cost vs buffering benefit, with arrows to growth/death/CIN/survival
  ──────  ───────────────────────────────────────────────────────────────────────────────────  ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
   C3      Tumor setting is broader than lower oxygen                                           In vivo vs in vitro stress-context panel: fixed known O2 in vitro versus effective resource stress in vivo
  ──────  ───────────────────────────────────────────────────────────────────────────────────  ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
   C4      Joint fits infer context-specific proliferation/CIN/survival axes                    Paired parameter/ratio forest plot: in vivo vs in vitro soft-coupled axes
  ──────  ───────────────────────────────────────────────────────────────────────────────────  ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
   C5      Dominant ploidy is resource-regime dependent                                         Fixed-O2 regime panel: low/intermediate/high O2 with dominant drivers annotated
  ──────  ───────────────────────────────────────────────────────────────────────────────────  ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
   C6      Multiple fitted solution regions may represent distinct biological regimes           Prior-aware landscape plot with clusters/regions and representative curve classes
  ──────  ───────────────────────────────────────────────────────────────────────────────────  ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
   C7      In vitro fits show oxygen reshapes ploidy without strong high-ploidy elimination     In vitro fit panel: growth plus ploidy trajectories, with weak growth-penalty parameters annotated
  ──────  ───────────────────────────────────────────────────────────────────────────────────  ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
   C8      In vitro ploidy reduction can arise from CIN plus survival buffering                 Mechanistic trajectory panel: high-ploidy cell missegregates, viable chromosome-loss daughters move downward

  That gives the manuscript a natural figure sequence:

  - Figure 1: Experimental system and data streams.
  - Figure 2: Model mechanism.
  - Figure 3: Separate in vitro and in vivo fits.
  - Figure 4: Joint in vivo/in vitro context differences.
  - Figure 5: Fixed-O2 resource-regime behavior.
  - Figure 6: Parameter landscape and candidate biological regimes.

  Supporting Nodes
  The secondary claims should mostly become subpanels or supplement panels:

   Node    Claim                                                                      Panel concept
  ━━━━━━  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
   S1      Mode predictability depends on O2                                          AUC by O2 level, top features annotated
  ──────  ─────────────────────────────────────────────────────────────────────────  ────────────────────────────────────────────────────────────────────────────────
   S2      Low-resource behavior is death dominated                                   Feature-importance or coefficient panel at low O2
  ──────  ─────────────────────────────────────────────────────────────────────────  ────────────────────────────────────────────────────────────────────────────────
   S3      Intermediate O2 is tradeoff-dependent                                      Multi-feature/interaction panel at intermediate O2
  ──────  ─────────────────────────────────────────────────────────────────────────  ────────────────────────────────────────────────────────────────────────────────
   S4      High-resource behavior depends on baseline missegregation and buffering    High-O2 feature panel, emphasizing p_mis_base and buffering terms
  ──────  ─────────────────────────────────────────────────────────────────────────  ────────────────────────────────────────────────────────────────────────────────
   S5      Continuous dominant-ploidy supports discrete mode analysis                 Scatter/side-by-side: continuous dominant ploidy vs binary mode classification
  ──────  ─────────────────────────────────────────────────────────────────────────  ────────────────────────────────────────────────────────────────────────────────
   S6      Full fixed-O2 curve classes are richer than binary mode                    Curve-class panel across O2 grid
  ──────  ─────────────────────────────────────────────────────────────────────────  ────────────────────────────────────────────────────────────────────────────────
   S7      Oxygen is incomplete proxy for resource limitation                         Conceptual limitation panel or schematic, likely supplement/discussion
  ──────  ─────────────────────────────────────────────────────────────────────────  ────────────────────────────────────────────────────────────────────────────────
   S8      CIN therapies may have ploidy-dependent effects                            Prediction schematic, likely discussion/supplement unless data included
  ──────  ─────────────────────────────────────────────────────────────────────────  ────────────────────────────────────────────────────────────────────────────────
   S9      Current in vitro fit does not prove ploidy-independent division            Ablation-design panel: current model vs ploidy-neutral division model

  Implementation Guardrail Nodes
  These should not sit in the biological graph as evidence. I would put them in a separate “Methods Validity” side graph:

  - fixed-O2 attractor versus finite-time trajectories;
  - prior-aware preprocessing requirement;
  - model-implied biological features;
  - run provenance / soft-coupling settings / warm-start initialization.

  Each can still get a panel, but mostly as supplement or methods panels.

  Practical Format
  I would create a claim-graph manifest with columns:

  node_id, parent_id, edge_type, claim_text, node_type, figure_panel_goal, visual_form, primary_or_supplement, later_evidence_status.

  The key discipline is that every node gets a visual form immediately. If a node cannot be visualized even conceptually, it is probably not a graph node; it is prose, framing, or a limitation.

## From Panel-First Claim Graph To First Graph Version

The first claim graph should be generated now, before the figure panels exist.

At this stage the graph is not an evidence audit and not a final manuscript figure. It is a figure-planning scaffold: it freezes the claims as nodes, assigns each node a panel contract, and makes explicit which panels need to exist for the manuscript logic to work. Figure panels then become the way each node is tested, refined, merged, moved to supplement, or dropped.

The workflow should therefore be:

1. Freeze the claim inventory for version 0.

   Use the current accepted claims as the source of truth. Do not yet decide whether the evidence is sufficient. The purpose of this pass is only to translate claims into a graph and panel plan.

2. Create a claim graph manifest.

   Make a table where each row is one graph node. Minimum columns:

   - `node_id`
   - `parent_id`
   - `edge_type`
   - `claim_text`
   - `node_type`
   - `figure_number`
   - `panel_id`
   - `panel_goal`
   - `visual_form`
   - `primary_or_supplement`
   - `evidence_status`

   For version 0, `evidence_status` should usually be `not_audited`. Evidence sufficiency is deliberately postponed.

3. Generate the version-0 claim graph from the manifest.

   The first generated graph should be a simple directed graph with the top-level thesis at the root, the main claims underneath, and the secondary or methods nodes attached where they clarify or support the main claims.

   This graph should answer:

   - What does the manuscript need to convince the reader of?
   - Which claims depend on which other claims?
   - Which figure or panel is responsible for each node?
   - Which nodes are central versus supporting or interpretive?

## Claim Graph Visualization

The version-0 claim graph should be visualized as a left-to-right directed graph.

The visual hierarchy should be:

1. **Root thesis**

   A single root node at the far left: resource limitation reshapes ploidy evolution through opposing chromosome-content costs and buffering advantages.

2. **Main claim nodes**

   Main claims should form the central spine of the graph. They should be grouped by the natural figure sequence rather than by the order in which analyses were performed:

   - Figure 1: experimental system and data streams.
   - Figure 2: model mechanism.
   - Figure 3: separate in vitro and in vivo fits.
   - Figure 4: joint in vivo/in vitro context differences.
   - Figure 5: fixed-O2 resource-regime behavior.
   - Figure 6: parameter landscape and candidate biological regimes.

3. **Supporting nodes**

   Secondary claims should sit to the right of, or underneath, the main claim they support. They should not compete visually with the main spine.

4. **Methods and interpretation guardrails**

   Guardrail nodes should be shown as a separate side band. They qualify interpretation but are not biological evidence nodes.

Each node should display:

- short node ID, such as `R0`, `C4`, or `S2`;
- one-line claim label;
- assigned figure/panel, such as `Fig. 4B`;
- node type, encoded by color or shape.

Recommended node types:

- root thesis: dark filled node;
- main biological claim: strong color;
- supporting analysis claim: lighter color;
- methods or interpretation guardrail: neutral gray;
- implementation/provenance note: dashed or muted node, preferably outside the main biological graph.

Edges should be visually encoded by edge type:

- `supports`: solid arrow;
- `mechanism`: solid arrow with a distinct color or label;
- `qualifies`: dashed arrow.

The graph should avoid dense cross-links. If a node appears to support many claims, prefer either:

- making it a shared upstream mechanism node, or
- assigning it to the strongest parent and listing secondary relationships in the manifest.

For version 0, the graph can be rendered as Mermaid inside `oxygen/figures/claim_graph_v0.md`. If the graph becomes too dense for Mermaid, switch to Graphviz/DOT and export both SVG and PNG. The visual graph is a navigation tool, not the final biological figure; final manuscript figures are generated from the panel backlog.

4. Generate a panel backlog from the same manifest.

   The panel backlog is the practical bridge from graph to figures. It should group nodes by the natural figure sequence:

   - Figure 1: Experimental system and data streams.
   - Figure 2: Model mechanism.
   - Figure 3: Separate in vitro and in vivo fits.
   - Figure 4: Joint in vivo/in vitro context differences.
   - Figure 5: Fixed-O2 resource-regime behavior.
   - Figure 6: Parameter landscape and candidate biological regimes.

   Each panel backlog entry should state the claim, the intended visual form, and what a reader should learn from the panel in one sentence.

5. Build first-pass panel sketches.

   These can be schematic or placeholder panels. The goal is not polish. The goal is to test whether the graph maps cleanly onto the figure sequence and whether any figure is overloaded or missing a needed claim.

6. Build data-backed panels.

   Once the panel sketches are coherent, replace placeholders with data-backed plots and model diagrams. This is where evidence sufficiency starts to matter, but it should be recorded as panel/node status rather than handled ad hoc.

7. Audit evidence after the first panel pass.

   Only after the graph and first panel draft exist should each node be assigned an evidence status such as:

   - `supported_main`
   - `supported_supplement`
   - `provisional`
   - `needs_analysis`
   - `drop_or_merge`

   This audit will likely move some nodes from main figures to supplement, merge overlapping claims, or turn some biological claims into interpretation guardrails.

## Version-0 Deliverables

The immediate deliverables should be:

1. `oxygen/figures/claim_graph_manifest.csv`

   A machine-readable claim graph and panel contract table.

2. `oxygen/figures/claim_graph_v0.md`

   A readable graph description, ideally including a Mermaid or Graphviz diagram generated from the manifest.

3. `oxygen/figures/figure_panel_backlog.md`

   A figure-by-figure checklist of panels required by the claim graph.

The claim graph should be regenerated whenever the manifest changes. The figure panels should be generated after the version-0 graph exists, because the graph defines what each panel is supposed to prove or clarify.
