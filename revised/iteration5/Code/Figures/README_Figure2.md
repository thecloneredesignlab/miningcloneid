# Figure 2: Stress Feedback and Progressive Presentation Build

## Review

- [Full Figure 2](../../Figures/assembled_fig2.png) and [vector PDF](../../Figures/assembled_fig2.pdf).
- [Animation](../../Figures/figure2_animation/figure2_progressive_build.mp4): 59.6 seconds, silent, 1920 x 1080, H.264.
- [Presenter-paced PowerPoint](../../Figures/figure2_animation/figure2_progressive_build.pptx): seven click-to-advance stages with speaker notes and fade transitions.
- [Storyboard PDF](../../Figures/figure2_animation/figure2_progressive_build.pdf) and [speaker notes](../../Figures/figure2_animation/speaker_notes.md).

The PowerPoint uses source-rendered slide images, not editable individual objects. The MP4 includes the moving-node and drawn-arrow transitions. Temporary video frames and review caches are not versioned.

## Figure Changes

Panel A incorporates the manually approved diagram edits: karyotype variation feeds selection and retention into adapted composition; a separate Post-MS survival node and the outer feedback path are removed; activation arrows and inhibitory T-bars distinguish the signed connections; and yellow CIN-associated daughter loss returns to the aggregate Cell death node. Panels B and C are unchanged.

The model equations, fitted parameters, simulations, and other manuscript figures are unchanged. This is a schematic presentation revision, not a new model or analysis.

## Animation Order

1. Recreate the slide-6 resource limitation / CIN / ploidy triangle in manuscript colors.
2. Resolve the resource-to-CIN edge through experienced resource-stress death hazard.
3. Split ploidy into karyotype variation and adapted population composition.
4. Close adapted composition -| death hazard. Hold for 12 seconds on the central takeaway: adaptation can reduce average experienced stress and CIN even at unchanged oxygen.
5. Add proliferation and its resource/composition dependencies.
6. Add CIN-associated daughter loss and broaden the aggregate label to Cell death.
7. Add WGD generation at constant probability per division, ending with the exact current `draw_panel_a()` function.

The simple resource-cost edge is temporarily set aside while the CIN branch is unpacked, then resolved through explicit growth and death effects. Peripheral slide annotations are omitted. WGD is not folded into the stress-induced-CIN shorthand.

## Scientific Interpretation

Only the resource-stress hazard drives inducible per-chromosome missegregation in the equations. CIN-associated nonviable-daughter production is a distinct division-linked loss flux; the aggregate Cell death node does not introduce an extra algebraic positive feedback into the CIN equation. Adaptation can lower population-average hazard and CIN at constant oxygen, but is not guaranteed in every fitted trajectory. Ploidy-dependent survival filtering remains part of selection and retention. WGD probability is constant per division, not directly oxygen dependent.

## Rebuild Without Fitting Data

Run from the repository root. The worker flag bypasses the full figure pipeline and its data requirements for this self-contained schematic:

```bash
FIGURE2_DRAW_WORKER=1 Rscript --vanilla revised/iteration5/Code/Figures/draw_Figure2.R
Rscript --vanilla revised/iteration5/Code/Figures/draw_Figure2_animation.R
Rscript --vanilla revised/iteration5/Code/Figures/validate_Figure2_animation.R
```

The animation command supports `--stills` for stage PNGs and the storyboard only. `FIGURE2_ANIMATION_OUTPUT_DIR` optionally changes the animation output folder. All paths are resolved relative to the script; the original presentation and local review package are not needed. The final animation stage uses the maintained `draw_Figure2.R`, not a separate copy of the model diagram.

Requirements: R with Cairo graphics, R packages `jsonlite`, `officer`, `xml2`, `zip`, and `png` (validation), plus `ffmpeg`, `ffprobe`, and `shasum` (validation). Static Figure 2 itself uses base R/grid. FFmpeg encoding uses two threads; intermediate frames are created in a temporary directory and removed automatically.

## Validation

The approved draft was checked visually at every stage, at selected transitions, and in the MP4 and rendered PowerPoint. The core-feedback slide is readable without overlap. Source and image comparisons verify that the animation endpoint matches the approved panel. The relocated scripts are rebuilt and checked independently in a clean branch worktree before integration.

The validator checks static figure regeneration, nonblank stage PNGs, MP4 dimensions/duration/frame rate, decoded holds against all seven stage images, 16:9 PowerPoint geometry, seven speaker-note pages, click-only fades, and byte-identical embedded stage images. Results and relative-path hash manifests are written beside the animation.
