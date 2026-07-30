# Revision difference record versus GitHub manuscript (`HypoxiaLTEEFigures/manuscript/ltee_hypoxia_model.tex`)

Original source identified from GitHub connector: `thecloneredesignlab/miningcloneid`, branch `HypoxiaLTEEFigures`, file `manuscript/ltee_hypoxia_model.tex`, SHA `aeb5c7a9bc0c1de642834e77c437c660a08a4a0b`.

## High-level differences

1. Removed unresolved placeholder text in the Introduction, including the CLONEID/sister-paper placeholder, the PubMed note placeholder, and the mutator-theory draft block that was not connected to the current results. The revised Introduction keeps the original scientific axis: resource stress, WGD/aneuploidy cost, CIN generation, and ploidy-dependent buffering.

2. Rewrote the Methods to preserve the original model logic while explicitly adding the evidence now available from `joint_pre.zip`: 500-seed t-SNE landscape clustering, subcluster-based joint warm starts, and dense fixed-O2 monotonicity classification.

3. Replaced overstrong result statements with output-audited statements. In vitro conclusions now state that the best fit supports oxygen-dependent distributional remodeling and broadening; it no longer claims perfect reconstruction of every high-ploidy replicate.

4. Expanded the in vivo Results with 500-seed fixed-O2 class counts and top-objective class mixing. This changes the interpretation from a single oxygen-response law to retained model uncertainty across several plausible response regimes.

5. Expanded the joint-fit Results with the actual best joint warm-start pair (`vi_seed366_C01Sc01` paired with `vt_seed10`) and natural-scale in vivo/in vitro parameter ratios from `joint_soft_coupling_projection.tsv`.

6. Added an explicit necrosis audit limitation. The original manuscript described necrosis as part of the observation model, but the current retained result tables have observed necrosis rows with no finite predicted necrosis fractions. The revised manuscript does not use necrosis as independent support until this is repaired.

7. Added a focused Discussion limitation section covering necrosis export, in vitro replicate mismatch, soft-coupling/bound saturation, and latent-O2 interpretation.


## Result-level differences

- Separate in vitro stability added: top 10 objective spread 0.0098 (0.254%).

- Separate in vivo stability added: top 10 objective spread 0.0555 (0.393%).

- In vitro branch evidence added: late 2N deprived distribution q25=44, q75=85, q99=107 at passage 38; 4N deprived predicted mean declines from 84.7 to 80.8.

- In vivo terminal evidence added: 2N predicted mean 46.1 versus observed cohort mean 44.3; 4N predicted means 48.9--50.8 versus observed cohort mean 50.0.

- Fixed-O2 response class counts added: complex nonmonotone 52.2%, inverted-U 29.8%, monotone-increasing 11.4%, u-shaped 5.0%, plateau 1.0%, monotone-decreasing 0.6%.

- Joint context ratios added: in best overall joint seed, `lam_max` ratio in vivo/in vitro=0.177, `p_misseg` ratio=11.12, `mu_hp` ratio=10.27, `alpha_o2` ratio=5.38, `gamma_growth` ratio=2.21.


## Files created by this revision

- `revised_manuscript/revised_ltee_hypoxia_model.tex`: completed LaTeX manuscript.

- `revised_manuscript/revised_ltee_hypoxia_model_readable.md`: plain-readable copy of the revised manuscript text.

- `analysis_tables/*.csv`: derived audit tables used to support the revision.
