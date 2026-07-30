# Expansion notes for revised manuscript

The following content was added or expanded because `joint_pre.zip` contains preprocessing evidence that was not fully represented in the original manuscript draft.

## Added/expanded content

1. **Landscape-informed joint warm-start rationale.** The revised manuscript explicitly states that the six joint warm starts came from the 500-seed in vivo t-SNE landscape: `vi_C01_Sc01`, `vi_C01_Sc02`, `vi_C02_Sc01`, `vi_C02_Sc02`, `vi_C03_Sc01`, and `vi_C03_Sc02`, paired with the best in vitro seed (`seed10`). This is supported by `joint_pre/landscape_subcluster/.../best_subcluster_summary.csv` and `top10/top10_index.tsv`.

2. **Fixed-O2 response-class uncertainty.** The revised manuscript adds a result paragraph summarizing the dense fixed-O2 classification: 52.2% complex nonmonotone, 29.8% inverted-U, 11.4% monotone-increasing, 5.0% u-shaped, 1.0% plateau-type, and 0.6% monotone-decreasing. This justifies treating O2-ploidy response as a family of plausible regimes rather than a single selected curve.

3. **Numerical evidence for joint context differences.** The revised manuscript adds the best overall joint seed and cross-pair median ratios for key soft-coupled parameters, including lower in vivo `lam_max`, higher in vivo `p_misseg`, higher in vivo `mu_hp`, and altered buffering parameters.

4. **Evidence limits.** The revised manuscript explicitly flags the necrosis export problem and avoids using finite necrosis predictions as independent support. It also softens the in vitro 2N severe-deprivation interpretation because one late observed replicate is higher than the model prediction.


## Why these additions are justified

The original draft already framed the work around separate in vivo, separate in vitro, joint soft-coupled fitting, and fixed-O2 attractor diagnostics. The new text does not add a new scientific direction; it documents the actual preprocessing evidence needed to substantiate those existing sections.

## Meaning for the manuscript

The strongest manuscript is not a single-response oxygen story. It is a model-resolution and uncertainty story: in vitro and in vivo fits require different stress-response parameterization, but the in vivo parameter landscape supports multiple oxygen-CIN-ploidy regimes.
