# Core supporting results after adding `joint_pre.zip`

This document audits the manuscript-supporting evidence using three inputs: the current manuscript in the GitHub `HypoxiaLTEEFigures` branch, the top result archive `top10.zip`, and the joint-fitting preprocessing archive `joint_pre.zip`. Paths below are relative to the extracted input roots.

## 1. Input scope and source hierarchy

- `top10/top10_index.tsv` indexes the retained top results: top 10 separate in vitro fits, top 10 separate in vivo fits, and top 10 joint optimizer seeds for each of six landscape-informed warm-start pairs.

- `joint_pre.zip` supplies the 500-seed parameter landscape, t-SNE cluster/subcluster assignment, and fixed-O2 dense-grid monotonicity classification used to justify why the six joint warm starts were selected.

- The revised manuscript should therefore treat the joint fits as landscape-informed follow-up fits, not as isolated single-seed examples.

## 2. Optimization stability and top-result index

| result_group | scope | n | objective_min | objective_median | objective_max | objective_spread | relative_spread_pct |
| --- | --- | --- | --- | --- | --- | --- | --- |
| fit_invitro_O2_buffering_500seed | top10 | 10 | 3.8525 | 3.8602 | 3.8623 | 0.0098 | 0.2535 |
| fit_invivo_O2_buffering_500seed | top10 | 10 | 14.1193 | 14.1539 | 14.1748 | 0.0555 | 0.3928 |
| fit_joint_multi_warmup_tsne_sigmaN0p1216_objmin_500seed_20260714_021540 | all_joint_top10_per_pair | 60 | 18.8523 | 19.1937 | 20.0192 | 1.1670 | 6.1900 |
| fit_joint_multi_warmup_tsne_sigmaN0p1216_objmin_500seed_20260714_021540 | fit_joint_tsne_vi_seed366_C01Sc01_vt_seed10 | 10 | 18.8523 | 18.8606 | 18.8666 | 0.0143 | 0.0761 |
| fit_joint_multi_warmup_tsne_sigmaN0p1216_objmin_500seed_20260714_021540 | fit_joint_tsne_vi_seed322_C02Sc02_vt_seed10 | 10 | 18.8901 | 18.9138 | 18.9256 | 0.0355 | 0.1880 |
| fit_joint_multi_warmup_tsne_sigmaN0p1216_objmin_500seed_20260714_021540 | fit_joint_tsne_vi_seed25_C02Sc01_vt_seed10 | 10 | 18.9705 | 18.9728 | 18.9730 | 0.0025 | 0.0132 |
| fit_joint_multi_warmup_tsne_sigmaN0p1216_objmin_500seed_20260714_021540 | fit_joint_tsne_vi_seed311_C03Sc02_vt_seed10 | 10 | 19.4145 | 19.4264 | 19.4342 | 0.0197 | 0.1016 |
| fit_joint_multi_warmup_tsne_sigmaN0p1216_objmin_500seed_20260714_021540 | fit_joint_tsne_vi_seed290_C01Sc02_vt_seed10 | 10 | 19.7913 | 19.7990 | 19.8055 | 0.0142 | 0.0718 |
| fit_joint_multi_warmup_tsne_sigmaN0p1216_objmin_500seed_20260714_021540 | fit_joint_tsne_vi_seed138_C03Sc01_vt_seed10 | 10 | 19.9782 | 20.0138 | 20.0192 | 0.0410 | 0.2053 |

**Conclusion supported.** The separate in vitro top 10 solutions are numerically tight (objective spread 0.0098; 0.254% relative to the best), and the separate in vivo top 10 solutions are also tight (objective spread 0.0555; 0.393% relative). The joint warm-start pairs, however, differ substantially in objective, so joint conclusions should prioritize the best-performing warm-start pair and report the pair-level spread.

**Why this is significant.** Tight separate-fit top 10 objectives indicate that the best separate solutions are not one-off optimizer accidents. The larger across-pair joint spread indicates that the starting in vivo landscape region matters biologically and numerically.

## 3. `joint_pre.zip` landscape clustering used for joint warm starts

### 3.1 Primary clusters

| dataset | cluster_id | n_seeds | objective_min | objective_median | selected_k | selected_average_silhouette |
| --- | --- | --- | --- | --- | --- | --- |
| invivo | vi_C01 | 99 | 14.1340 | 15.0195 | 3 | 0.8157 |
| invivo | vi_C02 | 385 | 14.1193 | 14.7938 | 3 | 0.8157 |
| invivo | vi_C03 | 16 | 14.6048 | 15.1217 | 3 | 0.8157 |
| invitro | vt_C01 | 108 | 4.2712 | 4.4076 | 2 | 0.7742 |
| invitro | vt_C02 | 392 | 3.8525 | 4.0308 | 2 | 0.7742 |

**Conclusion supported.** The in vivo best-fit seed landscape separates into three primary clusters (`vi_C01`, `vi_C02`, `vi_C03`; selected k=3, average silhouette 0.8157), whereas the in vitro landscape separates into two primary clusters (`vt_C01`, `vt_C02`; selected k=2, average silhouette 0.7742).

**Why this is significant.** The joint warm starts are justified by a structured 500-seed preprocessing analysis rather than arbitrary hand-picking. The top 10 in vitro seeds all fall in `vt_C02_Sc02`; the top 10 in vivo seeds span `vi_C01_Sc01`, `vi_C02_Sc01`, and `vi_C02_Sc02`, meaning the in vivo fit landscape remains biologically heterogeneous even among low-objective solutions.

### 3.2 Top-10 seed cluster membership

| result_group | seed | rank | objective | dataset | primary_cluster_id | subcluster_id |
| --- | --- | --- | --- | --- | --- | --- |
| fit_invitro_O2_buffering_500seed | seed10 | 1 | 3.8525 | invitro | vt_C02 | vt_C02_Sc02 |
| fit_invitro_O2_buffering_500seed | seed132 | 2 | 3.8533 | invitro | vt_C02 | vt_C02_Sc02 |
| fit_invitro_O2_buffering_500seed | seed81 | 3 | 3.8541 | invitro | vt_C02 | vt_C02_Sc02 |
| fit_invitro_O2_buffering_500seed | seed294 | 4 | 3.8594 | invitro | vt_C02 | vt_C02_Sc02 |
| fit_invitro_O2_buffering_500seed | seed337 | 5 | 3.8598 | invitro | vt_C02 | vt_C02_Sc02 |
| fit_invitro_O2_buffering_500seed | seed106 | 6 | 3.8605 | invitro | vt_C02 | vt_C02_Sc02 |
| fit_invitro_O2_buffering_500seed | seed317 | 7 | 3.8610 | invitro | vt_C02 | vt_C02_Sc02 |
| fit_invitro_O2_buffering_500seed | seed140 | 8 | 3.8610 | invitro | vt_C02 | vt_C02_Sc02 |
| fit_invitro_O2_buffering_500seed | seed285 | 9 | 3.8618 | invitro | vt_C02 | vt_C02_Sc02 |
| fit_invitro_O2_buffering_500seed | seed464 | 10 | 3.8623 | invitro | vt_C02 | vt_C02_Sc02 |
| fit_invivo_O2_buffering_500seed | seed25 | 1 | 14.1193 | invivo | vi_C02 | vi_C02_Sc01 |
| fit_invivo_O2_buffering_500seed | seed366 | 2 | 14.1340 | invivo | vi_C01 | vi_C01_Sc01 |
| fit_invivo_O2_buffering_500seed | seed292 | 3 | 14.1372 | invivo | vi_C02 | vi_C02_Sc01 |
| fit_invivo_O2_buffering_500seed | seed392 | 4 | 14.1406 | invivo | vi_C01 | vi_C01_Sc01 |
| fit_invivo_O2_buffering_500seed | seed90 | 5 | 14.1524 | invivo | vi_C02 | vi_C02_Sc01 |
| fit_invivo_O2_buffering_500seed | seed391 | 6 | 14.1553 | invivo | vi_C01 | vi_C01_Sc01 |
| fit_invivo_O2_buffering_500seed | seed264 | 7 | 14.1558 | invivo | vi_C02 | vi_C02_Sc01 |
| fit_invivo_O2_buffering_500seed | seed109 | 8 | 14.1724 | invivo | vi_C02 | vi_C02_Sc01 |
| fit_invivo_O2_buffering_500seed | seed322 | 9 | 14.1724 | invivo | vi_C02 | vi_C02_Sc02 |
| fit_invivo_O2_buffering_500seed | seed26 | 10 | 14.1748 | invivo | vi_C02 | vi_C02_Sc01 |

## 4. In vitro result support

### 4.1 Key paths

- `top10/fit_invitro_O2_buffering_500seed/seed10/fit_summary.tsv`

- `top10/fit_invitro_O2_buffering_500seed/seed10/invitro_lineage_summary.tsv`

- `top10/fit_invitro_O2_buffering_500seed/seed10/invitro_distribution_quantiles.tsv`


### 4.2 Branch snapshots from the best in vitro seed (`seed10`)

| cohort | branch | snapshot | passage_index | oxygen_pct | predicted_growth_rate | predicted_mean_kary_N | observed_growth | observed_mean_kary_N | observed_n_kary | observed_n_flow |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 2N | deprived | first | 2 | 20.5000 | 1.0554 | 46.6295 | nan | nan | 0 | 0 |
| 2N | deprived | last | 38 | 0.0000 | 0.3786 | 66.7229 | 0.4809 | nan | 0 | 200 |
| 2N | deprived | first_observed_kary | 3 | 20.5000 | 1.0556 | 46.7180 | nan | 48.1500 | 20 | 0 |
| 2N | deprived | last_observed_kary | 34 | 0.0000 | 0.3707 | 64.2664 | 0.4702 | 88.0500 | 20 | 0 |
| 2N | normoxia_like | first | 1 | 20.5000 | 0.9648 | 46.5914 | nan | 47.3000 | 20 | 0 |
| 2N | normoxia_like | last | 1 | 20.5000 | 0.9648 | 46.5914 | nan | 47.3000 | 20 | 0 |
| 2N | normoxia_like | first_observed_kary | 1 | 20.5000 | 0.9648 | 46.5914 | nan | 47.3000 | 20 | 0 |
| 2N | normoxia_like | last_observed_kary | 1 | 20.5000 | 0.9648 | 46.5914 | nan | 47.3000 | 20 | 0 |
| 4N | deprived | first | 2 | 20.5000 | 1.1124 | 84.3970 | nan | nan | 0 | 0 |
| 4N | deprived | last | 37 | 0.0000 | 0.3877 | 80.7575 | 0.4850 | nan | 0 | 200 |
| 4N | deprived | first_observed_kary | 5 | 20.5000 | 1.1128 | 84.6553 | nan | 86.1250 | 16 | 0 |
| 4N | deprived | last_observed_kary | 32 | 0.0000 | 0.3877 | 81.5801 | 0.4050 | 81.3000 | 20 | 0 |
| 4N | normoxia_like | first | 1 | 20.5000 | 0.9765 | 84.3383 | nan | 88.9000 | 20 | 0 |
| 4N | normoxia_like | last | 1 | 20.5000 | 0.9765 | 84.3383 | nan | 88.9000 | 20 | 0 |
| 4N | normoxia_like | first_observed_kary | 1 | 20.5000 | 0.9765 | 84.3383 | nan | 88.9000 | 20 | 0 |
| 4N | normoxia_like | last_observed_kary | 1 | 20.5000 | 0.9765 | 84.3383 | nan | 88.9000 | 20 | 0 |

### 4.3 Selected predicted distribution quantiles

| cohort | passage_index | segment_id | oxygen_pct | 0.01 | 0.25 | 0.75 | 0.99 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 2N | 1 | 20.5 | 20.50 | 40 | 44 | 48 | 55 |
| 2N | 3 | 20.5_20.5_20.5 | 20.50 | 40 | 44 | 48 | 88 |
| 2N | 16 | 20.5_20.5_2_1_0.5_0.3_0.2_0.1_0 | 0.00 | 38 | 44 | 50 | 101 |
| 2N | 32 | 20.5_20.5_2_1_0.5_0.3_0.2_0.1_0_0_0_0_0_0_0_0_0_0_0 | 0.00 | 35 | 44 | 82 | 107 |
| 2N | 34 | 20.5_20.5_2_1_0.5_0.3_0.2_0.1_0_0_0_0_0_0_0_0_0_0_0_0_0 | 0.00 | 34 | 44 | 84 | 107 |
| 2N | 37 | 20.5_20.5_2_1_0.5_0.3_0.2_0.1_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0 | 0.00 | 32 | 44 | 84 | 107 |
| 2N | 38 | 20.5_20.5_2_1_0.5_0.3_0.2_0.1_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0 | 0.00 | 32 | 44 | 85 | 107 |
| 4N | 1 | 20.5 | 20.50 | 73 | 81 | 88 | 97 |
| 4N | 3 | 20.5_20.5_20.5 | 20.50 | 73 | 81 | 88 | 98 |
| 4N | 16 | 20.5_20.5_2_1_0.5_0.3_0.2_0.1_0 | 0.00 | 67 | 78 | 89 | 103 |
| 4N | 32 | 20.5_20.5_2_1_0.5_0.3_0.2_0.1_0_0_0_0_0_0_0_0_0_0_0 | 0.00 | 62 | 75 | 87 | 105 |
| 4N | 34 | 20.5_20.5_2_1_0.5_0.3_0.2_0.1_0_0_0_0_0_0_0_0_0_0_0_0_0 | 0.00 | 62 | 75 | 87 | 106 |
| 4N | 37 | 20.5_20.5_2_1_0.5_0.3_0.2_0.1_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0 | 0.00 | 61 | 74 | 87 | 106 |

**Conclusion supported.** The separate in vitro fit supports oxygen-dependent growth suppression and broadening/remodeling of chromosome-number distributions. In the 2N deprived branch, the predicted distribution broadens markedly by late passages, with q25 remaining near 44 while q75 reaches 85 and q99 reaches 107 at passage 38. In the 4N deprived branch, the predicted mean declines from 84.66 at the early 4N matched point to 80.76 by the last deprived snapshot.

**Caution.** The observed 2N deprived A19 replicate means are split (66.85 and 88.05), while the model predicts 64.27 for the shared segment. Therefore the manuscript should say the fit captures broad high-ploidy expansion/remodeling, not that it perfectly recapitulates every high-ploidy replicate.

## 5. In vivo result support

### 5.1 Key paths

- `top10/fit_invivo_O2_buffering_500seed/seed25/fit_summary.tsv`

- `top10/fit_invivo_O2_buffering_500seed/seed25/terminal_ploidy_fit.tsv`

- `top10/fit_invivo_O2_buffering_500seed/seed25/burden_fit.tsv`

- `top10/fit_invivo_O2_buffering_500seed/seed25/necrosis_fit.tsv`


### 5.2 Terminal ploidy fit in the best separate in vivo seed

| harvest | cohort | obs_n | obs_mean_N | pred_mean_N | pred_mode_N | pred_mode_fraction | pred_frac_N_le44 | pred_frac_N_ge66 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| SUM159-2N-0-O_harvest | 2N | 604 | 45.9272 | 46.0876 | 44 | 0.0818 | 0.4443 | 0.0059 |
| SUM159-2N-0-RL_harvest | 2N | 547 | 47.7349 | 46.0876 | 44 | 0.0818 | 0.4443 | 0.0059 |
| SUM159-2N-0-RR_harvest | 2N | 2274 | 40.4763 | 46.0876 | 44 | 0.0818 | 0.4443 | 0.0059 |
| SUM159-2N-0-R_harvest | 2N | 303 | 43.0693 | 46.0876 | 44 | 0.0818 | 0.4443 | 0.0059 |
| SUM159-4N-0-0_harvest | 4N | 1097 | 53.6618 | 48.8653 | 47 | 0.0740 | 0.2427 | 0.0136 |
| SUM159-4N-0-L_harvest | 4N | 636 | 49.8381 | 50.8012 | 48 | 0.0676 | 0.1568 | 0.0296 |
| SUM159-4N-0-RR_harvest | 4N | 529 | 47.5369 | 48.8653 | 47 | 0.0740 | 0.2427 | 0.0136 |
| SUM159-4N-0-R_harvest | 4N | 512 | 49.0684 | 50.8012 | 48 | 0.0676 | 0.1568 | 0.0296 |

Cohort-level summary:

| cohort | obs_mean | pred_mean | pred_frac_le44 | pred_frac_ge66 |
| --- | --- | --- | --- | --- |
| 2N | 44.3019 | 46.0876 | 0.4443 | 0.0059 |
| 4N | 50.0263 | 49.8333 | 0.1998 | 0.0216 |

**Conclusion supported.** The best separate in vivo seed places terminal viable chromosome-number distributions near the diploid range: predicted cohort means are 46.09 for 2N-initiated tumors and approximately 48.9-50.8 for 4N-initiated tumors, close to observed cohort averages. High-ploidy mass (N >= 66) is small in the predicted terminal distributions.

**Why this is significant.** This supports the manuscript claim that the tumor context favors lower chromosome-number states relative to the initial 4N lineage. However, it does not by itself identify one unique oxygen-response mechanism; the fixed-O2 preprocessing shows multiple response classes.

### 5.3 Burden fit audit, positive observations only

| cohort | harvest | n | log_rmse | obs_last | pred_last |
| --- | --- | --- | --- | --- | --- |
| 2N | SUM159-2N-0-O_harvest | 10 | 0.6411 | 1729.6470 | 1968.4619 |
| 2N | SUM159-2N-0-RL_harvest | 10 | 0.2687 | 2622.3600 | 1968.4619 |
| 2N | SUM159-2N-0-RR_harvest | 10 | 0.2854 | 3375.3520 | 1968.4619 |
| 2N | SUM159-2N-0-R_harvest | 10 | 0.4060 | 1683.1360 | 1968.4619 |
| 4N | SUM159-4N-0-0_harvest | 11 | 0.5488 | 3914.1440 | 3125.9299 |
| 4N | SUM159-4N-0-L_harvest | 14 | 1.0972 | 5899.2780 | 1589.4592 |
| 4N | SUM159-4N-0-RR_harvest | 14 | 0.9490 | 2423.0250 | 3125.9299 |
| 4N | SUM159-4N-0-R_harvest | 13 | 0.2390 | 2096.9760 | 1589.4592 |

## 6. Fixed-O2 attractor and response-class evidence from `joint_pre.zip`

### 6.1 Class counts across all 500 in vivo seeds

| curve_class | n_seed | fraction_seed |
| --- | --- | --- |
| complex_nonmonotone | 261 | 0.5220 |
| inverted_u_shaped | 149 | 0.2980 |
| monotone_increasing | 57 | 0.1140 |
| u_shaped | 25 | 0.0500 |
| single_transition_increase_then_plateau | 5 | 0.0100 |
| monotone_decreasing | 3 | 0.0060 |

### 6.2 Class counts among the top 25 in vivo objectives

| curve_class | n_seed | fraction | top_n |
| --- | --- | --- | --- |
| complex_nonmonotone | 9 | 0.3600 | 25 |
| inverted_u_shaped | 10 | 0.4000 | 25 |
| monotone_increasing | 2 | 0.0800 | 25 |
| u_shaped | 4 | 0.1600 | 25 |

### 6.3 Top-10 in vivo fixed-O2 response classes

| seed | rank | objective_top10 | curve_class | monotonicity_reliability | dominant_mean_ploidy_o2_0 | dominant_mean_ploidy_o2_0p5 | dominant_mean_ploidy_o2_1 | dominant_mean_ploidy_o2_5 | min_spectral_gap |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| seed25 | 1 | 14.1193 | inverted_u_shaped | caution_small_gap | 2.1718 | 3.0148 | 1.0000 | 1.0001 | 0.0003 |
| seed366 | 2 | 14.1340 | inverted_u_shaped | unreliable_small_gap | 2.1140 | 2.7755 | 1.0016 | 1.0022 | 0.0004 |
| seed292 | 3 | 14.1372 | inverted_u_shaped | unreliable_small_gap | 2.4025 | 2.9307 | 1.0021 | 1.0164 | 0.0006 |
| seed392 | 4 | 14.1406 | inverted_u_shaped | caution_small_gap | 4.0640 | 1.0000 | 1.0000 | 1.0001 | 0.0000 |
| seed90 | 5 | 14.1524 | monotone_decreasing | caution_small_gap | 1.0003 | 1.0001 | 1.0001 | 1.0002 | 0.0063 |
| seed391 | 6 | 14.1553 | inverted_u_shaped | caution_small_gap | 3.6441 | 4.2834 | 1.0002 | 1.0001 | 0.0000 |
| seed264 | 7 | 14.1558 | inverted_u_shaped | unreliable_small_gap | 2.2517 | 3.9970 | 1.0001 | 1.0018 | 0.0009 |
| seed109 | 8 | 14.1724 | complex_nonmonotone | caution_small_gap | 1.3640 | 1.0001 | 1.0001 | 2.4634 | 0.0000 |
| seed322 | 9 | 14.1724 | complex_nonmonotone | caution_small_gap | 2.2666 | 2.2940 | 2.4566 | 2.4094 | 0.0004 |
| seed26 | 10 | 14.1748 | inverted_u_shaped | caution_small_gap | 3.9771 | 4.4686 | 1.0000 | 1.0000 | 0.0000 |

**Conclusion supported.** The fixed-O2 attractor analysis should be presented as model uncertainty rather than a single selected O2-ploidy law. Across all 500 seeds, 52.2% are complex nonmonotone, 29.8% inverted-U, 11.4% monotone-increasing, and smaller fractions fall into u-shaped, plateau, or monotone-decreasing classes. Among the top 25 objectives, inverted-U and complex nonmonotone remain dominant but mixed.

**Why this is significant.** This prevents overclaiming that oxygen alone monotonically selects for or against WGD. It is a core result because it directly constrains the biological interpretation of the in vivo landscape.

## 7. Joint soft-coupled result support

### 7.1 Best joint seed per warm-start pair

| pair | seed | objective |
| --- | --- | --- |
| fit_joint_tsne_vi_seed366_C01Sc01_vt_seed10 | seed472 | 18.8523 |
| fit_joint_tsne_vi_seed322_C02Sc02_vt_seed10 | seed54 | 18.8901 |
| fit_joint_tsne_vi_seed25_C02Sc01_vt_seed10 | seed497 | 18.9705 |
| fit_joint_tsne_vi_seed311_C03Sc02_vt_seed10 | seed18 | 19.4145 |
| fit_joint_tsne_vi_seed290_C01Sc02_vt_seed10 | seed155 | 19.7913 |
| fit_joint_tsne_vi_seed138_C03Sc01_vt_seed10 | seed122 | 19.9782 |

### 7.2 Best overall joint seed soft-coupled parameter ratios

| parameter | vivo_projected_natural | vitro_projected_natural | ratio_vivo_to_vitro | delta_transformed |
| --- | --- | --- | --- | --- |
| O2_crit | 0.22251 | 1.58621 | 0.14028 | -0.85302 |
| mu_hp | 0.05135 | 0.00500 | 10.26956 | 1.01155 |
| p_misseg | 0.11736 | 0.01055 | 11.12024 | 1.04611 |
| k_o_mis | 0.00103 | 0.00011 | 9.11335 | 0.95968 |
| buffer_smax | 0.90852 | 1.00000 | 0.90852 | -0.09148 |
| buffer_beta | 0.16539 | 1.59113 | 0.10394 | -0.98320 |
| buffer_n_exp | 10.00000 | 3.15757 | 3.16700 | 0.50065 |
| n_O | 2.17248 | 1.86857 | 1.16264 | 0.30391 |
| alpha_o2 | 2.69193 | 0.50000 | 5.38387 | 0.73109 |
| gamma_growth | 3.30759 | 1.50000 | 2.20506 | 1.80759 |
| lam_max | 0.23892 | 1.34852 | 0.17717 | -0.75161 |
| p_mis_base | 0.00002 | 0.00002 | 1.03210 | 0.01372 |
| p_wgd | 0.00025 | 0.00024 | 1.03465 | 0.01479 |
| gamma_mu | 2.63751 | 1.20000 | 2.19792 | 1.43751 |

### 7.3 Median ratios across the best seed from each of the six warm-start pairs

| parameter | n | ratio_median | ratio_q25 | ratio_q75 | vivo_median | vitro_median |
| --- | --- | --- | --- | --- | --- | --- |
| lam_max | 6 | 0.17713 | 0.11879 | 0.18149 | 0.23886 | 1.34852 |
| p_misseg | 6 | 16.83941 | 12.48715 | 39.95127 | 0.17771 | 0.01055 |
| mu_hp | 6 | 16.57977 | 10.28930 | 28.04270 | 0.08290 | 0.00500 |
| buffer_smax | 6 | 0.91118 | 0.90521 | 0.95216 | 0.91118 | 1.00000 |
| buffer_beta | 6 | 0.07081 | 0.03644 | 0.09781 | 0.11266 | 1.59113 |
| buffer_n_exp | 6 | 2.87166 | 0.85677 | 3.16700 | 9.06744 | 3.15757 |
| alpha_o2 | 6 | 4.43231 | 2.80213 | 6.13078 | 2.21616 | 0.50000 |
| gamma_growth | 6 | 1.70356 | 1.54597 | 1.85824 | 2.55535 | 1.50000 |
| p_mis_base | 6 | 32.16765 | 2.99085 | 1035.41903 | 0.00026 | 0.00001 |
| p_wgd | 6 | 0.98901 | 0.86598 | 1.03032 | 0.00028 | 0.00033 |
| O2_crit | 6 | 0.16812 | 0.11486 | 0.29399 | 0.25927 | 1.58621 |
| n_O | 6 | 1.55964 | 1.21808 | 1.90858 | 2.91430 | 1.86857 |
| gamma_mu | 6 | 2.18229 | 2.13742 | 2.22368 | 2.65811 | 1.20000 |

**Conclusion supported.** The joint fit supports three context-specific axes: lower effective proliferative ceiling in vivo, stronger stress-linked CIN/missegregation in vivo, and a different post-missegregation survival filter. In the best overall joint seed, `lam_max` is lower in vivo than in vitro (0.239 vs 1.349; ratio 0.177), while `p_misseg` and `mu_hp` are higher in vivo (ratios 11.12 and 10.27). Across the best seed from each warm-start pair, median ratios remain directionally similar for `lam_max`, `p_misseg`, `mu_hp`, `alpha_o2`, and `gamma_growth`.

**Caution.** Several parameters are at or near bounds in the joint projection, and the Welsch penalty often saturates. The manuscript should interpret these as model-resolved hypotheses, not direct measurements of oxygen, CIN, or survival.

## 8. Necrosis output audit

| result_group | n_rows | n_obs | n_pred_finite | all_pred_nan_rate |
| --- | --- | --- | --- | --- |
| fit_invivo_O2_buffering_500seed | 60 | 60 | 0 | 1.0000 |
| fit_joint_multi_warmup_tsne_sigmaN0p1216_objmin_500seed_20260714_021540 | 360 | 360 | 0 | 1.0000 |

**Conclusion.** Necrosis is currently not usable as direct support in the top-result files: all 60 separate in vivo rows and all 360 joint top-result rows have observed necrosis values but no finite predicted necrosis fractions. The manuscript can describe necrosis as an intended observation model only if the implementation/export problem is fixed, or it must state that necrosis was not successfully audited from the current output package.

## 9. Recommended manuscript-level claims that are adequately supported

1. Matched in vitro and in vivo SUM159 lineage data motivate a model in which resource stress can alter proliferation, chromosome missegregation, WGD generation, and daughter survival.

2. Separate in vitro fits support oxygen-dependent growth suppression and distributional broadening/remodeling rather than simple high-ploidy killing.

3. Separate in vivo fits support terminal convergence toward lower chromosome-number states, especially for 4N-initiated tumors.

4. The 500-seed landscape supports multiple in vivo response regimes; fixed O2 does not select a unique monotone outcome.

5. Joint soft-coupled fits support context-specific differences, especially lower in vivo proliferation ceiling and higher in vivo stress-linked missegregation, but these are model-resolved effects and should be framed as hypotheses.


## 10. Claims that should be weakened or held until repaired

1. Any claim that terminal necrosis independently constrains the fitted trajectories should be withheld until finite necrosis predictions are exported and reconciled.

2. Any claim of a single oxygen-to-ploidy rule should be replaced by a response-class uncertainty statement.

3. Any claim that the in vitro fit perfectly recapitulates the late 2N deprived high-ploidy replicate should be softened; the current fit captures broadening but misses one high observed mean.
