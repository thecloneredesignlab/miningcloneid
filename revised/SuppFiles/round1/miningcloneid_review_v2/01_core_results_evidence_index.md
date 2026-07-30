# 核心结果与证据索引（纳入 joint_pre.zip 的重新审查版）

## 1. 审查范围、结果层级与“显著性”的含义

本索引以三类材料为唯一结果依据：

1. Repo 分支 `HypoxiaLTEEFigures` 中的手稿、模型代码、配置、后处理与报告脚本；
2. `top10.zip` 中的独立 *in vitro* top-10、独立 *in vivo* top-10，以及 6 组 joint warm-start 配对各自的 top-10，共 60 个 selected joint solutions；
3. `joint_pre.zip` 中的 500-seed 参数景观、t-SNE/cluster/subcluster 前处理、warm-start 选择记录、full-500 fixed-O2 dominant-eigenvector 结果、曲线分类、spectral-gap 可靠性和 objective crosswalk。

下文“显著性”表示**科学意义及其对手稿主结论的重要程度**，不是统计学显著性。当前 top-10/top-60 是按目标函数选出的近优解，不是后验样本，也没有构成独立生物学重复或正式置信区间。因此，本索引只使用“跨 selected solutions 方向稳定”“函数级结果一致”“证据强/中/弱”等术语，不虚构 `P` 值。

便携复核表位于交付包的 `analysis/` 和 `analysis_tables/` 目录。关键汇总为 `analysis/analysis_summary.json`。

---

## 2. 最核心且可以保留的结果

### 结果 1：4N 严重资源剥夺支系的染色体数下降在 *in vitro* top-10 中稳定

**原始结果位置**

- `top10/fit_invitro_O2_buffering_500seed/seed*/fit_summary.tsv`
- `top10/fit_invitro_O2_buffering_500seed/seed*/invitro_lineage_summary.tsv`
- `top10/fit_invitro_O2_buffering_500seed/seed*/invitro_daily_counts.tsv`
- `top10/fit_invitro_O2_buffering_500seed/seed*/report/fit_report_seed*.html`
- 便携复核：`analysis/invitro_top10_robustness.csv`

**定量证据**

- top-10 objective：`3.852535–3.862302`；
- 最长 severe-deprivation 4N 支系进入该阶段时的预测平均染色体数：约 `84.21–84.71`；
- 终末预测平均染色体数：`80.13–80.89`；
- 下降：`3.54–4.17` 条染色体，中位数 `3.86`。

**支持结论**

在当前模型与数据约束下，4N 支系在严重资源剥夺阶段出现方向一致的温和 chromosome-number decline，不依赖单个 optimizer seed。

**为何有科学意义**

这是论文“高染色体状态在资源压力下可以通过代际重塑而下降”的直接结果基础。它把单一最佳解的示例提升为 selected near-optimal solutions 中可复现的趋势。

**证据等级与边界**

- 证据等级：较强；
- 边界：top-10 不是统计置信区间，也不能证明现实细胞中发生了特定染色体逐条丢失路径。

---

### 结果 2：上述 4N 下降不需要大量直接 hypoxia-origin killing

**原始结果位置**

- `top10/fit_invitro_O2_buffering_500seed/seed*/invitro_daily_counts.tsv`
- `top10/fit_invitro_O2_buffering_500seed/seed*/invitro_lineage_summary.tsv`
- Repo 的 live、hypoxia-dead、CIN-dead compartment 更新逻辑
- 便携复核：`analysis/invitro_top10_robustness.csv`

**定量证据**

- top-10 中最大 direct hypoxia-origin dead fraction：`1.739–1.815%`；
- 最大 CIN-associated dead fraction：`32.15–33.40%`。

**支持结论**

模型并不需要通过直接杀死大部分 4N 活细胞来实现平均染色体数下降。主要选择发生在 chromosome-transition products 的生成、存活与移除层面。

**为何有科学意义**

这一区分将机制从“低氧直接淘汰高倍体”改为“资源压力改变错误分离后代的生成与过滤”，是手稿的核心机制贡献。

**关键限制**

`CIN-associated dead` 不是纯生物学死亡指标。它同时包括：

- post-missegregation nonviable daughter mass；
- 超出 `N_MIN=22`、`N_MAX=154` 状态网格的 daughter/WGD mass；
- 其他由 finite-grid boundary routing 产生的质量。

因此，32–33% 不能直接等同于组织学坏死或实验测得的细胞死亡。

---

### 结果 3：2N deprived lineage 产生高染色体数成分，但模型没有重现后期分支分化

**原始结果位置**

- `top10/fit_invitro_O2_buffering_500seed/seed*/invitro_distribution_summary.tsv`
- `.../invitro_observed_kary.tsv`
- `.../invitro_ploidy_loglik.tsv`
- `.../report/fit_report_seed*.html`
- 便携复核：`analysis/invitro_top10_robustness.csv`、`analysis_tables/invitro_shared_segment_branch_divergence.csv`

**定量证据**

- late deprived 2N 的 top-10 预测平均染色体数：`63.72–64.40`；
- 对应两个实验分支的观测均值：`66.85` 与 `88.05`；
- 当前 lineage-level simulation 对这两个后期分支给出共享预测，无法产生其显著分化。

**支持结论**

模型支持 deprived 2N lineage 中出现 high-chromosome/WGD-like component，并同时保留 lower-chromosome mode；但不支持“所有 2N deprived trajectories 均先整体升高再整体回落”，也不能解释后期两个实验分支的不同终点。

**为何有科学意义**

它明确了模型能解释的是 variation generation 和混合分布，而不是 branch-specific fate。该边界防止把分布内出现低染色体后代误写成总体 depolyploidization。

---

### 结果 4：独立 *in vivo* 拟合可联合使用 burden、terminal chromosome-number 与 necrosis，但绝对拟合质量仅为中等

**原始结果位置**

- `top10/fit_invivo_O2_buffering_500seed/seed*/fit_summary.tsv`
- `.../burden_fit.tsv`
- `.../terminal_ploidy_fit.tsv`
- `.../necrosis_fit.tsv`
- `.../report/fit_report_seed*.html`
- 便携复核：`analysis/invivo_top10_fit_quality.csv`

**定量证据**

- objective：`14.1193–14.1748`；
- all-observation log-burden RMSE：`0.663–0.680`；
- tumor-balanced log-burden RMSE：`0.628–0.644`；
- terminal mean chromosome-number MAE：`2.41–3.71`；
- terminal distribution Wasserstein-1：`4.18–5.24` chromosomes；
- terminal distribution total variation：`0.681–0.714`；
- reconstructed necrosis-fraction MAE：`0.064–0.103`。

**支持结论**

模型抓住 cohort-level burden 与 terminal mean chromosome-number 的总体尺度，并能够在同一个目标函数中使用三种数据流；但不能声称对每只小鼠的时间轨迹或完整终末分布形状实现高精度重现。

**为何有科学意义**

这使手稿把模型正确定位为 mechanism-inference framework，而不是 individual-tumor forecasting model。

---

### 结果 5：当前 *in vivo* 结构主要是 cohort/lineage-level，而非 mouse-specific model

**结果位置**

- Repo config：`harvest_init_multiplier: FALSE`；
- `top10/fit_invivo.../seed*/run_effective_args.tsv`；
- `.../burden_fit.tsv` 与 `.../terminal_ploidy_fit.tsv`；
- 便携复核：`analysis_tables/invivo_common_day_heterogeneity.csv`。

**结构证据**

相同 starting-ploidy cohort 且相同 harvest day 的小鼠共享同一模型预测，个体差异主要进入 observation residual。模型没有 mouse-specific latent oxygen、initial-size multiplier 或其他 individual random effect。

**支持结论**

拟合参数反映共享 lineage/context 机制，不是每只肿瘤的个体环境测量。

**为何有科学意义**

它限制了因果与个体化语言：可以提出 tumor-context hypothesis，不能把每条 latent O2 trajectory 当作实测的小鼠组织氧状态。

---

### 结果 6：独立 *in vivo* 近优解存在明显 practical non-identifiability 与边界依赖

**结果位置**

- `top10/fit_invivo.../seed*/best_params.tsv`
- `.../parameter_table_input.csv`
- `analysis/independent_parameter_values_long.csv`
- `analysis/independent_parameter_boundary_hits.csv`
- `analysis_tables/invivo_top10_parameter_range_summary.csv`

**定量/结构证据**

objective 仅相差约 `0.0555`，但多个 headline parameters 跨越宽范围。`sigma_burden` 在 10/10 selected *in vivo* fits 中达到上界；其他参数也在部分解中触及上下界。

**支持结论**

当前数据允许多个机制组合产生近似相同的拟合质量。单个参数数值不能被解释为独立、精确的生物学测量。

**为何有科学意义**

这是理解 fixed-O2 多响应类、joint warm-start dependence 和 pair-specific solutions 的基础。模型不确定性必须在 Discussion 中作为中心结果，而不是附带技术问题。

---

## 3. joint_pre.zip 新增并改变解释的结果

### 结果 7：full-500 参数景观存在清晰的低维分区，但 subcluster 生物学分离度有限

**原始结果位置**

`joint_pre.zip` 中与 parameter landscape 相关的：

- primary-cluster assignment/summary；
- selected-k silhouette tables；
- subcluster assignment/summary；
- warm-start representative selection tables；
- t-SNE coordinates 与 standardized parameter features。

便携复核：

- `analysis/landscape_primary_cluster_summary.csv`
- `analysis/landscape_primary_cluster_silhouette.csv`
- `analysis/landscape_subcluster_summary.csv`
- `analysis/landscape_subcluster_silhouette.csv`

**定量证据**

- *in vivo* primary t-SNE clustering：`k=3`，cluster sizes `99/385/16`，embedding-space average silhouette `0.816`；
- *in vitro* primary t-SNE clustering：`k=2`，sizes `108/392`，embedding-space average silhouette `0.774`；
- 用 standardized original features 在每个 *in vivo* primary cluster 内做 subclustering，其 silhouettes 为 `0.179`、`0.096`、`0.330`；
- *in vitro* subcluster silhouettes 为 `0.137`（vt_C01）和 `0.690`（vt_C02）。

**支持结论**

这些 clusters 可作为 parameter-space coverage strata，但不能直接命名为稳定生物亚型。primary silhouette 来自非线性 t-SNE embedding，不能替代原始高维空间中的生物分离验证；用于选择 *in vivo* warm starts 的原始特征 subclusters 分离度多为弱到中等。

**为何有科学意义**

它决定了 joint fitting 的代表性：六个 warm starts 是探索性覆盖方案，而不是经独立验证的六种生物机制。

---

### 结果 8：joint warm-start 设计在 *in vitro* 侧覆盖不平衡

**原始结果位置**

- `joint_pre.zip` 的 warm-start representative/crosswalk tables；
- `top10/fit_joint_multi_warmup.../fit_joint_tsne_vi_seed*_vt_seed10/`；
- 便携复核：`analysis/joint_warm_start_fixed_o2_crosswalk.csv`、`analysis_tables/joint_pre_invitro_subcluster_representatives_and_coverage.csv`。

**定量证据**

- 六个 *in vivo* representative seeds：`366, 290, 25, 322, 138, 311`；
- 每个为其 exploratory subcluster 内 objective-minimum seed，不是 medoid；
- 六组全部与同一个 *in vitro* seed 10 配对；
- seed 10 属于 `vt_C02_Sc02`，该 subcluster 包含 `376/500` *in vitro* fits；其他 *in vitro* subclusters 未进入 joint warm-start panel。

**支持结论**

joint panel 探索了六个 *in vivo* regions 对一个固定 *in vitro* solution 的响应，不是 *in vivo* × *in vitro* landscape 的平衡 factorial comparison。

**为何有科学意义**

所有 context difference 都可能同时依赖所选 *in vivo* region 和单一 *in vitro* anchor。方向一致性可以作为候选机制，但不能宣称已覆盖 culture parameter uncertainty。

---

### 结果 9：joint fitting 对六个 warm-start region 的 objective 明显敏感

**原始结果位置**

- `top10/fit_joint_multi_warmup.../<pair>/seed*/fit_summary.tsv`
- `.../joint_components.tsv`
- `.../joint_warmup_initial_values.tsv`
- 便携复核：`analysis/joint_pair_summary.csv`、`analysis/joint_selected60_fit_summary.csv`

**定量证据**

六组 best joint objectives：

- vi366：`18.8523`；
- vi322：`18.8901`；
- vi25：`18.9705`；
- vi311：`19.4145`；
- vi290：`19.7913`；
- vi138：`19.9782`。

全 60 个 selected solutions 的范围为 `18.8523–20.0192`。pair 间差异大于多数 pair 内 top-10 spread。

**支持结论**

joint optimizer 保留多个局部候选区域，当前结果不能证明唯一 global optimum。

**为何有科学意义**

这是决定 effect size 是否可信的核心优化诊断。可以保留跨所有 selected solutions 的方向级结果，但不能把精确 ratio 当作唯一生物效应量。

---

### 结果 10：*in vivo* effective proliferative ceiling 较低，是 60/60 selected joint solutions 中最稳定的方向之一

**结果位置**

- `top10/.../<pair>/seed*/joint_soft_coupling.tsv`
- `.../invitro_effective_params.tsv`
- `analysis/joint_soft_coupling_parameter_summary.csv`
- `analysis/joint_soft_coupling_all60.csv`

**定量证据**

- `lam_max` 的 *in vivo/in vitro* ratio：`0.099–0.222`；
- 中位数：`0.177`；
- 60/60 均为 *in vivo* 较低。

**支持结论**

在当前参数化下，tumor context 需要更低的 effective maximum proliferation rate。

**为何有科学意义**

这是 joint analysis 中跨所有 selected solutions 方向最稳定的 context axis。其解释应是综合 growth constraint，不能唯一归因于 oxygen。

---

### 结果 11：stress-linked missegregation response 在 *in vivo* 更强，且函数级比较支持该方向

**结果位置**

- `top10/.../joint_soft_coupling.tsv`
- Repo 的 `p_mis(N,O2)` 函数；
- `analysis/joint_missegregation_function_all60.csv`
- `analysis/joint_missegregation_function_comparison_all60.csv`

**定量证据**

- `p_misseg` ratio 中位数 `16.84`，范围 `11.12–47.38`；
- `k_o_mis` ratio 中位数 `17.11`；
- 将完整 response function 在 `O2=0.1%` 代入：
  - N=44 时，*in vivo/in vitro* per-chromosome missegregation ratio 中位数 `18.76`；
  - N=88 时，中位数 `19.14`。

**支持结论**

相较 culture，selected tumor solutions 表现出更强的 nonzero-stress missegregation response。baseline missegregation 本身不够稳定，不应作为 headline conclusion。

**为何有科学意义**

它直接支持“tumor environment 改变 CIN generation”这一主机制轴；函数级验证避免仅凭某个参数方向做错误解释。

---

### 结果 12：altered-daughter survival 的正确方向与原稿相反

**结果位置**

- `top10/.../joint_soft_coupling.tsv`
- Repo survival function：`s_N=s_max exp[-beta_buf (44/N)^{n_exp}]`；
- `analysis/joint_survival_function_all60.csv`
- `analysis/joint_survival_function_summary_all60.csv`

**定量证据**

60 个 selected joint solutions 的函数级中位数：

- *in vivo*: `s44=0.807`, `s88=0.911`, `s88/s44=1.12`, absolute gradient `0.096`；
- *in vitro*: `s44=0.204`, `s88=0.837`, `s88/s44=4.11`, absolute gradient `0.633`；
- 60/60 solutions 中，*in vivo* 在 N=44 与 N=88 的 absolute survival 都更高；
- 0/60 solutions 中，*in vivo* ploidy gradient 比 *in vitro* 更强。

**支持结论**

正确结论是：**altered daughters have higher absolute fitted survival in vivo, but survival is less strongly dependent on mother-cell ploidy in vivo; culture has the stronger ploidy-dependent buffering gradient.**

**为何有科学意义**

这是本次重新审查最重要的科学方向更正。原稿根据单个 buffering parameters 推断 “more ploidy-dependent in vivo” 会给出相反的功能解释。非线性模型必须比较完整 derived function。

---

### 结果 13：Welsch soft-coupling 多数已接近饱和，不能把“大差异仍存在”解释为强 shrinkage 下的精确识别

**结果位置**

- Repo joint backend 与 README 中 Welsch 公式；
- `top10/.../joint_soft_coupling.tsv`
- `analysis/joint_soft_coupling_parameter_summary.csv`
- `analysis/joint_selected60_fit_summary.csv`

**定量证据**

- `sigma=0.65`, `c=0.4`；
- 单参数 penalty cap：`0.5*c^2=0.08`；
- 14 参数总 penalty 已占其理论 cap 的 `75.8–89.5%`，中位数 `85.4%`；
- `lam_max`, `p_misseg`, `k_o_mis`, `buffer_beta`, `gamma_mu` 等 headline splits 在全部或近全部 selected solutions 中达到/接近 saturation；
- 60/60 local refinements 被尝试，但 0/60 accepted。

**支持结论**

selected ratios 是多个局部候选解中的方向性结果。Welsch saturation 后继续增大 context separation 几乎不再增加 penalty，因此不能说这些差异是在持续强 Gaussian shrinkage 下被精确识别出来的。

**为何有科学意义**

它限制了所有精确 fold-change 的可信度，也解释了为何 warm-start region 可以保留显著不同的 context splits。

---

### 结果 14：full-500 fixed-O2 analysis 显示多种 response shapes，但多数解的 dominant mode 分离不足

**原始结果位置**

- `joint_pre.zip` 的 `fixed_o2_ploidy_monotonicity_curves.tsv`；
- `fixed_o2_ploidy_monotonicity_by_seed.tsv`；
- class counts、class crosswalk、run arguments 和 classifier validation；
- Repo `fixed_o2_ploidy_monotonicity.R` 与 `fix_o2_simulation.R`；
- 便携复核：`analysis/fixed_o2_by_seed_500.csv`、`analysis/fixed_o2_class_counts.csv`、`analysis/fixed_o2_reliability_counts.csv`。

**定量证据**

- `500 seeds × 201 O2 values = 100,500` analytical attractor rows；
- raw curve classes：
  - complex nonmonotone `261`（52.2%）；
  - inverted-U `149`（29.8%）；
  - monotone increasing `57`（11.4%）；
  - U-shaped `25`（5.0%）；
  - increase-then-plateau `5`（1.0%）；
  - monotone decreasing `3`（0.6%）；
- spectral-gap reliability：
  - reliable `141`（28.2%）；
  - caution `168`（33.6%）；
  - unreliable `191`（38.2%）；
- caution + unreliable 合计 `71.8%`；
- 6 个 joint *in vivo* warm starts 中 0 个 reliable；3 个 caution，3 个 unreliable。

**支持结论**

多种 oxygen–ploidy response regimes 的确存在于 fitted solution set；但大多数解的 dominant eigenmode 没有足够清晰的 spectral separation。曲线类别应被视为 conditional asymptotic diagnostics，不是已验证生物亚型。

**为何有科学意义**

`joint_pre.zip` 使 Figure 6 的 class counts 和 spectral-gap caveat 可以独立复核，也显著降低了对六个 joint warm starts 作为稳定长时机制代表的信心。

---

### 结果 15：Figure 6 objective overlay 使用的不是 fitted MAP objective

**结果位置**

- `joint_pre.zip` 的 fixed-O2 by-seed table 与 objective-rank table；
- Repo `o2ipa_choose_objective(..., objective_source="auto")`；
- `analysis/fixed_o2_objective_definition_audit.csv`
- `analysis/fixed_o2_objective_spearman_correlations.csv`

**定量证据**

- 500/500 rows 的 selected `objective_source` 为 `raw_likelihood`；
- 该 raw diagnostic 是 burden+ploidy raw likelihood，未包含完整 MAP 中相同尺度的 necrosis、prior 和 normalized data structure；
- raw likelihood 与 separate-fit MAP objective 的 Spearman correlation：`0.429`；
- `objective_data` 与 MAP 的 Spearman correlation：`0.915`。

**支持结论**

Figure 6 的 objective overlay 不能被称为“best MAP fits”或用于选择 preferred MAP response class，除非明确重算或重新标注 objective definition。

**为何有科学意义**

错误的 objective label 会改变读者对 response class 与拟合质量关系的判断，是结果图中的实质性方法问题。

---

### 结果 16：primary parameter clusters 与 fixed-O2 response class 仅中度相关

**结果位置**

- `joint_pre.zip` 的 cluster assignments 与 fixed-O2 class assignments；
- `analysis/fixed_o2_cluster_class_association.csv`
- `analysis/fixed_o2_landscape_cluster_crosswalk.csv`

**定量证据**

- Cramér's V：`0.290`。

**支持结论**

parameter-space cluster 与 long-term oxygen-response class 不能互换使用。一个 parameter cluster 内可以包含不同 response shapes，反之亦然。

**为何有科学意义**

这阻止了把 exploratory cluster labels 直接转化为“biological response types”。

---

## 4. 代码—结果一致性与发布阻断项

### 结果 17：necrosis likelihood 内部有效，但 export 损坏，且原稿 equation 与实现相差一倍

**结果位置**

- `top10/fit_invivo.../seed*/necrosis_fit.tsv`；
- joint seed directories 中的相同文件；
- `.../burden_fit.tsv` 的 harvest-time dead/total burden；
- `.../fit_summary.tsv` 的 necrosis objective；
- Repo in-vivo backend 的 necrosis loss；
- 便携复核：`analysis/invivo_top10_fit_quality.csv`、`analysis_tables/invivo_necrosis_reconstruction.csv`。

**定量证据**

- 10/10 independent *in vivo* 与 60/60 joint `necrosis_fit.tsv` 中 predicted fields 为 `NaN`；
- 用 `burden_fit.tsv` 的 terminal dead/total ratio、相同 clipping/logit 和 `sigma=0.75` 可精确重构 reported objective；
- 无 `1/2` 的重构 objective：`0.631–1.118`，与 `fit_summary.tsv` 一致；
- 加 `1/2` 后为 `0.316–0.559`；
- 原稿 equation 含 `1/2`，实现不含，存在 factor-of-two mismatch。

**支持结论**

necrosis term 确实参与拟合；问题位于方法公式与 post-fit export。修订手稿已将 equation 对齐实际实现，但结果导出必须修复。

**为何有科学意义**

这是可审计性与 objective definition 的高优先级问题。若不修复，读者不能从标准导出文件追踪逐样本 necrosis residual。

---

### 结果 18：Figure 4 AUC exact ordering 仍无法从供应文件独立重算

**结果位置**

- Repo 图与 AUC 分析脚本；
- `joint_pre.zip` 文件清单；
- `analysis/joint_pre_AUC_file_inventory.txt`。

**证据**

`joint_pre.zip` 包含 fixed-O2 curves/classes、cluster tables 和多种 crosswalk，但没有发现可明确对应 Figure 4 B–D 的 seed-level AUC input/result table。

**支持结论**

手稿可以把这些 panels 描述为 descriptive one-feature screens；不能声称本次交付已独立复核其精确 feature ranking。

**发布要求**

公开每个 reference O2、outcome definition、feature matrix、seed IDs、AUC values、方向和生成命令。

---

### 结果 19：结果 provenance 尚不足以精确重建最终图表

**结果位置**

- selected *in vivo* 与 joint `run_provenance.tsv`；
- pair-level provenance logs；
- *in vitro* selected directories 的文件清单；
- `analysis/provenance_audit.csv`。

**证据**

- selected independent *in vivo* runs 记录 commit `54c5d500...` 且 working tree 为 dirty；
- selected joint runs 记录 commit `3d498599...` 且为 dirty；
- selected *in vitro* seed directories 缺少同等粒度 `run_provenance.tsv`；
- independent 与 joint 来源于不同 commit；
- pair-level provenance 是多任务追加日志，不能单凭一个 commit 精确恢复所有 task 的工作树。

**支持结论**

结果具有内部分析价值，但正式发布前必须从 clean tagged commit 重新生成或封存 exact dirty diff。

**为何有科学意义**

模型结论依赖代码实现细节；无法精确恢复代码状态会削弱所有 figure-to-result traceability。

---

## 5. 最终证据排序

### 最稳健、可作为主结论

1. 4N severe-deprivation chromosome decline 在 *in vitro* top-10 中稳定；
2. 该下降不需要大量 direct hypoxia-origin killing；
3. selected joint solutions 中 *in vivo* effective proliferative ceiling 一致较低；
4. selected joint solutions 中 nonzero-stress missegregation response 一致更强；
5. 完整 survival function 显示 *in vivo* absolute survival higher、但 ploidy dependence weaker。

### 可作为探索性/条件性结论

1. 500-seed fixed-O2 solution set 包含多种 response shapes；
2. parameter landscape 存在 exploratory clusters；
3. *in vivo* 资源轴可用于模型内部排序；
4. one-feature AUC 可用于候选 driver screening。

### 当前不应作为确定结论

1. 唯一的 oxygen–ploidy response class；
2. 六个 warm starts 代表六个真实生物亚型；
3. 精确的 joint fold-change 是唯一可识别 effect size；
4. t-SNE cluster 间距离、面积或视觉分离代表真实生物距离；
5. CIN-associated dead fraction 等于实验细胞死亡；
6. fitted effective O2 等于直接组织 pO2；
7. Figure 4 AUC 的精确 feature ranking 已被本次复核；
8. Figure 6 的 raw-likelihood objective overlay 等于 MAP fit ranking。
