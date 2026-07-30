# 手稿扩展与实质性更正说明（纳入 joint_pre.zip）

## 总体原则

完成版手稿保留原稿的主问题与总体结构：资源限制一方面对高染色体状态施加 proliferation/death cost，另一方面可能增加 CIN/WGD variation；最终 ploidy trajectory 由 variation generation 与 altered-daughter filtering 的平衡决定。没有增加新实验、外部队列或新的模型模块。

本次扩展全部来自：

- `top10.zip` 中现有 separate/joint fitting outputs；
- `joint_pre.zip` 中已完成的 parameter-landscape 与 fixed-O2 preprocessing；
- Repo 当前代码对 output semantics、objective、soft coupling 和 fixed-O2 analysis 的定义；
- 对现有 nonlinear response functions 和 likelihood components 的直接重算。

以下逐项列出相对 Repo 原稿新增或实质改写的内容、原因与意义。

---

## 1. 增加完整 Abstract

**增加内容**

新增研究背景、数据流、模型结构、主要结果、corrected survival interpretation、full-500 reliability 和 joint uncertainty 的结构化摘要。

**为何增加**

Repo 原稿没有完整 Abstract，且正文中存在未收束的结果方向。

**意义**

形成可投稿的完整稿件入口，并在摘要中提前限制 latent oxygen、warm-start dependence 和 identifiability 的解释。

---

## 2. 清理 Introduction 中所有编辑占位与未完成理论提纲

**删除/整合内容**

- `XX:` sister-paper/CLONEID 占位；
- 未处理 PubMed 链接提醒；
- mitochondrial adaptation 编辑备注；
- mutator theory、stress-induced mutagenesis、fitness-valley 的提纲式草稿。

**为何修改**

这些内容不是完成的结果论证，部分需要新文献和新理论扩展，违反“不要添加新内容”的要求。

**意义**

Introduction 聚焦当前模型实际表示并由结果支持的机制。

---

## 3. 将 oxygen 明确限定为 *in vivo* latent effective resource coordinate

**修改内容**

在 Introduction、Methods、Results 和 Discussion 中统一说明：

- *in vitro* oxygen 是 assigned condition；
- *in vivo* effective O2 是 mean-field latent model variable；
- 它可吸收 perfusion、nutrient competition、density、waste、stroma 等未建模约束；
- 不能与组织实测 pO2 一一对应。

**为何修改**

原稿部分段落容易把 fitted O2 当作直接氧测量或把所有 context difference 唯一归因于 oxygen。

**意义**

避免过度因果解释，并正确表达模型变量的可观测性层级。

---

## 4. 增加 top-10 robustness，而不把 top-10 当作置信区间

**增加内容**

加入 separate *in vitro* 与 *in vivo* top-10 的 objective、trajectory、fit-quality 与 dead-compartment ranges；加入 60 个 selected joint solutions 的方向稳定性。

**为何增加**

原稿主要引用单一 seed 的精确数值。top10.zip 的核心用途是判断结论是否依赖单个 optimizer seed。

**意义**

把结论表述从“best-fit example”提升为“selected near-optimal solutions 中方向稳定”，同时明确这些范围不构成统计置信区间。

---

## 5. 更正 *in vitro* 4N direct-death 数值与机制语言

**原稿问题**

原稿以单一解写 direct hypoxia dead `at or below 1.7%`，并容易把较大的 dead-buffer component 解释为纯生物死亡。

**修订**

- 改为 top-10 范围 `1.74–1.82%`；
- 增加 CIN-associated dead fraction `32.15–33.40%`；
- 明确其含 nonviable daughters 与 finite-grid boundary routing；
- 保留“ploidy reshaping does not require strong direct high-ploidy elimination”。

**意义**

保留最重要机制结论，同时消除不准确绝对上限与 compartment overinterpretation。

---

## 6. 实质改写 *in vitro* 2N deprived trajectory

**原稿问题**

原稿称 2N deprived lineage 经过 transient high-ploidy/WGD expansion 后整体回到较低 ploidy。

**修订**

- 改为形成 high-chromosome component，同时保留 lower-chromosome mode；
- 明确 top-10 late prediction `63.72–64.40`；
- 明确两个观测分支均值 `66.85` 与 `88.05`；
- 明确当前 shared-lineage model 未重现 late branch divergence；
- 删除总体 trajectory 已经回落的强结论。

**意义**

将“分布中存在 lower-chromosome descendants”与“总体 mean depolyploidizes”严格区分。

---

## 7. 增加 *in vivo* 绝对拟合质量与 cohort-level limitation

**增加内容**

报告 burden RMSE、terminal mean-N MAE、Wasserstein-1、total variation 和 necrosis MAE，并说明相同 cohort/harvest-time 小鼠共享 prediction。

**为何增加**

原稿只强调机制和目标函数，没有充分展示 individual trajectory/distribution mismatch。

**意义**

把模型准确定位为 cohort-level mechanism inference，而不是 individual-tumor predictor。

---

## 8. 增加 practical non-identifiability 与 boundary dependence

**增加内容**

描述近似 objective 对应宽参数范围、多个 bound hits，以及 10/10 selected *in vivo* fits 中 `sigma_burden` 到达上界。

**为何增加**

这是理解 multiple response classes 与 joint local regions 的核心结果，原稿不足。

**意义**

阻止把 individual parameters 解释为唯一精确生物测量，并为 profile-likelihood/functional uncertainty 的审稿要求提供依据。

---

## 9. 更正 joint parameter architecture：14 个 biological parameters 全部 soft-coupled

**原稿问题**

Supplementary Methods 将 `alpha_o2` 与 `gamma_growth` 写为 hard-shared，并保留 “I disagree” 编辑备注；Results 又写 14 个参数均 context-specific。

**修订**

- 统一为 14 个参数的 transformed center–delta representation；
- 明确本次 results 中没有 active hard-shared biological parameters；
- 明确 union transformed bounds、feasibility check 和 no clipping。

**意义**

消除 Methods、Results、config、README 与实际 `joint_soft_coupling.tsv` 之间的矛盾。

---

## 10. 更正 soft-coupling penalty：quadratic 改为 bounded Welsch

**原稿问题**

原稿 equation 使用 Gaussian/quadratic penalty。

**修订**

加入实际 Welsch 公式：

`0.5*c^2*(1-exp[-(|delta|/(sigma*c))^2])`

并说明 `sigma=0.65`, `c=0.4`, 单参数 cap `0.08`。

**意义**

Welsch saturation 会实质改变 large context splits 的解释。大分离不再持续受到二次惩罚，不能使用“强 shrinkage 下仍分开”的旧语言。

---

## 11. 增加 joint optimization 与 warm-start coverage 诊断

**增加内容**

- 六组 best objectives `18.852–19.978`；
- 60/60 local refinement attempted，0/60 accepted；
- total soft penalty 占 cap 的 `75.8–89.5%`；
- 六个 *in vivo* objective-minimum representatives；
- 六组全部绑定同一 *in vitro* seed 10；
- 未覆盖其他 *in vitro* subclusters。

**为何增加**

这些信息决定 joint solutions 是否可被称为 global/representative mechanisms。

**意义**

将精确 effect-size 结论降级为 local candidate regimes 中稳定的方向性 hypothesis。

---

## 12. 更正 altered-daughter survival 的科学方向

**原稿问题**

原稿根据 `s_max`, `beta_buf`, `n_exp` 的单独方向推断 *in vivo* survival “more ploidy-dependent”。

**修订**

直接计算完整 nonlinear function：

- *in vivo*: `s44=0.807`, `s88=0.911`, ratio `1.12`, gradient `0.096`；
- *in vitro*: `s44=0.204`, `s88=0.837`, ratio `4.11`, gradient `0.633`。

完成版在 Abstract、Results、Figure 5 caption 和 Discussion 中统一改为：

**absolute altered-daughter survival is higher in vivo, but ploidy dependence is weaker in vivo; culture has stronger ploidy-dependent buffering.**

**意义**

这是最重要的结果方向更正，避免非线性参数误读。

---

## 13. 将 missegregation comparison 从单参数提升为完整函数比较

**增加内容**

除 `p_misseg` 与 `k_o_mis` ratio 外，将 Eq. `p_mis(N,O2)` 在 `O2=0.1%`、N=44/88 代入，得到约 `18.8–19.1×` 的 median *in vivo/in vitro* ratio。

**为何增加**

`p_misseg` 与 `k_o_mis` 同时变化，单看一个参数不能完整表征 stress response。

**意义**

函数级方向比单参数 fold-change 更接近模型实际生物预测。

---

## 14. 增加 joint_pre 的 parameter-landscape preprocessing 方法与限制

**增加内容**

- 500 seeds/context；
- 14 parameters standardized；
- pooled initial/best t-SNE；
- primary cluster sizes 与 silhouette；
- subclustering 在 standardized original features；
- *in vivo* subcluster silhouettes `0.179/0.096/0.330`；
- representatives 是 objective minima，不是 medoids；
- t-SNE clusters 是 exploratory strata，不是 biological subtype validation。

**为何增加**

joint warm starts 的来源以前只在 figure caption 中简写，缺少可复核方法与 representativeness caveat。

**意义**

读者可以理解六个 joint pairings 的采样逻辑和覆盖不足。

---

## 15. 增加 full-500 fixed-O2 classification 与 spectral-gap reliability

**增加内容**

- 500 × 201 = 100,500 analytical attractor evaluations；
- `shape_first_v3` finite-difference classification；
- 六类 raw curve counts；
- reliable/caution/unreliable counts `141/168/191`；
- 71.8% caution 或 unreliable；
- 六个 joint warm starts 中没有 reliable seed。

**为何增加**

`joint_pre.zip` 首次提供足够表格独立复核 Figure 6 class counts 和 reliability。

**意义**

它把“multiple regimes”从视觉描述提升为数值结果，同时要求把大多数 small-gap curves 的机制解释降级。

---

## 16. 更正 fixed-O2 class assignment 方法

**原稿问题**

原稿写曲线由 “regression-smoothed response shape” 分组。

**修订**

改为实际记录的 `shape_first_v3` finite-difference rule，包含 signed steps、reversal fraction、amplitude 和 terminal plateau。

**意义**

确保 Methods 与实际代码一致，避免读者误以为使用回归平滑模型或其统计不确定性。

---

## 17. 增加 Figure 6 objective-definition audit

**增加内容**

- fixed-O2 pipeline 的 `objective_source=raw_likelihood`（500/500）；
- raw likelihood 与 MAP Spearman `rho=0.429`；
- `objective_data` 与 MAP `rho=0.915`；
- Figure 6 objective overlay 需要明确重算或 relabel。

**为何增加**

原稿把该 overlay 表述为 fit objective/high-quality fits，实际并非完整 MAP objective。

**意义**

防止使用错误 ranking 得出 preferred response class。

---

## 18. 保留 Figure 4 AUC 为 descriptive screen，但明确无法独立重算 exact ranking

**修改内容**

Results 与 caption 不再把 one-feature AUC 当作 causal driver proof；新增说明 joint_pre/top10 未含明确的 seed-level AUC input/result tables。

**为何修改**

本次可复核 fixed-O2/cluster outputs，但仍缺 Figure 4 B–D exact AUC tables。

**意义**

区分“分析设计可见”与“精确数值已复核”。

---

## 19. 更正 necrosis loss equation 并新增 export audit

**原稿问题**

原稿 equation 含 `1/2`，但 Repo 实现和 `fit_summary.tsv` 的 objective 不含；所有 selected `necrosis_fit.tsv` predicted fields 为 `NaN`。

**修订**

- equation 改为实际实现的 mean standardized squared logit residual，无额外 `1/2`；
- 增加 reconstructed no-half objective `0.631–1.118` 与 reported value 一致；
- 明确 export bug 与 likelihood active 是两个不同问题。

**意义**

使 Methods 与 fitted objective 一致，并指出必须修复的可审计性问题。

---

## 20. 增加 finite-grid boundary sensitivity limitation

**增加内容**

明确 out-of-grid daughter/WGD mass 被路由到 CIN-associated dead compartment，并要求在将该 compartment 生物学化之前做 grid sensitivity 或分离输出。

**意义**

防止 numerical truncation 被误解释为 biological lethality。

---

## 21. 增加 provenance 与 reproducibility section

**增加内容**

- selected independent *in vivo* 与 joint 来自不同 commits；
- 两组均记录 dirty worktree；
- selected *in vitro* 缺同等 provenance file；
- final release 需要 clean commit、resolved YAML、parameter table、seed manifest、environment/container、figure-to-source map 和 checksums。

**意义**

将 computational provenance 变成正式发布标准，而不是隐藏在运行目录中的工程细节。

---

## 22. Figure captions 全面改为与证据一致

**修改内容**

- Figure 3：加入 top-10 数值、branch-divergence limitation 和 CIN/boundary caveat；
- Figure 4：AUC 改为 descriptive screen，并说明 exact tables 缺失；
- Figure 5：更正 survival direction、14 soft-coupled Welsch、warm-start reliability 与 coverage；
- Figure 6：改为 finite-difference class rule、spectral-gap counts、raw-likelihood objective label。

**意义**

确保快速阅读图注不会得到与正文相反或超出证据的结论。
