# 修订版手稿相对 Repo 原稿的完整差异记录

## 1. 对照对象

- 原稿：Repo 分支 `HypoxiaLTEEFigures` 中 `manuscript/ltee_hypoxia_model.tex`；
- 修订稿：本交付包 `03_revised_manuscript.tex`；
- 对照原则：记录语义、方法、数值、结构和结论方向的全部实质性变化，而不是逐字符 patch。未改变的 LaTeX 图文件、参数表与 bibliography 路径予以保留。

---

## 2. Front matter

### 标题

保留原题：

`Resource limitation rewires chromosome instability and ploidy evolution across in vitro and in vivo cancer models`

删除原稿中的 alternative-title 编辑备注。

### 作者与单位

保留作者顺序与共同第一作者/通讯作者结构；统一 Moffitt 机构名称和标点。

### Abstract

原稿没有完成的摘要。修订稿新增完整 Abstract，覆盖：

- WGD/aneuploidy 研究问题；
- matched SUM159 2N/4N *in vitro*/*in vivo* design；
- five mechanistic components；
- 14-parameter Welsch soft coupling；
- *in vitro* 4N decline 与 direct-death range；
- joint `lam_max` 与 missegregation direction；
- corrected survival-function direction；
- full-500 spectral-gap reliability；
- warm-start/local-refinement uncertainty。

---

## 3. Introduction

### 删除内容

删除全部未完成编辑占位：

- sister paper/CLONEID 的 `XX` 句；
- PubMed link reminder；
- mitochondrial adaptation reminder；
- Mutator Theory/Stress-Induced Mutagenesis/fitness-valley 的提纲式文本；
- 未整合、未完成的 bold summary 草稿格式。

### 重组内容

原稿较长的文献综述被重组为五个连续逻辑段：

1. context-dependent WGD/aneuploidy prevalence；
2. high chromosome content 的 energetic/homeostatic cost；
3. resource stress 同时促进 CIN variation；
4. oxygen 作为 tractable but incomplete model variable；
5. generation–filtering tension 与本文研究问题。

### 新增限制

明确 *in vivo* fitted oxygen 是 latent effective resource coordinate，不是直接 pO2 measurement。

---

## 4. Main Methods

原稿 main Methods 只有概要段，详细方法分散在 Supplementary Information。修订稿将投稿所需的核心可复现信息前移到 main Methods，同时保留 supplementary detail。

### Study design and data streams

保留八个 untreated tumors、matched 2N/4N lineages 和 *in vitro* passage data；增加：

- terminal histologic necrosis 的条件性使用；
- paired inclusion；
- joint 六个 *in vivo* warm starts 全部绑定一个 *in vitro* seed 的事实。

### State representation

保留 viable chromosome-state distribution 和两类 dead pools；增加：

- `N_MIN=22`, `N_MAX=154`；
- burden 与 viable single-cell observations 的不同语义；
- CIN-dead 中包含 boundary routing。

### Resource response equations

保留原始 Hill stress、growth、death、missegregation 和 survival 函数，但压缩重复说明并统一 notation。

### In vivo oxygen dynamics

保留 mean-field supply-demand model；增加当前 fit flags：

- untreated；
- dynamic oxygen feedback；
- chromosome-dependent growth；
- no harvest-specific multiplier；
- no crowding。

### Observation models

增加完整 objective hierarchy 与 inclusion semantics。特别更正：

- necrosis equation 删除原稿额外的 `1/2`，与代码和 reported objective 对齐；
- 明确 day-0/nonpositive burden exclusion；
- 明确 *in vitro* 三项 likelihood 是各自归一化后组合。

### Joint soft coupling

这是最重要的方法更正之一。

原稿：

- 将 `alpha_o2` 与 `gamma_growth` 写为 hard-shared；
- 其余 12 个参数 soft-coupled；
- 保留作者编辑备注 “I disagree”；
- penalty 写成 quadratic Gaussian。

修订稿：

- 14 个 active biological parameters 全部 soft-coupled；
- 无 active hard-shared biological parameter；
- `vivo=center+delta/2`, `vitro=center-delta/2`；
- union transformed bounds；
- infeasible reconstruction rejected, no clipping；
- 实际 Welsch penalty；
- `sigma=0.65`, `c=0.4`, per-parameter cap `0.08`。

### Optimization and robustness

新增 separate top-10 与 joint six-pair top-10 的用途说明；明确 selected sets 不是 posterior/confidence interval。

### Parameter-landscape preprocessing

原稿只在 Results/figure caption 中概括“three clusters/two subclusters”。修订稿新增：

- 500 seeds/context；
- 14-parameter standardization；
- pooled initial/best t-SNE；
- primary k-means cluster sizes 与 embedding-space silhouette；
- subclustering 在 original standardized features；
- subcluster silhouettes；
- representatives 是 objective minima 不是 medoids；
- all six share in-vitro seed10；
- exploratory coverage, not validated subtype。

### Fixed-O2 methods

原稿写 “regression-smoothed response shape”。修订稿更正为：

- 500 × 201 analytical dominant-eigenvector evaluations；
- `shape_first_v3` finite-difference rule；
- spectral-gap reliability thresholds；
- unreliable curves 统一给 final unreliable label；
- objective overlay 使用 raw likelihood 而非 MAP；
- AUC exact tables 缺失；
- projections 不是 independent interventions。

---

## 5. Results: matched design/model overview

保留原稿“resource limitation imposes opposing pressures”的核心逻辑；删除无法由结果直接支撑的宽泛因果语言。强调 matched lineage background 只减少 model-system confounding，并不把 context difference 自动变成单一环境变量的因果效应。

---

## 6. Results: *in vitro*

### 4N severe-deprivation

原稿：单一 seed，mean `84.3→80.8`，direct dead `≤1.7%`。

修订稿：

- top-10 objective `3.8525–3.8623`；
- decline `3.54–4.17`，median `3.86`；
- terminal `80.13–80.89`；
- direct hypoxia dead `1.74–1.82%`；
- CIN-associated dead `32.15–33.40%`；
- 明确 CIN compartment 不是纯实验 death。

结论保留但更精确：decline does not require strong direct elimination。

### 2N deprived lineage

原稿：描述为 transient high-ploidy expansion 后 shift toward lower ploidy。

修订稿：

- 改为 high-chromosome component 与 lower mode 共存；
- 明确 late predictions `63.72–64.40`；
- 明确 observations `66.85/88.05`；
- 明确 model does not reproduce late branch divergence；
- 删除整体 trajectory 已回落的强结论。

---

## 7. Results: separate *in vivo*

原稿主要讨论 latent oxygen、AUC 和 cluster。修订稿先增加 absolute fit-quality：

- objective range；
- all-observation 与 tumor-balanced burden RMSE；
- terminal mean-N MAE；
- Wasserstein-1；
- total variation；
- necrosis MAE。

新增 cohort-level prediction limitation、parameter non-uniqueness、bound hits 和 `sigma_burden` upper-bound issue。

原稿 AUC 的“strongest separators/drivers”语言降级为 descriptive one-feature associations；明确 exact seed-level AUC tables 未供应。

---

## 8. Results: joint context differences

### Proliferation

将原稿定性“lower in vivo”扩展为：

- 60/60 direction consistent；
- ratio range `0.099–0.222`；
- median `0.177`；
- 参数是 effective composite growth constraint。

### Missegregation

将原稿单参数解释扩展为：

- `p_misseg` median ratio `16.84`；
- `k_o_mis` median ratio `17.11`；
- full response at O2=0.1%, N=44/88：median `18.76/19.14`；
- baseline missegregation 不作为稳定 headline。

### Survival direction

原稿：more ploidy-dependent survival in vivo。

修订稿：完全反向更正为：

- absolute survival higher in vivo at N=44 and N=88；
- vivo `s88/s44=1.12`, vitro `4.11`；
- vivo absolute gradient `0.096`, vitro `0.633`；
- culture has stronger ploidy-dependent buffering。

### Optimization/warm-start

新增独立 subsection：

- pair best objectives；
- 0/60 local refinements accepted；
- Welsch cap occupancy；
- saturated headline splits；
- asymmetric warm-start coverage；
- no reliable warm-start fixed-O2 curve；
- no unique global optimum/effect size claim。

---

## 9. Results: full-500 fixed-O2

原稿：多个 response classes，class/objective overlap。

修订稿新增精确可复核数值：

- 100,500 evaluations；
- 261 complex、149 inverted-U、57 increasing、25 U、5 plateau、3 decreasing；
- 141 reliable、168 caution、191 unreliable；
- final interpretation counts；
- 71.8% caution/unreliable；
- cluster–class Cramér's V `0.290`；
- raw objective 与 MAP Spearman `0.429`；
- objective overlay 不选择 preferred MAP class。

---

## 10. Figure captions

### Figure 1

保留 matched data-stream overview；删除过度详细但未在图中可核对的表述。

### Figure 2

保留 opposing-pressure model；明确 fixed-O2 dominant eigenvector 是 long-term composition 的定义。

### Figure 3

加入 top-10 robustness、late branch mismatch、direct vs CIN dead ranges 和 boundary caveat。

### Figure 4

AUC 改为 descriptive screen；明确 exact AUC tables 未供应；clusters 标为 exploratory strata。

### Figure 5

- 14 soft-coupled Welsch；
- six objective-minimum *in vivo* representatives all paired with seed10；
- no reliable warm-start；
- survival direction corrected；
- ratios conditioned on warm starts/bounds/saturation。

### Figure 6

- regression-smoothed 改 finite-difference `shape_first_v3`；
- 加 141 reliable/191 unreliable；
- panel B 明确 projection not intervention；
- panel D 改为 raw-likelihood diagnostic，不是 MAP；
- 加 Spearman `0.429`。

---

## 11. Discussion

原稿主要将 joint parameters 写成三个 biological hypotheses。修订稿保留前两个方向，重写第三个，并增加五类限制：

1. nonlinear derived-function interpretation；
2. full-500 spectral-gap uncertainty；
3. low-silhouette subclusters 和 one-sided *in vitro* coverage；
4. practical non-identifiability、bounds、rejected local refinements 和 Welsch saturation；
5. cohort-level structure、finite-grid routing、latent oxygen；
6. necrosis export + factor-of-two mismatch；
7. fixed-O2 objective label 与 AUC tables；
8. dirty provenance。

治疗相关段落保留，但改为 model implication，不写成已验证 therapeutic prediction。

---

## 12. Supplementary Information

原稿 supplementary 包含大量完整方程，但也存在重复、typo、未完成实验段、错误 joint methods 和不一致 equation。修订稿保留足够重建模型的结构，新增并统一以下专题：

- *in vivo* experiment subset；
- *in vitro* lineage representation；
- dead compartments；
- finite-grid boundary handling；
- likelihood normalization/inclusion；
- landscape preprocessing/representative selection；
- full-500 fixed-O2 classification/reliability；
- joint initialization/bounds/penalty；
- optimization diagnostics；
- output audit/provenance；
- existing parameter-table inputs。

删除重复推导和与当前代码不一致的旧 hard-sharing/quadratic equations。

---

## 13. 保留的 Repo integration paths

修订稿继续使用：

- `oxygen/figures/assembled_fig1.png` 至 `assembled_fig6.png`；
- `oxygen/figures/iteration1/fig6d_selection_diagnostic.png`；
- `manuscript/tables/tao_model_parameters.tex`；
- `manuscript/tables/tao_fixed_parameters.tex`；
- `manuscript/references_Zotero_TaoLi.bib`。

因此，可将 `03_revised_manuscript.tex` 放入 Repo 后进行 production compile，但 Figure 5/6 内容本身需要按 corrected direction/objective label 重新生成或核对。
