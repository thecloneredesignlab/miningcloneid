# 代码—结果—手稿一致性验证报告

## Overall assessment: Needs revision before external sharing

该分析包包含足够信息支持若干方向稳定的机制结论，并且 `joint_pre.zip` 补齐了 full-500 fixed-O2 和 parameter-landscape 的主要中间表。然而，joint optimization coverage、practical identifiability、spectral-gap reliability、objective definition、necrosis export 和 provenance 仍存在会影响论文结论或可复现性的高优先级问题。

---

## 1. 验证对象

- Repo：`thecloneredesignlab/miningcloneid`, branch `HypoxiaLTEEFigures`；
- manuscript：`manuscript/ltee_hypoxia_model.tex`；
- result package：`top10.zip`；
- joint preprocessing：`joint_pre.zip`；
- independent reanalysis：`scripts/reanalyse_results.py` 及 `analysis/` 输出；
- revised manuscript：`03_revised_manuscript.tex`。

---

## 2. 数据与结果目录完整性

### Separate *in vitro*

10/10 selected seeds 包含：

- `fit_summary.tsv`；
- natural/transformed parameters；
- daily counts；
- lineage/distribution summaries；
- growth/karyotype/flow likelihood outputs；
- HTML report。

### Separate *in vivo*

10/10 selected seeds 包含：

- `fit_summary.tsv`；
- `burden_fit.tsv`；
- `terminal_ploidy_fit.tsv`；
- `necrosis_fit.tsv`；
- parameter/config/provenance files；
- HTML report。

### Joint

6 pairs × 10 selected seeds 包含：

- `fit_summary.tsv`；
- `joint_components.tsv`；
- `joint_soft_coupling.tsv`；
- effective context parameter outputs；
- warm-start and initialization records；
- burden/ploidy/necrosis outputs；
- HTML report。

### joint_pre

包含：

- 500-seed landscape inputs/assignments；
- primary/subcluster summaries and silhouettes；
- representative-selection crosswalks；
- full-500 fixed-O2 curves and per-seed classification；
- class counts、reliability、objective-rank tables；
- validation and run-argument tables。

未发现明确对应 Figure 4 B–D 的 seed-level AUC input/result table。

---

## 3. 关键计算 spot checks

### 3.1 *In vitro* 4N severe-deprivation metrics

**Status: Verified**

直接从 lineage/daily-count outputs 重算：

- decline `3.54–4.17` chromosomes；
- max hypoxia-dead `1.74–1.82%`；
- max CIN-dead `32.15–33.40%`。

### 3.2 *In vivo* burden fit

**Status: Verified**

对 included positive observations 在 log scale 重算 RMSE，并按 tumor mean-squared-error equal weighting 重算 tumor-balanced RMSE。结果与便携复核表一致。

### 3.3 Terminal chromosome-number fit

**Status: Verified**

从 per-state terminal outputs 重算 sample mean error、1D Wasserstein-1 和 total variation。结果显示 mean-level fit 中等，但 distribution shape mismatch 较大。

### 3.4 Necrosis likelihood

**Status: Partially verified; export defect found**

- `necrosis_fit.tsv` predicted fields 在所有 selected independent/joint results 中为 missing；
- 从 terminal burden row 的 total-dead/total ratio 重构 predicted necrosis；
- clipping、logit、`sigma=0.75` 后，不含 `1/2` 的 mean standardized squared residual 精确匹配 reported objective；
- 原稿含 `1/2` 的 equation 比实现小一倍。

结论：likelihood active；standard export broken；original equation wrong。

### 3.5 Joint objective decomposition

**Status: Verified at selected-solution level**

`objective_total` 可由 weighted *in vivo*、weighted *in vitro*、soft penalty 和 constraint terms 对齐。current joint restrictions 为 inactive/zero contribution。

### 3.6 Soft-coupling semantics

**Status: Verified**

结果与 Repo 当前实现一致：14 parameters soft-coupled，center–delta reconstruction，union bounds，Welsch penalty。原稿 hard-shared/quadratic description 不一致。

### 3.7 Survival function direction

**Status: Verified and original narrative contradicted**

完整 nonlinear function 显示：

- higher absolute survival in vivo at N=44/88；
- much weaker ploidy gradient in vivo；
- stronger ploidy-dependent buffering in vitro。

### 3.8 Fixed-O2 class/reliability counts

**Status: Verified from supplied full tables**

500 seeds、201 O2 values、100,500 rows；class counts 与 reliability counts 可独立重算。

### 3.9 Fixed-O2 objective ranking

**Status: Discrepancy found**

pipeline's automatic selected objective 是 raw burden+ploidy likelihood，不是 full MAP。raw-vs-MAP Spearman `0.429`，因此 Figure 6 objective overlay 的 label/interpretation 不成立。

---

## 4. 方法与统计风险

### High: Practical identifiability

近优 objective 对应宽参数范围与多个 bound hits。当前仅依赖 top-ranked multistarts，不足以为 individual parameters 建立 likelihood-based uncertainty。

**Required fix**

对 headline parameters 和 derived functions 进行 profile likelihood、bootstrap、bounded multistart envelope 或等价 practical-identifiability analysis；优先报告 function-level uncertainty。

### High: Joint warm-start coverage

六个 *in vivo* objective-minimum representatives 均配同一个 *in vitro* seed，未覆盖其他 culture subclusters。

**Required fix**

至少进行 representative *in vivo* × representative *in vitro* cross-design，或证明扩展 starts 不改变方向/函数区间。

### High: Local convergence/global optimum

60/60 selected joint local refinements 均未 accepted；pair best objective 差异明显。

**Required fix**

诊断 convergence code、gradient/finite-difference behavior、bounds、initialization；增加 independent starts 和 alternative penalty scales；不能使用 global optimum language。

### High: Fixed-O2 spectral reliability

只有 28.2% separate fits reliable，六个 joint warm starts 无一 reliable。

**Required fix**

主图与正文必须区分 reliable/caution/unreliable；对 warm starts 报告 relaxation/spectral-gap uncertainty；避免把 small-gap class 当作稳定机制。

### High: Necrosis export/equation

标准表 predicted values missing，原 equation 与 implementation 相差一倍。

**Required fix**

修复 export，加入 per-seed automated reconstruction test，并统一 code/report/manuscript equation。

### Medium–High: Welsch saturation

大多数 total penalty 接近 cap，多个 headline splits 饱和。

**Required fix**

报告 standardized deltas、penalty saturation 和 alternative `sigma/c` sensitivity；不要把 exact fold-change 解释为强 regularization 下的唯一估计。

### Medium–High: Finite-grid boundary routing

CIN-dead 混合 biological nonviability 与 numerical boundary loss。

**Required fix**

单独输出 `misseg_nonviable` 与 `boundary_dropped`，并对 N_MIN/N_MAX 做 grid sensitivity。

### Medium: t-SNE cluster interpretation

primary silhouette 是 embedding-space measure；subclusters 中多个 silhouette 很低。

**Required fix**

将 clusters 定义为 exploratory strata；报告 t-SNE hyperparameters/seeds/stability；在 original transformed feature space 验证 cluster stability，不使用 cluster area/distance 作生物解释。

### Medium: Cohort/replicate heterogeneity

mouse-specific 和 late culture branch-specific effects 未表示。

**Required fix**

明确 limitation；若要解释 individual heterogeneity，需要 hierarchical/random-effect 或 branch-specific latent model。

---

## 5. Visualization review

### Figure 3

原 caption 对 2N trajectory 的“随后回落”过强；应显示 mixed distribution 和 branch mismatch。4N dead decomposition 应把 direct hypoxia loss 与 CIN/boundary loss分开。

### Figure 4

AUC panels 缺 seed-level source tables。应在 supplement 发布 exact inputs/results，并将 AUC 定义为 one-feature discrimination，不是 causal importance。

### Figure 5

survival panel/narrative 必须按函数级结果重画或重标：absolute survival higher in vivo, gradient weaker in vivo。Caption 应说明 one-sided *in vitro* warm-start coverage 和 Welsch saturation。

### Figure 6

- curve classification 应标 `shape_first_v3` finite-difference rule；
- reliable/caution/unreliable 应可视化；
- objective overlay 必须重算为 MAP 或明确标为 raw likelihood；
- no warm-start is reliable 应在 joint-selection context 中可见。

---

## 6. Provenance and reproducibility

**High-severity finding**

- independent *in vivo* 与 joint selected runs 来自不同 commits；
- 两者均记录 dirty working tree；
- selected *in vitro* 缺同等 provenance file；
- Figure 4 AUC intermediate tables 未供应。

**Required release bundle**

1. clean tagged commit 或 exact dirty patch；
2. resolved YAML/config and parameter tables；
3. all seed IDs and exact selection rules；
4. software/package versions or container digest；
5. raw inputs and checksums；
6. full intermediate tables underlying every panel；
7. figure-to-source table；
8. one command/workflow that regenerates all figures and manuscript tables。

---

## 7. Readiness assessment

### Directionally supported

- low direct-death requirement for 4N decline；
- lower effective proliferative ceiling in vivo；
- stronger nonzero-stress missegregation response in vivo；
- corrected absolute-survival/gradient relationship。

### Requires caveated presentation

- all precise joint fold changes；
- parameter clusters and response classes；
- latent O2 trajectories；
- CIN-associated dead fraction；
- individual-tumor fit claims。

### Release blockers

1. survival direction/figure consistency；
2. joint method description；
3. fixed-O2 objective label；
4. spectral-gap reliability disclosure；
5. necrosis export/equation；
6. identifiability and optimization uncertainty；
7. clean provenance and missing AUC tables。
