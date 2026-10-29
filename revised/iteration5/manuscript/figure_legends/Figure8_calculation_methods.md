# Figure8 及相关补图：参数、有限时间模拟、传代与统计汇总方法

执行流程更新（2026-10-17）：当前完整入口为 `run_Figure8.R`，覆盖准备、模拟、数据生成、汇总、校验和绘图；数据使用统一无版本号目录 `data/Figures/Figure8`。见[完整流程说明](../../Code/Figures/Figure78_pipeline_ZH.md)。下文科学公式保持不变；其他章节中的旧 run 名和 helper 名仅追溯历史发布记录，不是新流程的运行依赖。正式数据尚未全量重新模拟，原数值算例仍标识历史参考结果。

编写及实现核对日期：2026-10-15（America/New_York）。

英文对应文件：[Figure8 calculation methods, English](/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration5/manuscript/figure_legends/Figure8_calculation_methods_EN.md)。中英文使用相同章节、公式、数值和复算代码；本文件记录当前实现，不提出新的模型或修改现有图。

最重要的统计层级是：单次模拟的群体平均倍性 → 同一拟合 seed 内的培养传代重复平均 → 多个拟合 seeds 之间的中值和 IQR。图上显示中值，并不意味着模型内部把群体平均倍性改成了细胞倍性的中位数。

## 1. 当前图形范围与版本选择

### 1.1 主图与补图对应表

下表中“原生拟合 p”表示每个 seed 自己的最终 `p_misseg`；“固定 p”表示只覆盖该原生参数，不固定有效错误分离概率。O₂ 范围为模型施加的百分比，连续曲线的步长为 0.1%。

| 图 | 内容和环境 | p 模式 | 起始群体 | 时间 | seeds 汇总 |
| --- | --- | --- | --- | --- | --- |
| Figure8A | PCAWG/DepMap 器官倍性及三个模型比较点 | 模型使用原生拟合 p | 模型为 2N 起始；培养使用拟合分布 | 模型第1000天 | 模型合并 C01 500+C02 500；实证按器官/队列 |
| Figure8B | 作者论文 Figure5：WGD 比例与器官氧分压 | 不适用 | 不适用 | 不适用 | 19 条原始作者数据行 |
| Figure8C | in vivo，O₂=0–5% | 原生拟合 p | 2N、4N、1:1 混合 | 第1000天 | 合并1000 |
| Figure8D | in vitro，O₂=0–5% | 原生拟合 p | 2N、4N、1:1 混合 | 第1000天 | 合并1000，先在 seed 内平均传代重复 |
| supp_fig8_1 | A 的固定 p 对照 | 0.01 | 同 A | 第1000天 | 模型合并1000 |
| supp_fig8_2 | A 的固定 p 对照 | 0.3 | 同 A | 第1000天 | 模型合并1000 |
| supp_fig8_3 | in vivo，O₂=0–5%，C01/C02 分开 | 原生拟合 p | 七种 Q | 第1000天 | 每家族500 |
| supp_fig8_4 | in vivo，O₂=0–5%，两个固定 p panel | 0.01、0.3 | 七种 Q | 第1000天 | 每个 panel 合并1000 |
| supp_fig8_5 | in vitro，O₂=0–5%，行是固定 p，列是 C01/C02 | 0.01、0.3 | 七种 Q | 第1000天 | 每家族500，先平均传代重复 |
| supp_fig8_6 | in vitro，O₂=0–20%，2×5 时间矩阵 | 原生拟合 p | 2N、4N、1:1 | 100–1000天，每100天 | 合并1000 |
| supp_fig8_7 | in vivo，O₂=0–5%，2×5 时间矩阵 | 原生拟合 p | 2N、4N、1:1 | 同上 | 合并1000 |
| supp_fig8_8 | 将 supp_fig8_6 按家族拆分，上 C01、下 C02 | 原生拟合 p | 同 supp_fig8_6 | 同上 | 每个家族块500 |
| supp_fig8_9 | 将 supp_fig8_7 按家族拆分，上 C01、下 C02 | 原生拟合 p | 同 supp_fig8_7 | 同上 | 每个家族块500 |

### 1.2 当前完整流程

`run_Figure8.R` 直接读取每个 seed 的原始完整参数/环境配置及作者实证原始文件，独立生成有限时间预测，并在同一个当前 Figure8 目录写出合并和家族分组绘图表。正确培养配置、随机传代、终点与时间快照都在当前 worker 中处理，同一情景各时间点来自同一轨迹。`draw_Figure8.R` 不需要旧 run、fallback 或归档源码即可绘制新生成数据。

现有保留数据由独立只读兼容路径支持。早期提琴图、C/D 分别代表家族、主图 C/D 固定 p，以及单时间点补图6等设计仍已废弃。下文版本名只标识历史发布记录，不是计算阶段依赖。

## 2. 实证数据与原始资料来源

### 2.1 本地作者数据包

作者提供的数据位于只读输入目录：

```text
/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration5/data/share_bundle_figures_3_5
```

与当前图直接相关的文件为：

| 相对于数据包根目录的文件 | 用途 |
| --- | --- |
| `data/Figure3_PCAWG_included_samples.csv` | 作者 Figure3 映射/纳入的2338条肿瘤样本记录 |
| `data/Figure3_DepMap_included_lines.csv` | 作者 Figure3 映射/纳入的1433条细胞系记录；1185条具有有限倍性值 |
| `originals/organ_lut.csv` | 作者保存的19行 Figure5 输入，控制当前 B 的点和 SD |
| `supplementary/Supplementary_Table_3.csv` | 核对 Figure5 的类别、WGD 比例、氧均值和 SD |
| `scripts/oxygen_plots.Rmd` | 作者图形和普通最小二乘回归的原始分析依据 |
| `originals/nowicka.txt` | 较早的20行氧/WGD资料，仅用于版本差异审计 |
| `documentation/METHODS.md` | 数据包的映射、缺失值、来源及解释限制 |

不能把保存的均值、SD 或四分位数当成逐样本氧测量重新采样，也不从图像反推不存在的观测值。

### 2.2 论文来源及当前重现范围

Figure8 的实证比较对应 Beck 等的 [Energetics of whole genome doubling and genomic instability，Cancer Letters，2025，DOI:10.1016/j.canlet.2025.217878](https://pubmed.ncbi.nlm.nih.gov/40517938/)。其 Figure3 使用倍性中值和四分位数；Figure5 展示 WGD 肿瘤比例与器官氧分压，误差条为氧均值±SD。

当前 A 使用作者样本映射和有限倍性值重新汇总；当前 B 重现19行作者资料，不另建“经清洗的独立器官队列”。本地包保留的其他初步平均倍性–器官氧分析、primary-only/donor-level 汇总，不是当前 B 的输入。

## 3. 拟合参数、配置与模型文件来源

### 3.1 1000 个完整参数集

Figure8 使用 C01、C02 各500个 joint-fitting seeds，共1000个。每个 seed 使用自己最终的完整自然尺度参数，而不是每家族一个最优 seed，也不是只使用14个 soft-coupling 参数。

拟合根目录为：

```text
/Users/4482173/Documents/GitHub/soft_couping_org/oxygen/results/fit_joint_unified_global_invitro_500seed_all_xxlarge_r442_exact_20260828_145253
```

| Family | 家族目录 | 使用的 seeds |
| --- | --- | --- |
| C01 | `fit_joint_tsne_vi_seed366_C01_vt_seed144` | 1–500 |
| C02 | `fit_joint_tsne_vi_seed25_C02_vt_seed144` | 1–500 |

每个 seed 的 `joint_best_params_long.tsv` 中，`shared_invivo` scope 提供26个 in vivo 参数，`invitro_effective` scope 提供33个 in vitro 参数。`complete_input_signature`、输入清单及逐参数 passthrough 表记录其身份与数值。

不按拟合终点去重；两个 seeds 即使参数或预测相同，也各占一个 ensemble 位置。不交换不同 seeds 的参数，不按是否收敛、存活方向或是否符合预期趋势再次筛选。若参数在固定 O₂ 场景不活跃，仍保留其原始值，而不是用共同值替换。

### 3.2 in vitro 的权威配置

修复后的培养配置严格读取每个 seed 的 `fit_result.rds$ctx$invitro_cfg`，没有回退到 in vivo `fit_config.rds`。完整33参数与其最终表核对，保留初始高斯分布参数。

现行1000个培养配置均记录：`DT=0.05` 天、`Crowding=TRUE`、`crowding="logistic"`、`K=10^12`。这是实际配置核对结果，不是因为看到同样的数值就省略了每个 seed 的配置传递。

### 3.3 原生 p 与有效概率

原生拟合 p 的图保留每个 seed 自己的 `p_misseg`。固定 p 补图只通过 `figure6_force_p_misseg()` 覆盖该自然尺度参数为0.01或0.3，保留该 seed 的 `p_mis_base`、`k_o_mis` 及其他参数。

模型实际使用的有效概率为：

\[
p_{\mathrm{mis,eff}}(n,O_2)
=\operatorname{clip}_{[0,1]}\left[
p_{\mathrm{mis,base}}
+p_{\mathrm{misseg}}\frac{\mu(n,O_2)}{\mu(n,O_2)+k_{o,\mathrm{mis}}}
\right].
\]

因此固定原生参数不等于固定所有染色体状态的有效概率；不同 O₂、环境、家族和 seed 的状态速率仍由模型计算。Figure8 没有施加 Figure7 的 β 阻断。

### 3.4 实际模型代码

冻结输入指向：

```text
/Users/4482173/Documents/GitHub/soft_couping_org/oxygen/code/O2_supply_demand_MAP
```

核心包括固定氧入口、`model/model_O2_supply_demand_MAP.R/.cpp` 和相关 util。`f78_engine()` 按字节复制所需闭包到 run 内 `.model_snapshot`，核对 SHA256 后加载。不是因为当前仓库也存在同名 oxygen 目录，就自动使用其文件。

新的准备步骤可以由 `FIGURE_MODEL_CODE_ROOT`、`FIGURE_JOINT_RESULT_ROOT` 指定来源；已经发布结果的实际来源以冻结清单、算子和 fingerprint 为准。本次文档不切换模型或重新编译。

## 4. 单次群体倍性与起始条件

### 4.1 群体平均倍性的单位

当前染色体网格为22–154条，共133个状态，`N_UNIT=22`。设 \(X_{sri}(t)\) 为 seed s、培养重复 r、状态 i 的存活细胞数，\(n_i\) 为染色体数：

\[
P_{sr}(t)=\frac{\sum_i(n_i/22)X_{sri}(t)}{\sum_iX_{sri}(t)}.
\]

单位为 N。这是当时存活群体的细胞数量加权平均倍性，不是初始细胞的平均、不包括死细胞，也不是细胞倍性的中位数。in vivo 没有随机传代重复，可省略 r。

实现保存归一化状态 \(x_i=X_i/\sum X\) 并单独追踪群体大小，因此同一数值可计算为 \(\sum_i(n_i/22)x_i\)。这避免长期增长的浮点溢出，不改变加权定义。

### 4.2 in vivo 的纯状态起始

in vivo 纯2N、纯4N分别在44、88条状态初始化为单位质量；不是从拟合的培养分布中抽取。每个 seed 在同一模型、同一 O₂ 下分别处理这些起始。

### 4.3 in vitro 的拟合高斯起始

培养中的“2N”“4N”指相应起始 cohort，而不是强制所有细胞恰为44或88条。每个 seed 使用自己的 `init_mean_2N`、`init_sd_2N` 或4N对应参数：

\[
g_{s,k,i}=\frac{\phi(n_i;m_{s,k},\sigma_{s,k})}
{\sum_j\phi(n_j;m_{s,k},\sigma_{s,k})},\qquad k\in\{2N,4N\}.
\]

\(\phi\) 是在染色体数网格上求值的正态密度；参数单位与该网格一致。实际 helper 为 `ivt_gaussian_initial_state()`：检查有限均值和正标准差，在现有有限网格上用 `dnorm` 求值并归一化。不是把密度直接当绝对细胞数，也不是对每个格子做连续区间积分。

### 4.4 混合比例 Q

\[
Q=\frac{\text{初始4N-cohort细胞数}}{\text{初始2N-cohort细胞数}},\qquad
f_2=\frac{1}{1+Q},\quad f_4=\frac{Q}{1+Q}.
\]

纯2N为 Q=0，纯4N为 Q=∞，2N:4N=1:1–1:5 对应 Q=1–5。给定两个各自归一化的起始向量，混合起始为：

\[
x_s(0;Q)=f_2x_{s,2}(0)+f_4x_{s,4}(0).
\]

三种主图起始只使用 Q=0、1、∞；七种补图起始还包括 Q=2、3、4、5。混合与纯起始的总接种量相同，1:1 不意味着接种两倍总细胞。Q 是初始 cohort 数量比例，不是终点谱系比例；1:1混合也不是一个纯3N状态。

## 5. in vivo：固定氧的有限时间传播

### 5.1 连续增长算子

模型提供每个状态的分裂率 λ、直接死亡率 μ、普通分裂存活子细胞权重 B 和 WGD 权重 W。行是子状态、列是母状态：

\[
J=B\operatorname{diag}[\lambda(1-p_w)]
+W\operatorname{diag}(\lambda p_w),\qquad
G=J-\operatorname{diag}(\lambda),\qquad M=G-\operatorname{diag}(\mu).
\]

子细胞数量和可存活性已经包含在 B/W 中，不能再重复乘子细胞数。固定 O₂、参数及环境后，连续 in vivo 系统为：

\[
\dot X=MX,\qquad X(t)=\exp(tM)X(0).
\]

主图取 t=1000天，时间矩阵取指定有限 t；没有把所有起始强行设置为主特征向量，也不以已经收敛为前提。

### 5.2 谱移位只用于数值稳定

对 \(a=\max\operatorname{Re}\operatorname{eig}(M)\)，计算 \(\exp[t(M-aI)]\)，将共同因子 \(\exp(at)\) 保存在 log mass 中。共同因子在归一化倍性中抵消，但在组合不同起始 cohort 的数量时必须保留。

使用特征值选择稳定性移位，不等于使用主特征向量作为稳态预测。时间快照按100天矩阵指数连续推进，保留每段的状态归一化因子与 log mass。

### 5.3 混合群体必须按终点数量组合

连续系统对初始向量线性。令纯起始扩增因子为 \(m_2(t),m_4(t)\)，对应群体平均倍性为 \(P_2(t),P_4(t)\)。则：

\[
P_Q(t)=\frac{f_2m_2P_2+f_4m_4P_4}{f_2m_2+f_4m_4}
=w_2P_2+(1-w_2)P_4,
\]

\[
w_2=\operatorname{logistic}[\ell_2-\ell_4-\log Q],\qquad
\ell_k=\log m_k,\quad 0<Q<\infty.
\]

纯起始极限分别取 w₂=1或0。不能用初始权重直接平均两个终点倍性，除非两个终点数量因子恰好相同。

主图 C/七Q的 in vivo 来源使用这个精确线性组合；时间矩阵直接同时传播纯2N、1:1、纯4N三个初始向量。两种方法经过有限时间一致性检验。某些混合曲线贴近2N曲线，可以由2N来源的终点数量优势造成，不需要假定混合被错误初始化成2N。

## 6. in vitro：有限时间动力学及随机传代

### 6.1 培养动力学不是恒定线性矩阵加稀释

每个 seed 使用权威培养配置和完整参数。当前固定 O₂ 生产路径中，实际总活细胞数为 \(\mathcal N(t)\)，密度反馈为：

\[
c(t)=\max[0,1-\mathcal N(t)/K_{\mathrm{cap}}],\qquad K_{\mathrm{cap}}=10^{12}.
\]

在0.05天步长上，绝对计数 Euler 更新为：

\[
X^+=X+\Delta t\left[c(t)GX-\operatorname{diag}(\mu)X\right],
\qquad \Delta t=0.05\ \mathrm{day}.
\]

普通分裂/WGD产出和母细胞移除乘 c，直接死亡不乘 c。每步同步更新归一化状态与总数量，每次传代后也重新计算密度。当前实现不以之前错误的恒定线性培养传播器替代该过程。

固定 O₂ 和下述阈值随机传代是批准的预测设计；保留拟合动力学和配置，不代表逐条复制原始实验的可变氧及全部传代日程。

### 6.2 接种量及阈值来源

现行保存的 canonical rule 来源于六条 lineage 的114次 maintenance passage 的等权几何汇总：接种量在整数化前为569799.21840352，扩增比约10.84729935146，目标量在整数化前为6180782.69225091。

具体地，对 lineage ℓ 的各次有效记录 j，先分别对接种量 S_ℓj 与扩增比 T_ℓj/S_ℓj 取几何均值，得到 S_ℓ、R_ℓ；再对六个 lineage 的均值等权取几何均值：

\[
S_* = \exp\left[\frac{1}{6}\sum_{\ell=1}^{6}\log S_\ell\right],\qquad
R_* = \exp\left[\frac{1}{6}\sum_{\ell=1}^{6}\log R_\ell\right],\qquad
N_{\mathrm{seed}}=\operatorname{round}(S_*),\quad
N_{\mathrm{threshold}}=\lceil S_*R_*\rceil.
\]

因此每个 lineage 等权，不是把114次记录直接混在一起、让传代次数较多的 lineage 自动占更高权重。目标量使用向上取整，接种量使用四舍五入；约5天来自各 lineage 传代时长中位数的中位数。

实际接种569799个活细胞，阈值6180783个活细胞。它们是既有、数据推导并已批准的共同实验协议，不是1000个 seed 各自拟合出来的参数，也不是任意挑选的单次细胞数。

### 6.3 每次传代的具体顺序

1. 从上次接种后的实际群体状态连续推进；每天20个 Euler 子步。
2. 在每个实际整数天的生长更新后检查总活细胞数。首次达到阈值时执行传代；未达到就继续等待。
3. 用当时真实总量，而不是恰好等于阈值的假想总量，将期望状态数量做系统随机整数化，总数取最近整数。
4. 从整数化的各染色体状态计数中，通过顺序超几何抽样实现多元超几何无放回抽样，共保留569799个细胞。
5. 更新状态组成并将总量重置为接种量，在同一实际时间轴继续推进。

规则中约5天是代表性实验时长，不是每5天强制传代；不重置或倒退模拟天数。无阈值跨越就不人为补一次传代。传代弃掉的细胞与直接生物学死亡是不同事件。

### 6.4 重复、随机数与确定性例外

每个 seed×O₂×p模式×起始群体有200次名义传代重复。随机流由 family、seed、O₂、p、Q/起始条件等信息确定；不同重复有独立子流，可复现抽样。优化 seed 编号与传代随机重复不是同一概念。

实现先计算保守上界：

\[
r_{\mathrm{ub}}=\max_j\left\{\max\left[\sum_iG_{ij},0\right]-\mu_j\right\}.
\]

若 \(1000r_{\mathrm{ub}}<\log(6180783/569799)\)，则在 c∈[0,1] 下无法达到传代阈值，所有名义重复相同；只执行一个确定性 Euler 轨迹，并分别保存 `executed_repeats=1` 与 `n_stochastic_repeats=200`。这不是把无传代培养改算成 in vivo。

### 6.5 培养混合必须直接模拟

每个培养混合从该 seed 的混合拟合高斯向量开始，直接运行密度反馈和阈值抽样。混合组成影响增长、阈值日期和样本，因此不能通过两个纯群体终点倍性的简单平均或 in vivo 的线性叠加代替。

修复 run v11 重新生成了培养预测及其实际传代；没有把旧 v4/v9/v10 的培养终点或传代强行拼入。v12 延用相同认证引擎，以保存每日值的连续轨迹取得所有快照；不是对每个快照重新启动模拟。

## 7. 中值与 IQR：三个统计层级

### 7.1 seed 内先平均传代重复

对固定条件，先得到每个培养重复的群体平均倍性 \(P_{sr}(t)\)，再对同一 seed 做算术平均：

\[
Y_s(t)=\frac{1}{R_s}\sum_{r=1}^{R_s}P_{sr}(t).
\]

重复之间等权，不按各次终点细胞数量再加权；细胞数量权重已经用于每次重复内部的群体平均。无传代的确定性例外取 R_s=1。in vivo 的 Y_s 就是该 seed 的确定性群体平均倍性。

### 7.2 seeds 之间再取中值和四分位数

对当前条件的 seed 集合 S：

\[
\mathrm{center}=Q_{0.5}\{Y_s:s\in S\},\qquad
\mathrm{lower}=Q_{0.25}\{Y_s:s\in S\},\qquad
\mathrm{upper}=Q_{0.75}\{Y_s:s\in S\}.
\]

合并图 S 包含1000个身份 `(pair_label, seed_number)`；家族图 S 只包含该家族500个身份。不是合并全部200000次培养重复后取分位数，也不是将细胞、重复和优化 seeds 混作同一种样本。

R 使用 type-7 分位数：对排序值 y，\(h=1+(n-1)p\)，在 \(y_{\lfloor h\rfloor}\) 与下一个值之间线性插值。图上的 IQR 条或阴影画 Q25 至 Q75，而 IQR 宽度是 Q75−Q25。

### 7.3 合并两家族的正确顺序

先拼接 C01 的500个 Y_s 和 C02 的500个 Y_s，再计算分位数。两个家族各占500个位置，每个 seed 等权。不平均家族中值或四分位数，不根据 objective、终点群体大小或拟合终点重复次数重新加权。

### 7.4 区间不是什么

模型 IQR 是优化拟合结果集合所对应预测的描述性离散范围，不是90%/95%置信区间、后验可信区间或生物学重复间差异。均值的 Monte Carlo 标准误另存为 `stochastic_mcse`，不替代图上的 seed IQR。

相同或几乎相同的预测可以使 IQR 很薄或为零；图按真实宽度显示，不加抖动或人为放大。即使多次优化落到同一参数，也不会因此去掉这些 seed。

## 8. Figure8A：实证与模型的倍性中值/IQR

### 8.1 实证样本汇总

按作者 Figure3 的 `figure3_label` 分组，PCAWG、DepMap 分别计算有限 `ploidy` 的中值和 Q25/Q75。两个队列共有18个映射器官类别，图中有36个实证组和3个模型组，共39行。

PCAWG输入有2338条唯一 `sample_id`；DepMap输入有1433条唯一 `model_id`，其中248条无有限倍性，1185条实际参与倍性汇总。不把缺失值当0，也不把没有倍性的纳入记录算入该组有效 n。

这条重现路径没有额外做 primary-only、纯度、WGD不确定性或新增QC筛选，也不先按 donor 求均值。作者 caption 的 primary-tumor 说法与本地实际纳入的部分 specimen types 不完全一致，本文明确保留该限制。DepMap 24Q2 是原资料归属，未据此声称独立验证了整个归档版本。

### 8.2 三个模型比较点

模型均为第1000天、2N-cohort起始，合并1000 seeds：in vivo O₂=0.5%，in vitro O₂=0.5%，in vitro O₂=20%。培养从拟合2N高斯分布开始，使用真实阈值传代。

实证和模型均画 N 单位的中值与 Q25–Q75，但其分布单位不同：实证是肿瘤/细胞系记录，模型是 seed 内重复平均后的预测。相同单位与画法不意味着两者是同一种不确定性或配对样本。

## 9. Figure8B：作者 Figure5 的重现

### 9.1 点、坐标与误差条

横坐标为作者给出的 WGD 肿瘤比例，不是平均倍性；纵坐标为报告的器官氧分压，图轴按原图语境标为 mmHg。19行 `organ_lut.csv` 与 Supplementary Table3 的类别、WGD、氧均值和 SD 一致；其中一行没有 SD，因此不绘制该行的误差条。

每条有 SD 的纵向误差条为：

\[
[\bar O-\mathrm{SD},\ \bar O+\mathrm{SD}].
\]

这不是IQR，也不是由测量次数计算的置信区间。若下界为负，保留原始均值±SD，不通过截断制造另一套统计量；不据此声称真实氧测量为负。

### 9.2 回归与置信带

使用未加权普通最小二乘：

\[
O_k=\alpha+\beta_{\mathrm{OLS}}w_k+\varepsilon_k,
\]

w 为 WGD 比例。`geom_smooth(method="lm", se=TRUE, level=.95)` 画拟合均值的默认95%置信带，不是新观测的预测区间，也不是器官测量 CI。SD条不作为拟合权重；这里的回归系数 β_OLS 与 Figure7 阻断 β 无关。

### 9.3 行依赖、映射与统计版本差异

为重现原图，保留重复 kidney/marrow 行、共享广义肿瘤器官分组的不同 brain/lung 氧资料，以及原作者表观 suprarenal→thyroid 映射问题。不能把19行称为19个独立器官队列；没有静默进行器官去重、解剖纠错或亚型/primary-only筛选。

本地19行重新计算 Pearson R=0.616415699813，双侧 p=0.00494393936871。论文 caption 报告 R=0.62、p=0.0013；较早20行资料给出约 R=0.6670465、p=0.0013154。它们显示数据版本/统计结果不一致，不能仅凭此确定原因，也不把 caption 的 p 值硬编码为重算值。

本地氧表没有独立的逐测量数据、测量 N、完整 SD 推导或明确单位元数据；mmHg 遵从原图的氧分压语境，不是模型百分比与器官氧之间的新校准。因此该回归是原图的描述性重现，不产生新的独立跨器官因果推论。

## 10. Figure8C/D：两环境的主曲线

C、D 均为第1000天，O₂=0–5%，51个坐标，三个起始，每坐标1000 seeds，共各153个中值/IQR点。所有 seeds 保留自己的原生 p。C 为无培养传代的 in vivo；D 为拟合配置下有阈值传代的 in vitro。

两图均为实线：2N蓝色 `#0072B2`、4N粉色 `#CC79A7`、1:1混合橙色 `#E69F00`，主图不标 Q。曲线颜色和共同坐标不改变底层统计量。

A 的低氧 in vivo 2N 点与 C 的相应条件一致；A 的低氧培养点与 D 一致；A 的20%培养点位于扩展补图中。D 是 supp_fig8_6 第1000天 panel 的精确0–5%数值子集。

## 11. supp_fig8_1–5：固定 p 与多混合比例

supp_fig8_1/2重复 A 的实证比较和三个模型条件，仅将原生 p 分别覆盖为0.01或0.3。实证记录不因模型 p 的改变而改变；模型仍各点1000 seeds，培养模型仍先平均传代重复。

supp_fig8_3为原生 p 的 in vivo 七Q曲线，C01/C02分开，每个坐标500 seeds。supp_fig8_4则为固定 p=0.01、0.3 的两个匹配 in vivo panel，每个坐标合并1000 seeds，不能误写成同样按家族拆分。

七Q为0、1、2、3、4、5、∞。固定0.01的数据来自 v7 混合比例结果；固定0.3来自 v8 对既有完整参数纯终点的扩展，均按第5节的终点数量权重组合。

supp_fig8_5为培养固定 p 的2×2矩阵：列C01/C02，行0.01/0.3，七Q，每个曲线点500 seeds。培养的所有混合均直接模拟。它有1428个中值/IQR点，即两家族×两p×七Q×51氧坐标。20%的固定 p 纯2N额外比较点用于补图1/2的 A 式比较，不属于补图5的曲线范围。

## 12. supp_fig8_6/7：有限时间矩阵

两个图都是按行优先排列的2×5矩阵：上排100、200、300、400、500天，下排600、700、800、900、1000天。每个元素画三个起始的氧–倍性中值/IQR。

supp_fig8_6为 in vitro 0–20%，201氧坐标，共6030000条 seed×时间×条件记录和6030个汇总点。supp_fig8_7为 in vivo 0–5%，51氧坐标，共1530000条记录和1530个汇总点。

in vivo 用100天谱移位矩阵指数连续推进三个初始向量。培养每个条件只运行一条完整1000天过程的每次重复，保存所需时间点；快照取在该整数天的任何传代完成之后，不在100、200等时刻重新接种或重置。

第1000天 panel 分别还原此前的603个培养汇总点和153个主图 C 汇总点。原单panel supp_fig8_6保留在当时时间矩阵 audit 的 `before/`，不是当前补图6。

## 13. supp_fig8_8/9：家族拆分矩阵

补图8直接把补图6的已有 seed 快照按 family 分组；补图9对应补图7。每图上方为C01、下方为C02，各家族块均是完整2×5时间矩阵，保持共同坐标、三种起始和颜色。

每个家族×时间×起始×氧坐标重新对该家族500个 Y_s 取 type-7分位数。因此补图8有12060个汇总点，补图9有3060个；没有新的模拟、p覆盖、传代抽样或 seed 去重。

从同一原始快照重新合并两家族可以还原补图6/7，但不能通过平均补图8/9中的两家族中值/IQR来还原它们。

## 14. 实际数值与独立核对算例

### 14.1 原生 p、第1000天、合并1000 seeds

下表由保存的 seed 终点独立重算；均为 N 单位。

| 环境 | O₂（%） | 起始/Q | Q25 | 中值 | Q75 |
| --- | ---: | --- | ---: | ---: | ---: |
| in vivo | 0.5 | 2N / 0 | 1.02247325325 | 1.94841629820 | 2.77625918248 |
| in vivo | 0.5 | 1:1 / 1 | 1.08828221615 | 2.09871161189 | 2.77625918248 |
| in vivo | 0.5 | 4N / ∞ | 2.38555689654 | 2.77376137893 | 2.77625918249 |
| in vitro | 0.5 | 2N / 0 | 1.05054978945 | 1.05113082191 | 1.05177089750 |
| in vitro | 0.5 | 1:1 / 1 | 1.06759770428 | 1.06844011929 | 1.06928354859 |
| in vitro | 0.5 | 4N / ∞ | 3.58926924619 | 3.58950572552 | 3.58960961060 |
| in vitro | 20 | 2N / 0 | 5.99349544150 | 5.99360786483 | 5.99370919149 |
| in vitro | 20 | 1:1 / 1 | 5.99349607836 | 5.99360367889 | 5.99370376412 |
| in vitro | 20 | 4N / ∞ | 5.99350403120 | 5.99360572860 | 5.99372045777 |

前三个 A 模型比较条件分别是上表的 in vivo 0.5%/2N、in vitro 0.5%/2N、in vitro 20%/2N；与 `Figure8A_display_summary.tsv` 一致。

### 14.2 本次文档核对的覆盖范围

独立重算了 C、D、补图6/7第1000天、补图8/9第1000天家族分组、补图5的全部分位数，共4002个汇总坐标；保存值最大差异约6.1×10⁻¹⁴，属于文本序列化/浮点精度范围，使用10⁻¹⁰检验容差。

另从实证原始样本重算36组 A 的中值/IQR，最大差异约4.9×10⁻¹⁵；重算19行 B 的 Pearson统计；对全部原生 p in vivo Q=1–5混合终点核对数量加权公式。

这些是从已有数值做的独立后处理核对，不是重新积分模型或重跑传代。seed 内均值定义和原先逐重复数值检查依据实际 worker 与原验证记录；本次没有重新检查每个生产任务的完整逐重复数组。

## 15. 只读复算代码

### 15.1 模型 seed 汇总与时间/家族子集

以下 R 代码只读取现存终点和汇总表，无模拟、Rcpp编译、随机数或文件写入。主表约1319000行；需有相应读取内存，但不需要完整生产任务 RDS。

```r
library(data.table)
workspace <- "/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration5"
base <- file.path(workspace, "data/Figures/Figure8")
vivo_run <- file.path(base, "figure8_fitted_family_20261006_v8")
culture_run <- file.path(base, "figure8_invitro_repair_20261007_v11")
time_run <- file.path(base, "figure8_time_snapshots_20261007_v12")
family_run <- file.path(base, "figure8_family_time_snapshots_20261008_v13")
vivo <- fread(file.path(vivo_run, "fitted_mixture_seed_endpoints.tsv"),
  select = c("pair_label", "seed_number", "O2_pct", "Q", "mean_ploidy"))
culture <- fread(file.path(culture_run, "culture_seed_endpoints.tsv"),
  select = c("pair_label", "seed_number", "O2_pct", "Q", "p_setting", "mean_ploidy"))
summarize <- function(z, keys) {
  z[, .(q25 = unname(quantile(mean_ploidy, .25, type = 7)),
        median = median(mean_ploidy),
        q75 = unname(quantile(mean_ploidy, .75, type = 7)), n = .N), by = keys]
}
check <- function(calculated, saved, keys, label) {
  z <- merge(calculated, saved, by = keys, suffixes = c(".calc", ".saved"))
  error <- max(vapply(c("q25", "median", "q75"), function(k) {
    max(abs(z[[paste0(k, ".calc")]] - z[[paste0(k, ".saved")]]))
  }, 0.0))
  stopifnot(nrow(z) == nrow(calculated), nrow(z) == nrow(saved), error < 1e-10)
  cat("PASS", label, nrow(z), "points; max error", error, "\n")
}
keys <- c("O2_pct", "Q")
C <- summarize(vivo[Q %in% c(0, 1, Inf)], keys)
D <- summarize(culture[p_setting == "fitted" & O2_pct <= 5], keys)
extended <- summarize(culture[p_setting == "fitted"], keys)
check(C, fread(file.path(culture_run, "Figure8C_curve_summary.tsv")), keys, "Figure8C")
check(D, fread(file.path(culture_run, "Figure8D_curve_summary.tsv")), keys, "Figure8D")
check(extended, fread(file.path(time_run, "supp_fig8_6_curve_summary.tsv"))[day == 1000],
      keys, "supp_fig8_6/day1000")
check(C, fread(file.path(time_run, "supp_fig8_7_curve_summary.tsv"))[day == 1000],
      keys, "supp_fig8_7/day1000")
family_keys <- c("pair_label", keys)
check(summarize(culture[p_setting == "fitted"], family_keys),
      fread(file.path(family_run, "supp_fig8_8_curve_summary.tsv"))[day == 1000],
      family_keys, "supp_fig8_8/day1000")
check(summarize(vivo[Q %in% c(0, 1, Inf)], family_keys),
      fread(file.path(family_run, "supp_fig8_9_curve_summary.tsv"))[day == 1000],
      family_keys, "supp_fig8_9/day1000")
fixed_keys <- c(family_keys, "p_setting")
fixed_saved <- fread(file.path(culture_run, "supp_fig8_5_curve_summary.tsv"))
fixed_saved[, p_setting := as.character(p_setting)]
check(summarize(culture[p_setting != "fitted" & O2_pct <= 5], fixed_keys),
      fixed_saved, fixed_keys, "supp_fig8_5")
print(rbindlist(list(cbind(context = "in vivo", C[abs(O2_pct - .5) < 1e-10]),
                    cbind(context = "in vitro", extended[O2_pct %in% c(.5, 20)]))),
      digits = 12)
```

### 15.2 实证分位数与 B 的 Pearson统计

```r
library(data.table)
workspace <- "/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration5"
src <- file.path(workspace, "data/share_bundle_figures_3_5")
run <- file.path(workspace, "data/Figures/Figure8/figure8_invitro_repair_20261007_v11")
pc <- fread(file.path(src, "data/Figure3_PCAWG_included_samples.csv"))
dm <- fread(file.path(src, "data/Figure3_DepMap_included_lines.csv"))
stopifnot(nrow(pc) == 2338L, nrow(dm) == 1433L,
          sum(is.finite(dm$ploidy)) == 1185L)
raw <- rbindlist(list(
  pc[is.finite(ploidy), .(site = figure3_label, ploidy, cohort = "PCAWG")],
  dm[is.finite(ploidy), .(site = figure3_label, ploidy, cohort = "DepMap")]))
calculated <- raw[, .(median = median(ploidy),
  q25 = unname(quantile(ploidy, .25, type = 7)),
  q75 = unname(quantile(ploidy, .75, type = 7))), by = .(site, cohort)]
saved <- fread(file.path(run, "empirical_tissue_ploidy_summary.tsv"))
z <- merge(calculated, saved, by = c("site", "cohort"),
           suffixes = c(".calc", ".saved"))
stopifnot(nrow(z) == 36L)
for (k in c("median", "q25", "q75")) {
  stopifnot(max(abs(z[[paste0(k, ".calc")]] - z[[paste0(k, ".saved")]])) < 1e-12)
}
author <- fread(file.path(src, "originals/organ_lut.csv"))
stopifnot(nrow(author) == 19L, sum(is.na(author$SD)) == 1L)
test <- cor.test(author$wgd, author$PtO2.normal.tissue, method = "pearson")
print(c(n = nrow(author), R = unname(test$estimate), p = test$p.value), digits = 12)
fit <- lm(PtO2.normal.tissue ~ wgd, data = author)
print(coef(fit), digits = 12)
```

## 16. 结果链条、代码与解释边界

### 16.1 数值来源索引

下列 run 均在本工作目录的 `data/Figures/Figure8/` 下；文件名是内部历史标识，不改变当前图号。

| Run | 当前用途及主要文件 |
| --- | --- |
| `figure8_empirical_full1000_20261001_v4` | 最初完整参数冻结、实证表及保留的 in vivo 比较；其旧培养配置不是当前权威 |
| `figure8_mixture_ratios_20261006_v7` | 固定0.01、七Q的 in vivo `mixture_ratio_seed_endpoints.tsv` 及汇总 |
| `figure8_fitted_family_20261006_v8` | 原生 p 的 `fitted_pure_seed_endpoints.tsv`、`fitted_mixture_seed_endpoints.tsv`、家族汇总；固定0.3扩展表 `supp_fig8_4_seed_endpoints.tsv` |
| `figure8_invitro_repair_20261007_v11` | 权威修正培养结果 `culture_seed_endpoints.tsv`，主 C/D、A比较表和补图5；1319000个培养 seed/条件终点，6100任务 |
| `figure8_time_snapshots_20261007_v12` | `supp_fig8_6/7_seed_snapshots.tsv.gz` 及对应 `*_curve_summary.tsv`；连续轨迹快照 |
| `figure8_family_time_snapshots_20261008_v13` | `supp_fig8_8/9_curve_summary.tsv`；从 v12 快照按家族分组，无新模拟 |

v11 的1319000个终点包括603000个原生 p 三起始0–20%结果、714000个固定 p 七起始0–5%结果和2000个固定 p 的20%纯2N比较点。v12共5040任务，生产验证记录明确区分终点、时间快照和绘图点数量。

各 run 的 `source_manifest.tsv`、`code_manifest.tsv` 和 fingerprint 绑定输入与代码；`seed_input_manifest.tsv` 绑定 family/seed 身份；培养修复的 `authoritative_fit_inputs.tsv` 记录逐 seed 配置来源及 SHA，`full_native_parameter_passthrough.tsv` 记录逐参数数值核对。它们与实际汇总表共同构成证据链，不能只凭目录时间或某个 `completion.tsv` 就推断计算定义。

### 16.2 实现索引

| 内容 | 代码 |
| --- | --- |
| 主入口和修复优先级 | [draw_Figure8.R](/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration5/Code/Figures/draw_Figure8.R)、[draw_Figure78_invitro_repair.R](/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration5/Code/Figures/draw_Figure78_invitro_repair.R) |
| 完整参数和起始分布 | [figure78_full_parameters.R](/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration5/Code/Figures/util/analysis/figure78_full_parameters.R)：`f78v4_freeze`、`f78v4_initial` |
| 原生 p 的有限时间 in vivo 与混合 | [figure8_fitted_family_20261006.R](/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration5/Code/Figures/util/analysis/figure8_fitted_family_20261006.R)：`f8f_pure_endpoint`、`f8f_expand`、`f8f_summary` |
| 固定0.01多Q | [figure8_mixture_ratios_20261006.R](/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration5/Code/Figures/util/analysis/figure8_mixture_ratios_20261006.R)：`f8q_compute` |
| 培养混合初始与跨seed汇总 | [figure8_pooled_passage_20261007.R](/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration5/Code/Figures/util/analysis/figure8_pooled_passage_20261007.R)：`f8p_initial`、`f8p_summary`；其早期传播函数不替代当前修复 worker |
| 权威配置、当前培养 worker | [figure78_invitro_repair.R](/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration5/Code/Figures/util/analysis/figure78_invitro_repair.R)：`f78i_correct_prepared`、`f78i_propagate`、`f78i_worker8` |
| Euler与阈值随机抽样 | [figure78_invitro_propagator.cpp](/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration5/Code/Figures/util/analysis/figure78_invitro_propagator.cpp)：`f78i_euler_cpp` |
| 修复后的汇总与A比较点 | [figure78_invitro_repair_postprocessing.R](/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration5/Code/Figures/util/analysis/figure78_invitro_repair_postprocessing.R)：`f78i_aggregate` |
| 原始实证A | [figure8_empirical_comparison.R](/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration5/Code/Figures/util/analysis/figure8_empirical_comparison.R)：`f8e_build`；其中其他旧B探索表不是当前B |
| 当前作者Figure5数据 | [figure78_revision_20261005.R](/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration5/Code/Figures/util/analysis/figure78_revision_20261005.R)：`f78r_empirical5` |
| A/B绘图 | [figure78_full_drawing.R](/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration5/Code/Figures/util/analysis/figure78_full_drawing.R)：`f78v4_tissue_plot`、`f78v5_author_figure5`；[figure8_ab_display_20261007.R](/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration5/Code/Figures/util/analysis/figure8_ab_display_20261007.R) |
| 连续时间快照 | [figure8_time_snapshots_20261007.R](/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration5/Code/Figures/util/analysis/figure8_time_snapshots_20261007.R)：`f8t_vivo`、`f8t_vitro`、`f8t_aggregate` |
| 家族拆分与type-7分位数 | [figure8_family_time_snapshots_20261008.R](/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration5/Code/Figures/util/analysis/figure8_family_time_snapshots_20261008.R)：`f8tf_summary`、`f8tf_prepare` |
| 原生p覆盖与路径 | [figure6_robustness.R](/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration5/Code/Figures/util/analysis/figure6_robustness.R)：`figure6_force_p_misseg`、`f6r_paths` |
| 数据推导的共同传代协议 | [figure6_full_range_q10.R](/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration5/Code/Figures/util/analysis/figure6_full_range_q10.R)：`f6g_canonical_passage_rule` |

关键配置如下；旧配置中的 panel 名称仍须服从第1节的现行入口与版本选择。

| 配置 | 用途 |
| --- | --- |
| [figure8_empirical_comparison.json](/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration5/Code/config/figure8_empirical_comparison.json) | 初始完整参数任务设置，不代表旧培养实现仍有效 |
| [figure8_mixture_ratios_20261006.json](/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration5/Code/config/figure8_mixture_ratios_20261006.json) | 固定0.01、七Q的 in vivo 扩展 |
| [figure8_three_start_20261006.json](/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration5/Code/config/figure8_three_start_20261006.json) | 固定0.3纯终点的父级三起始来源；七Q扩展由后续 v8 实现 |
| [figure8_fitted_family_20261006.json](/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration5/Code/config/figure8_fitted_family_20261006.json) | 原生 p in vivo 来源、家族及混合比例 |
| [figure78_invitro_repair_20261007.json](/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration5/Code/config/figure78_invitro_repair_20261007.json) | 当前培养配置修复、O₂范围、传代重复与汇总 |
| [figure8_time_snapshots_20261007.json](/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration5/Code/config/figure8_time_snapshots_20261007.json) | 补图6/7的时间、氧范围及快照定义 |
| [figure8_family_time_snapshots_20261008.json](/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration5/Code/config/figure8_family_time_snapshots_20261008.json) | 补图8/9的家族拆分和type-7分位数 |

### 16.3 验证记录与本地可复算范围

既有生产记录位于以下 audit 目录，不由本次文档重写：

```text
/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration5/audit/figure78_invitro_repair_20261007
/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration5/audit/figure8_time_snapshots_20261007
/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration5/audit/figure8_family_time_snapshots_20261008
```

包括权威参数/配置、原实验 lineage 参考重放、步长核对、记录快照不改变终点/传代、全部seed身份/数量/参数、终点与时间矩阵一致性、绘图和发布校验。家族拆分的既有记录为23项数值、44项chart、52项PDF/PNG和2页视觉复核，共121项；它们是原有生产验证，不是本次重新跑模型的结果。

本地保留完整 seed 终点、压缩快照和汇总表，但没有完整收集 v11/v12 的生产任务 RDS。未完成 checkpoint 不能作为已完成任务使用。因此第15节直接使用可用完整表复算，不声称本次重验了全部逐重复数组；如需该层审计，应从原生产归档取得，而不是偷偷补跑。本次没有 SSH/HPC提交、重新拟合、抽样或模型模拟。

### 16.4 生物学与数值解释边界

- 修正培养配置不保证模型重现观察到的缺氧2N升倍性。原生 p、0.5%氧、1000天的培养2N与1:1预测仍接近下边界，不能为满足趋势而修改结果。
- 当前网格22–154对应1–7N。原修复诊断中，低氧2N/1:1的终点质量高度集中于22–24条；高氧部分结果接近更高倍性。没有新增扩展网格敏感性实验，不能宣称平台与网格无关。
- 保存的 `minimum_log_population` 是各重复终点中的最小值，不是整条轨迹的最小值；终点边界质量也不能代表所有时间的边界质量。
- 分数形式的期望细胞数量和比例不等于已经观测到该数量的单细胞。原有分母、人口量和边界诊断应与曲线一起理解。
- 优化seeds不是生物学重复或后验样本；传代重复不是额外拟合seeds。模型IQR、实证IQR、氧SD和OLS置信带不能互换。
- 模型施加的O₂百分比没有校准成器官mmHg；原图的映射和非独立行限制仍保留。任何模型–实证一致性或因果解释都需要额外证据。
