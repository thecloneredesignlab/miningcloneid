# 四个 seed 的 eFAST 后续试运行：完成、核验与解读

## 主要结论

1. **计算完整成功，但 ploidy ST 尚未整体收敛。** 两个 smoke、280/280 主轨迹、四个 seed 汇总和两个最终比较作业全部完成，退出码均为 0:0。成功运行与指数收敛是两个不同检查。
2. **seed25/138：N4097→8193、M=4，并未带来全面改善。** 五个重点参数的 R90 通过数分别从 3/5→1/5、5/5→3/5；均值接近不代表五次重复一致。
3. **seed101/174：固定 N4097，M=8 不是更稳定的整体替代。** R90 通过数分别从 1/5→0/5、1/5→1/5；seed101 的高氧和 seed174 的低氧仍明显不稳定。
4. **目前不建议直接扩展到全部 500 个 seed 的更高 N 或统一改用 M8。** 优先检查已有响应的频谱、相位与氧过渡区域，并以相同模型和参数分布的独立方差估计设计做交叉核对。若继续 eFAST，只做有明确判别目标的少量 seed 试运行。

## 1. 本轮范围、数据与核验

- 独立 in-vivo 模型；四个拟合最优解的邻域。每个参数区间为最优值±全局自然尺度参数跨度的 10%，裁剪到全局界限；均匀/对数均匀变换保持不变。
- 14 个有效参数同时变化；O2 从 0 到 5%，步长 0.025%，共 201 点；每个 seed/config 五个相位重复。两个输出为 dominant mean ploidy 和 leading-eigenvalue net live-cell growth rate。
- seed25、138：N8193/M4；seed101、174：N4097/M8。每组 140 条完整焦点参数轨迹；合计 345,840,600 次主算子求值，两个输出共享一次求解。
- 计算代码为 `9cea6b22feee20cae9dea2cc64467c8b2199bc05`，执行时的投递凭据提交为 `0bc4237a3d598e7efc404f9f25132887c63753d4`。SALib 1.5.2、Python 3.9.25、原 SIF 和模型源文件均保持一致。
- 已核对每个任务的 manifest、相位、N/M、参数顺序、父设计、样本 SHA256、缓存 SHA256、输出完成凭据、恢复证明及汇总来源。280 个缓存和全部样本文件在 HPC 就地核验；完整原始输出 gzip 除恢复源文件外未重新计算哈希，其他输出以原完成凭据与缓存交叉核对。这不等于重新运行模型或独立重算全部原始 FFT。
- 本地收集 64 个文件、13,589,169 字节，所有 SHA256 匹配。原始轨迹保留在 HPC；新派生表和图保存在本目录。
- 从 12 份五相位数组独立核对 **586,656 项**数值或分类，最大数值差 `1.0547e-15`；涵盖 2,016 行重复诊断、672 行配置差异、135,072 行氧曲线摘要、672 行过程/排名摘要。
- 12 份指数数组没有未定义单元。八份可用方差数组共 225,120 个单元均有限且为正；最小方差约 `1.24e-8`。本轮问题不能解释为把 NaN 置零或完全恒定输出。

原始完成清单为 `completed_results_manifest.tsv`，HPC 审计为各配置的 `slurm/completed_hpc_audit.json`，独立核验为 `analysis_validation.json`。

## 2. 如何判定稳定：R90、D90 与局部最大值

对每个参数、每个 O2 点：

- `R(O2) = max(ST across 5 phases) - min(ST across 5 phases)`；`R90` 是该曲线在指定 O2 窗口的第 90 百分位，阈值 **0.10**。
- `D(O2) = abs(mean_test(ST) - mean_reference(ST))`；`D90` 阈值 **0.05**。这些是指数的绝对差单位，不是参数扰动百分比或 p 值。
- 本轮 reference 均为 N4097/M4。25/138 的 D 是同 M 的 N 分辨率变化；101/174 的 D 是同 N 的 M 设置变化。后者不能称为 N 收敛检查。
- `联合 p90`：测试配置 R90 和 D90 同时通过。`严格最大值联合通过`：测试配置 Rmax≤0.10 且 Dmax≤0.05。M 稳健性的更强要求还包括 reference 的 R90 通过，另列于 `both_configurations_R90_and_D90_pass`。
- p90 可允许少量 O2 点出现很大差异。必须同时看 Rmax、Dmax 及位置；达到 p90 阈值不等于整条 O2 曲线均稳定。
- 五次重复的范围不是置信区间，两个 N/M 配置不是嵌套样本，相同相位种子也不是同一组参数向量。这里的阈值是固定的实用诊断，不能证明估计无偏或达到某个置信水平。

三个窗口分别为全范围 [0,5]、低氧 [0,1]、较高氧 [3,5]；中间过渡区域不能因为窗口端点看起来稳定而被忽略。

## 3. 五个重点参数：全 0–5% 网格

表中 R90 参考→新配置；D90 为两配置五相位均值差。`p90通过/局部仍失败`表示 p90 通过、但至少一个最大值未通过。

### 3.1 seed25、138：固定 M=4，提高 N

| seed | 参数 | R90：参考→新配置 | D90 | 测试配置联合判定 |
|---|---|---|---|---|
| seed25 | buffer_beta | 0.0413 → 0.0478 | 0.0155 | 严格通过 |
| seed25 | buffer_smax | 0.1342 → 0.1233 | 0.0217 | 未通过 |
| seed25 | p_mis_base | 0.1125 → 0.1547 | 0.0284 | 未通过 |
| seed25 | mu_hp | 0.0750 → 0.1096 | 0.0519 | 未通过 |
| seed25 | O2_crit | 0.0806 → 0.1404 | 0.0451 | 未通过 |
| seed138 | buffer_beta | 0.0891 → 0.0780 | 0.0271 | 严格通过 |
| seed138 | buffer_smax | 0.0525 → 0.3080 | 0.1177 | 未通过 |
| seed138 | p_mis_base | 0.0423 → 0.2870 | 0.1585 | 未通过 |
| seed138 | mu_hp | 0.0535 → 0.0604 | 0.0317 | p90通过/局部仍失败 |
| seed138 | O2_crit | 0.0585 → 0.0988 | 0.0378 | p90通过/局部仍失败 |

### 3.2 seed101、174：固定 N4097，改变 M

| seed | 参数 | R90：参考→新配置 | D90 | 测试配置联合判定 |
|---|---|---|---|---|
| seed101 | buffer_beta | 0.0773 → 0.1091 | 0.0494 | 未通过 |
| seed101 | buffer_smax | 0.4938 → 0.5057 | 0.2555 | 未通过 |
| seed101 | p_mis_base | 0.5449 → 0.6369 | 0.2663 | 未通过 |
| seed101 | mu_hp | 0.5553 → 0.5903 | 0.2761 | 未通过 |
| seed101 | O2_crit | 0.2031 → 0.2101 | 0.0291 | 未通过 |
| seed174 | buffer_beta | 0.1842 → 0.1594 | 0.0771 | 未通过 |
| seed174 | buffer_smax | 0.2273 → 0.2810 | 0.0732 | 未通过 |
| seed174 | p_mis_base | 0.1426 → 0.0899 | 0.0464 | p90通过/局部仍失败 |
| seed174 | mu_hp | 0.1132 → 0.2035 | 0.0237 | 未通过 |
| seed174 | O2_crit | 0.0790 → 0.4656 | 0.1068 | 未通过 |

## 4. seed25：均值多数接近，重复稳定性仍不足

- 五个重点参数的 R90 从 **3/5→1/5**；只有 `buffer_beta` 同时严格通过重复和跨 N 差异检查。
- `buffer_beta`：R90 `0.0413→0.0478`，虽然略增，但仍在阈值内；D90 `0.0155`、Dmax `0.0303`，整条曲线的 ST 在本次检查中稳定。
- `buffer_smax`：R90 `0.1342→0.1233`，略有改善，但仍未通过；D90 `0.0217`、Dmax `0.0391`。这正是**均值已接近、相位重复仍不一致**的例子。
- `p_mis_base`：R90 `0.1125→0.1547`，变差；D90 `0.0284`、Dmax `0.0453`。不能仅凭橙/青色均值曲线重合宣布收敛。
- `mu_hp`：R90 `0.0750→0.1096`，从 p90 通过变为未通过；D90 `0.05188`、Dmax `0.06109`，跨 N 差异也略超阈值。
- `O2_crit`：R90 `0.0806→0.1404`，重复差异增加；D90 `0.04512` 通过但 Dmax `0.05430` 于 1.70% O2 未通过。
- 低氧五参数联合 p90 为 **3/5**：主要问题是 `mu_hp`、`O2_crit`。高氧也为 **3/5**：主要问题换成 `p_mis_base`、`buffer_smax`。因此不能把未通过全部归因于单个中间氧尖峰。
- 全 14 参数：新配置 R90 为 **10/14**、D90 **13/14**、联合 p90 **10/14**、严格联合 **7/14**。本轮均值的整体分辨率变化有限，但不能给出全参数全范围 ST 已收敛的结论。

## 5. seed138：N4097 的表面稳定在 N8193 下没有保持

- 五参数 R90 从 **5/5→3/5**；只有 `buffer_beta` 严格联合通过。`mu_hp` 和 `O2_crit` 虽联合 p90 通过，仍有明显局部失败。
- `buffer_beta`：R90 `0.0891→0.0780`，D90 `0.02714`、Dmax `0.04274`；它的 **ST** 在本轮检查中稳定。
- `buffer_smax`：R90 `0.0525→0.3080`，D90 `0.1177`、Dmax `0.1845`，最大均值差在 **2.575% O2**。
- `p_mis_base`：R90 `0.0423→0.2870`，D90 `0.1585`、Dmax `0.2968`，最大均值差在 **2.25% O2**。两者的 N8193 曲线在约 1.5–2.8% O2 再次出现大范围相位差异。
- `mu_hp`：R90 `0.0604`，但 Rmax `0.3499`（1.625%）；Dmax `0.1232`（1.775%）。`O2_crit`：R90 `0.0988`，但 Rmax `0.4883`（1.975%）、Dmax `0.2344`（1.925%）。这里 p90 隐藏的局部问题很大。
- 低氧五参数联合 p90 **5/5**，但严格联合 **3/5**；高氧 [3,5] 五参数 **5/5 严格通过**。所以**局部高氧结论比全 O2 范围结论更可靠**，主要困难在中间过渡区域。
- 全 14 参数 R90 从 **11/14→3/14**；新配置联合 p90 **3/14**，严格联合 **1/14**。N2049→4097 的改善不是单调趋势，N8193 不能作为“已解决 seed138”的证据。
- 还需单独注意 **S1 不等于 ST**：`buffer_beta` 的 ploidy S1 R90 从 `0.0946→0.4608`、D90 `0.1293`。ST 稳定并不自动保证该参数的 S1 也稳定。

## 6. seed101：M8 没有解决高氧严重不稳定

- 五个重点参数 M4→M8 的 R90 **全部增大**，通过数 **1/5→0/5**；全 14 参数 M8 R90 通过数也为 **0/14**。
- `p_mis_base`、`buffer_smax`、`mu_hp` 的 M4/M8 D90 分别为 `0.2663`、`0.2555`、`0.2761`，Dmax 分别达 `0.4043`、`0.3363`、`0.3969`。它们的高氧 ST 水平对频率设置高度敏感。
- `buffer_beta`：R90 `0.0773→0.1091`；D90 `0.04943` 虽刚好通过，Dmax `0.05719` 未通过。不能把边界通过读成强稳健证据。
- `O2_crit`：全范围 D90 只有 `0.02912`，但 Rmax `0.7649`（4.95%），Dmax `0.1572`（4.75%）。高氧窗口内 R90 `0.7497`、D90 `0.08985`，均失败。宽全域统计淡化了高氧局部的严重不稳定。
- 低氧五参数联合 p90 **4/5**、严格联合 **1/5**；高氧五参数 R90 和 D90 **全部失败**。这一 seed 不宜用高氧 ST 精确数值或过程排名作为主要机制证据。
- ploidy S1 也有问题：`buffer_beta` R90 `0.6252→0.6477`、D90 `0.1486`；`buffer_smax` D90 `0.06569`；`O2_crit` 新 R90 `0.10075`。因此不能笼统声称这一 seed 只影响 ST、不影响 S1。

## 7. seed174：少数全域 p90 改善掩盖低氧失败

- 五参数 R90 通过数保持 **1/5**，但通过者从 M4 的 `O2_crit` 换成 M8 的 `p_mis_base`。只有 `p_mis_base` 达到测试配置联合 p90；**没有任何重点参数严格联合通过**。
- `p_mis_base` 的全域 R90 `0.1426→0.0899`、D90 `0.04635`，表面改善；但 Rmax `0.5425`（0.375% O2）、Dmax `0.1819`（0.475%）。它不能被标为整条曲线可靠。
- `buffer_beta` R90 `0.1842→0.1594`，略降但仍失败；`buffer_smax` 从 `0.2273→0.2810`，变差。
- `mu_hp` R90 `0.1132→0.2035`；`O2_crit` 从 `0.0790→0.4656`，Rmax `0.6723` 于 0.55% O2。M8 的低氧重复范围明显变宽。
- 低氧五参数 **R90、D90、联合检查均为 0/5**；高氧联合 p90 **3/5**、严格联合 **2/5**。
- 全 14 参数 R90 **9/14→4/14**；新配置联合 p90 **4/14**、严格联合 **0/14**。若要求 M4 和 M8 的重复检查均通过，再加均值差通过，五个重点参数为 **0/5**，全 14 参数也仅 **3/14**。

## 8. 三个氧窗口的重点参数汇总

每格为通过的参数数目/5；“联合”使用测试配置重复与均值差；严格联合同时要求两个最大值通过。

| seed | O2窗口 | R90 | D90 | 联合 p90 | 严格联合 |
|---|---|---|---|---|---|
| seed25 | 0–5% | 1/5 | 4/5 | 1/5 | 1/5 |
| seed25 | 0–1% | 3/5 | 5/5 | 3/5 | 3/5 |
| seed25 | 3–5% | 3/5 | 5/5 | 3/5 | 3/5 |
| seed138 | 0–5% | 3/5 | 3/5 | 3/5 | 1/5 |
| seed138 | 0–1% | 5/5 | 5/5 | 5/5 | 3/5 |
| seed138 | 3–5% | 5/5 | 5/5 | 5/5 | 5/5 |
| seed101 | 0–5% | 0/5 | 2/5 | 0/5 | 0/5 |
| seed101 | 0–1% | 4/5 | 5/5 | 4/5 | 1/5 |
| seed101 | 3–5% | 0/5 | 0/5 | 0/5 | 0/5 |
| seed174 | 0–5% | 1/5 | 2/5 | 1/5 | 0/5 |
| seed174 | 0–1% | 0/5 | 0/5 | 0/5 | 0/5 |
| seed174 | 3–5% | 3/5 | 5/5 | 3/5 | 2/5 |

## 9. Growth 与 S1：能提供什么对照

- N8193/M4 两 seed 的 growth S1/ST 均为 **28/28**：重复、跨配置差异、严格最大值联合均通过。
- N4097/M8 两 seed 的 growth S1/ST 均为 **28/28 p90 联合通过**；growth ST 的严格最大值联合为 **26/28**。两个例外都是 `O2_crit`：seed101 Rmax `0.1063`（0.30%），seed174 Rmax `0.1125`（0.675%），均轻度超 0.10；其 Dmax 均≤0.05。
- ploidy S1 总体比 ST 稳定：N 组联合 p90 **27/28**；M 组 **25/28**。但第 5、6 节列出的 `buffer_beta` 等例外必须单独保留。
- growth 的稳定支持两输出不能混为一谈：leading eigenvalue 与 leading-eigenvector-derived ploidy 有不同的数值/响应行为。这是观察到的差别，不能仅凭这一对照确定 ploidy 不稳定的具体原因。

## 10. 数值恢复是否能解释这些差异

- N8193/M4：seed138 有 1 条轨迹的 1 个点经过恢复；对 ploidy ST 的最大影响 `1.22e-15`。
- N4097/M8：seed101 的 5 条轨迹/6 点、seed174 的 8 条轨迹/9 点经过恢复；对 ploidy ST 的最大影响 `9.33e-13`，对 ploidy S1 最大影响 `1.14e-12`。
- 所有恢复证明和原始恢复源文件/矩阵哈希均通过；最大双精度残差约 `6.94e-18`。50/100 位结果在所记录精度内一致，记录的零差异不是数学上误差严格为零的证明。
- 最大单点 ploidy 修正约 `3.48e-9`；这些影响远小于观察到的 `0.05–0.7` 指数差异，不能解释当前 ST 不稳定。
- 因此已知恢复不是主要原因。仍不能从上述检查直接证明是频谱混叠、模型分岔或某一种生物学相互作用；这些是下一步需要鉴别的机制/数值解释。

## 11. 图怎么读

### 11.1 原始 `figures/ploidy_ST_comparison.pdf`

每组一张，两行分别是该组两个 seed，五列是 `buffer_beta`、`buffer_smax`、`p_mis_base`、`mu_hp`、`O2_crit`。

- 横轴 O2%，纵轴 ST；ST 高表示参数本身及其涉及的交互作用对该输出方差贡献大。ST 不提供增加参数后 ploidy 上升还是下降的方向，也不能把多个 ST 直接相加解释为互斥方差份额。
- 灰线是 N2049/M4 历史背景；橙线为 N4097/M4；青线为本组新配置。线是五相位均值，半透明色带是 min–max **完整重复范围**，不是置信区间。
- 浅蓝背景为低氧 [0,1]，浅橙背景为 [3,5]。背景只标氧窗口，与曲线色带含义不同。
- 同时看均值线与色带：seed25 的均值大多接近但色带仍较宽；seed138 的中间氧尖峰再现；seed101 高氧色带宽、均值随 M 改变；seed174 低氧青色色带变宽。
- `growth_ST_comparison.pdf` 用同样方法阅读，输出改为增长率的 ST。

### 11.2 新增 `figures/ploidy_ST_diagnostic_overview.pdf`

每组各一张，A 为参考配置 R90，B 为新配置 R90，C 为两配置均值差 D90；两行 seed，五列重点参数。**格内数字不是 ST 本身，是误差/离散诊断。**

- 灰：p90 和最大值都通过。
- 蓝色带 `*`：p90 通过，但局部最大值未通过。
- 橙：p90 未通过。
- 深灰 NA：不完整或未定义；本轮没有这一类。

例如 seed138 的 `mu_hp`、`O2_crit` 在 B/C 中为蓝色星号，表示“整体大部分网格通过，局部仍不稳定”；seed174 的 `p_mis_base` 也是同类问题。

六张 PDF 均渲染检查；原 HPC 四张对比 PDF 及其 PNG 保持原始哈希不变。

## 12. 对 Noemi 机制假设的影响

这里沿用前轮的描述性判据：低氧 death 过程的平均参数指数是否最大，missegregation 和 buffering 的过程均值是否在高氧增加。过程均值不是严格的分组 Sobol/ST 指数；death 为 `mu_hp/gamma_mu`，`O2_crit/n_O` 单列 shared oxygen stress，不能直接全部算作 death。

- 四个 seed/配置组合的 ploidy S1 和 ST 均为 **0/4** 满足三个条件；每个指数下的 20 个相位结果也都是 **0/20** 同时满足。这个比例只描述选定诊断 seed，不能推广为 500 个 seed 或总体概率。
- 没有一个新配置在这一过程均值判据下表现为低氧 death 最大。高氧 buffering 增强在 seed25、138 的 S1/ST 均保留；seed101 只在 ST 而非 S1 中保留，且其 ST 很不稳定；seed174 的 buffering 高减低为负。
- ploidy ST 的高减低符号在 N 组 27/28 参数-seed 对、M 组 24/28 对保持一致。但幅度、相位稳定性和过程排名并非都稳定，符号一致不能替代数值收敛。
- seed25 的 S1 missegregation 高减低只从 `-0.00142` 变到 `+0.00067`，属于接近零的小幅翻转，不能拿它作为强“高氧增强”证据。seed138 的 ST 低氧最高参数从 `O2_crit` 变为 `buffer_beta`，也说明不能只看一个配置的排名。

**可写出的结论：部分拟合邻域支持高氧 buffering 的重要性，但没有得到跨这四个诊断邻域一致的“低氧 death 主导→高氧 missegregation/buffering 主导”转移；未收敛的 ploidy ST 更不支持用精确幅度或排名强化这一机制叙事。** 这不等于证明死亡过程在原模型或完整拟合集合中不重要；Figure 4B 的相关方向信息应继续保留。

## 13. 下一步建议：先诊断估计器，不做 500-seed 扩容

### 建议顺序

1. **先用已有响应做针对性诊断。** 选 seed25 的 `p_mis_base/mu_hp/O2_crit`，seed138 的 1.5–2.8% O2，seed101 的 3–5%，seed174 的 0–1%；比较每个相位的响应分布、频谱功率分布及 ST 对分界频率的依赖。保留全部 14 参数的原设计，不能把删除参数后的运行当作相同问题。
2. **优先安排独立估计方法的交叉核对。** 在同一固定 O2 算子、同一最优解邻域和同一分布下，可采用 Sobol/Jansen 等独立采样设计，检查 ST 的量级及 O2 依赖是否一致；同样需要重复与分辨率诊断。这是建议，本轮没有新增或改动任何相关任务。
3. **若继续 eFAST，维持少量 seed，并固定 M 后再提高 N。** seed25/138 可把 N16385/M4 作为进一步判别，但当前没有证据保证下一次翻倍会解决问题。seed101/174 如要检验 M8 的 N 收敛，应比较 N4097/M8→N8193/M8；不能拿 N8193/M4 去替代这一比较。下一轮仍建议五次重复，不按是否通过删 seed。
4. **先使用已经稳定的结论。** seed25 的 buffer_beta ST、seed138 较高氧窗口的重点 ST、绝大部分 growth 指数可按各自诊断范围表述；全域不稳定参数作为诊断与限制报告。S1、ST、growth 分开写。

### 为什么 M8 不一定更稳

SALib 中 M 决定 Fourier 谐波求和数量，也影响设计频率。N4097/M4 的主频为 512、互补频率预算 64；N8193/M4 为 1024/128；N4097/M8 为 256/16。三者均满足 N>4M²，M8 的 16 仍足以给其余 13 参数分配不同频率。满足约束只说明设计合法，不保证非线性模型中的估计收敛。采样与分析都必须使用同一个 M。

公式与实现依据：[SALib 1.5.2 sampler](https://raw.githubusercontent.com/SALib/SALib/v1.5.2/src/SALib/sample/fast_sampler.py)、[SALib 1.5.2 analyzer](https://raw.githubusercontent.com/SALib/SALib/v1.5.2/src/SALib/analyze/fast.py)。本次仅把频率分配变化作为可能需要检查的因素，没有据此断言已证明混叠。

### 时间代价

- N8193/M4：seed25 单任务中位 **5.10 h**、最长 **5.27 h**；seed138 中位 **4.89 h**、最长 **5.52 h**。
- N4097/M8：seed101 中位 **2.54 h**、最长 **3.75 h**；seed174 中位 **2.45 h**、最长 **2.53 h**。中位时间与其 N4097/M4 基线接近。
- 按样本数线性外推（不是实测）：N16385/M4 的单任务中位约 **9.8–10.2 h**、历史最长外推 **10.5–11.0 h**，接近 12 h 限制，不能保证不会超时。N8193/M8 中位约 **4.9–5.1 h**、历史最长外推约 **5.1–7.5 h**。

**本轮没有投递新任务。** 上述仅是下一阶段选项，不是已执行计划。

## 14. 可供稿件使用的英文简述

Five-phase eFAST follow-up runs did not establish uniform convergence of ploidy total-effect indices. Increasing N from 4097 to 8193 at M=4 did not consistently reduce phase dependence in seeds 25 and 138, with substantial residual instability in the intermediate-oxygen response of seed 138. At fixed N=4097, increasing M from 4 to 8 did not provide a generally more stable alternative for seeds 101 and 174. Instability remained pronounced at higher oxygen in seed 101 and lower oxygen in seed 174. Growth-rate indices were substantially more stable. These selected-neighborhood diagnostics support cautious, output- and oxygen-specific interpretation rather than a universal mechanistic switch or precise ranking of unresolved ploidy total effects.

## 15. 文件与重现

- 原始对比图：各配置 `figures/ploidy_ST_comparison.pdf`、`growth_ST_comparison.pdf`。
- 新诊断图：各配置 `figures/ploidy_ST_diagnostic_overview.pdf`。
- 根目录 `ploidy_ST_key_parameter_details.tsv`：四 seed×五参数×三窗口，共 60 行；`ploidy_ST_parameter_details.tsv`：168 行。
- 各配置的全部输出诊断表、机制/相位表、运行时间表、数值恢复表和方差表均可追溯到收集清单与 NPZ。
- 报告脚本 `report_resolution_nm.py`；收集器 `collect_resolution_nm_results.py`；防误判测试 `test_resolution_nm_report.py`。独立统计助手为历史 `report_resolution_pilot.py`；冻结计算代码没有修改。

```bash
# HPC 就地审计与生成紧凑结果清单；不会提交任务
python3 oxygen/code/O2_supply_demand_MAP/analysis/eFAST/collect_resolution_nm_results.py audit --repo-root "$PWD"
python3 oxygen/code/O2_supply_demand_MAP/analysis/eFAST/collect_resolution_nm_results.py manifest --root oxygen/results/eFAST/neighborhood10pct/pilot/N8193_M4_and_N4097_M8

# 完成清单内文件同步到本地后核验；报告需要 NumPy>=2、Matplotlib
python3 oxygen/code/O2_supply_demand_MAP/analysis/eFAST/collect_resolution_nm_results.py verify-local --root oxygen/results/eFAST/neighborhood10pct/pilot/N8193_M4_and_N4097_M8
python3 oxygen/code/O2_supply_demand_MAP/analysis/eFAST/report_resolution_nm.py --root oxygen/results/eFAST/neighborhood10pct/pilot/N8193_M4_and_N4097_M8
python3 -m unittest discover -s oxygen/code/O2_supply_demand_MAP/analysis/eFAST -p test_resolution_nm_report.py -v
```

核验范围不包括重新运行模型或对全部原始输出做独立 FFT；原始输出留在 HPC，必要时可用于后续频谱定位。
