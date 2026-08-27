# `selected_day` 与 `duration_days` 的时间尺度不一致问题

核查日期：2026-07-30

## 结论摘要

当前体外模型中存在一个具有实质科学影响的时间对齐问题：

- 实验 `passage_duration` 是一次真实传代从接种到收获的记录间隔；
- 模型 `duration_days` 是根据实验记录建立的单段最大模拟时间窗；
- 模型 `selected_day` 是该时间窗内预测活细胞数最接近目标终末细胞数的模拟日。

三者定义不同。但是，模型在计算预测增长率和向下一次传代传递染色体组成时，实际使用的是 `selected_day`，而不是实验 `passage_duration` 或模型 `duration_days`。因此，从链式状态传播的角度看，`selected_day` 成了模型的“有效传代持续时间”。

当 `selected_day` 早于真实传代日时，模型仍执行相同次数的传代，却在更短的累计时间内完成这些传代。最严重的例子是 seed340 的 2N O1/O2 谱系：实验累计时间为 122 天，而模型累计 `selected_day` 只有 40 天，缩短 82 天，即压缩了 **67.2%** 的实验周期。这个结果不能解释为轻微的终点取值差异，而是模型时间尺度发生了实质改变。

此外，当前程序在进入下一次 passage 时，无条件把上一段 `selected_day` 的染色体组成缩放到下一段记录的 `initial_cells`。当预测活细胞数因负增长而低于下一次接种所需数量时，这一步不是稀释，而是人为增加活细胞总数。模型因此可以在没有足够细胞完成真实传代的情况下继续整个 passage 序列。该问题与时间压缩相互作用，是 seed340 仍能传播终点 ploidy 状态的重要结构性原因。

## 1. 实验如何进行

体外实验使用近二倍体（2N）和近四倍体（4N）的 SUM159 谱系。对每个起始倍体：

- control 分支维持在 20.5% O₂；
- O1 和 O2 是两条平行的低氧谱系，依次经历 2%、1%、0.5%、0.3%、0.2%、0.1% 和 0% O₂，随后主要维持在目标 0% O₂；
- 每次传代的实验操作单元为：按记录的初始细胞数接种，在指定 O₂ 条件下培养，经过记录的传代间隔后收获并计数，再取其中一部分细胞接种下一次传代；
- 每次传代可关联初始和终末校正细胞数、增长率，以及部分传代的单细胞染色体数和 G0/G1 DNA 含量分布。

O1 和 O2 在实验中是分别传代的平行谱系。现有记录不足以进一步确定它们是否应被视为严格意义上的生物学重复，因此本文件只称其为“平行谱系”。

实验设计的稿件依据见：

- `revised/iteration1/manuscript/ltee_hypoxia_model_round3_integrated.tex` 第 364–371 行。

## 2. 实验中的时间和细胞尺度

数据源：

`/Users/4482173/Documents/GitHub/soft_coupling/oxygen/ploidyOxygen/data/fit_objects/fit_data.Rds`

当前 `fit_data.Rds` 有 131 个条目，其中 114 个正式 A-passage 条目具有完整的 `passage_duration`、`initial_cells` 和 `final_cells`。其余 17 个 seed/landmark 条目没有完整传代时间和细胞数，因此不进入下述累计时间计算。

### 2.1 时间尺度

| 实验谱系 | 有完整记录的传代次数 | 原始实验累计时间 |
|---|---:|---:|
| 2N control | 12 | 29 天 |
| 4N control | 12 | 29 天 |
| 2N O1 | 23 | 122 天 |
| 2N O2 | 23 | 122 天 |
| 4N O1 | 22 | 125 天 |
| 4N O2 | 22 | 121 天 |

单次传代间隔并不恒定：

- 取值为 1、2、3、4、5、6、7、9、10、12 或 14 天；
- 范围为 1–14 天；
- 中位数为 5 天。

现有拟合对象只保存每次传代的间隔天数，没有原始绝对日历日期。因此，本文件所说的“实验累计时间”是同一实验谱系各次 `passage_duration` 的总和，而不是由绝对日期重新计算的日历跨度。

### 2.2 细胞数量尺度

在 114 个完整传代记录中：

- 初始细胞数范围为 267,844–1,532,740；
- 终末细胞数范围为 1,605,234–13,953,612；
- 114/114 次传代的真实终末细胞数均大于同次初始细胞数；
- 在 108 个可连接的相邻传代中，108/108 次的前一传代真实终末细胞数均足以提供下一传代记录的接种细胞数。

因此，真实实验中的传代并不需要通过“提前收获”来满足下一次接种。记录的传代时间是实验操作时间尺度，而不是由模型根据细胞数反推的自由变量。

## 3. 三个时间概念

### 3.1 实验 `passage_duration`

`passage_duration` 是数据中记录的单次传代持续时间，单位为天。它描述细胞从本次接种到本次收获所实际经历的培养和 O₂ 暴露时间。

模型通过以下代码读取该字段：

- `oxygen/code/O2_supply_demand_MAP/util/o2_supply_demand_map_invitro_lineage_utils.R` 第 3–11 行。

### 3.2 模型 `duration_days`

`duration_days` 是模型为每个 passage segment 设置的最大模拟时间窗。

模型不是始终逐条使用原始 `passage_duration`。如果同一个模拟 segment 同时对应 O1 和 O2 记录，则代码对两条记录取中位数：

\[
duration\_days_s
=
\operatorname{median}
\left(
passage\_duration_{s,\mathrm{O1}},
passage\_duration_{s,\mathrm{O2}}
\right).
\]

`initial_cells` 和 `final_cells` 也采用同样的 segment 中位数汇总。随后模型建立从第 0 天到 `duration_days` 的逐日报告网格。

相关代码：

- `oxygen/code/O2_supply_demand_MAP/util/o2_supply_demand_map_invitro_objective_utils.R` 第 7–28 行；
- `oxygen/code/O2_supply_demand_MAP/util/o2_supply_demand_map_invitro_lineage_utils.R` 第 39–47、98–135 行；
- `oxygen/code/O2_supply_demand_MAP/util/o2_supply_demand_map_invitro_lineage_simulation_utils.R` 第 106–114 行。

底层数值模拟使用 `DT = 0.05` 天，但 `selected_day` 是从逐日输出点中选择的。`duration_days` 决定模拟最多运行到哪一天，本身并不保证该日状态会用于终点比较或下一次传代。

### 3.3 模型 `selected_day`

对每个 segment，模型只考虑大于 0 的模拟日，并选择预测活细胞数与目标细胞数绝对差最小的一天：

\[
selected\_day_s
=
\underset{t\in\{1,\ldots,duration\_days_s\}}{\arg\min}
\left|
N^{\mathrm{pred}}_{\mathrm{live},s}(t)
-
N^{\mathrm{target}}_s
\right|.
\]

目标值优先使用该 segment 汇总后的 `final_cells`；只有终末细胞数缺失时才回退到下一 segment 的 `initial_cells`。

这个选择规则有两个重要性质：

1. `selected_day` 不是实验记录的传代日，而是依赖当前模型参数和预测生长曲线的数据驱动选择；
2. 即使模型在整个 `duration_days` 时间窗内都没有达到目标细胞数，代码仍会选择“距离最近”的一天，不会判定该 segment 不可达。

相关代码：

- `oxygen/code/O2_supply_demand_MAP/util/o2_supply_demand_map_invitro_lineage_simulation_utils.R` 第 175–224、274–288 行。

## 4. `selected_day` 如何取代真实传代时间

### 4.1 下一次传代的状态来自 `selected_day`

模型提取 `selected_day` 当天的活细胞染色体组成。下一次 passage segment 开始时，这个组成被按下一段观察到的初始细胞总数重新缩放：

\[
\boldsymbol{x}_{s+1}(0)
=
\frac{\boldsymbol{x}_{s}(selected\_day_s)}
{\sum_N x_{s,N}(selected\_day_s)}
\times
N^{\mathrm{obs}}_{\mathrm{initial},s+1}.
\]

因此，虽然模型没有单独维护一个显式的全局日历，但链式状态实际经历的累计时间是：

\[
T_{\mathrm{effective}}
=
\sum_s selected\_day_s,
\]

而不是：

\[
T_{\mathrm{experiment}}
=
\sum_s passage\_duration_s.
\]

相关代码：

- `oxygen/code/O2_supply_demand_MAP/util/o2_supply_demand_map_invitro_lineage_simulation_utils.R` 第 238–258 行。

### 4.2 预测增长率也使用 `selected_day`

模型预测增长率为：

\[
\widehat g_s
=
\frac{
\log N^{\mathrm{pred}}_{\mathrm{live},s}(selected\_day_s)
-
\log N^{\mathrm{pred}}_{\mathrm{live},s}(0)
}{
selected\_day_s
}.
\]

分母不是实验 `passage_duration`，也不是模型 `duration_days`。因此，`selected_day` 同时决定：

- 用于增长拟合的时间长度；
- 用于染色体数和流式终点比较的状态；
- 传给下一次 passage 的染色体组成。

相关代码：

- `oxygen/code/O2_supply_demand_MAP/util/o2_supply_demand_map_invitro_summary_utils.R` 第 42–75 行。

## 5. 实验累计时间与模型累计 `selected_day`

缩短百分比按各实验谱系的真实累计时间计算：

\[
\text{缩短百分比}
=
\frac{
T_{\mathrm{experiment}}-T_{\mathrm{selected}}
}{
T_{\mathrm{experiment}}
}
\times 100\%.
\]

结果来源：

- `oxygen/results/fit_invitro_O2_buffering_500seed/seed10/invitro_lineage_summary.tsv`
- `oxygen/results/fit_invitro_O2_buffering_500seed/seed340/invitro_lineage_summary.tsv`

| Seed | 谱系 | 实验累计时间 | 模型累计 selected time | 缩短 | 缩短占实验周期 |
|---|---|---:|---:|---:|---:|
| seed10 | 2N control | 29 天 | 23 天 | 6 天 | 20.7% |
| seed10 | 4N control | 29 天 | 23 天 | 6 天 | 20.7% |
| seed10 | 2N O1/O2（每条线） | 122 天 | 120 天 | 2 天 | 1.6% |
| seed10 | 4N O1 | 125 天 | 106 天 | 19 天 | 15.2% |
| seed10 | 4N O2 | 121 天 | 106 天 | 15 天 | 12.4% |
| seed340 | 2N control | 29 天 | 29 天 | 0 天 | 0.0% |
| seed340 | 4N control | 29 天 | 28 天 | 1 天 | 3.4% |
| **seed340** | **2N O1/O2（每条线）** | **122 天** | **40 天** | **82 天** | **67.2%** |
| seed340 | 4N O1 | 125 天 | 111 天 | 14 天 | 11.2% |
| seed340 | 4N O2 | 121 天 | 111 天 | 10 天 | 8.3% |

用户原表中的八行数值均正确。为覆盖全部正式实验组，上表补入了 seed10 和 seed340 的 4N control。

### 5.1 严格区分模型 `duration_days` 与原始实验时间

2N O1 和 O2 每次传代的原始持续时间相同，因此两条线的累计 `passage_duration` 和模型累计 `duration_days` 都是 122 天。

4N O1 和 O2 有两个 passage 的持续时间不同：

- A15：O1 为 12 天，O2 为 10 天，模型取中位数 11 天；
- A19：O1 为 7 天，O2 为 5 天，模型取中位数 6 天。

所以，4N O1 和 O2 的原始实验累计时间分别为 125 天和 121 天，但模型对两条线共享的累计 `duration_days` 是 123 天。

严格比较模型内部的 `duration_days` 和 `selected_day`：

| Seed | 模型谱系 | 累计 duration_days | 累计 selected_day | selected 相对 duration 缩短 |
|---|---|---:|---:|---:|
| seed10 | 2N control | 29 天 | 23 天 | 6 天（20.7%） |
| seed10 | 4N control | 29 天 | 23 天 | 6 天（20.7%） |
| seed10 | 2N O1/O2 | 122 天 | 120 天 | 2 天（1.6%） |
| seed10 | 4N O1/O2（共享模型段） | 123 天 | 106 天 | 17 天（13.8%） |
| seed340 | 2N control | 29 天 | 29 天 | 0 天（0.0%） |
| seed340 | 4N control | 29 天 | 28 天 | 1 天（3.4%） |
| **seed340** | **2N O1/O2** | **122 天** | **40 天** | **82 天（67.2%）** |
| seed340 | 4N O1/O2（共享模型段） | 123 天 | 111 天 | 12 天（9.8%） |

因此，4N 总时间差应分解为两个步骤：

| Seed | 实验谱系 | 原始实验时间 | 中位数汇总后的 duration | selected time | 汇总造成的变化 | 选择造成的变化 | 总变化 |
|---|---|---:|---:|---:|---:|---:|---:|
| seed10 | 4N O1 | 125 | 123 | 106 | −2 | −17 | −19 天 |
| seed10 | 4N O2 | 121 | 123 | 106 | +2 | −17 | −15 天 |
| seed340 | 4N O1 | 125 | 123 | 111 | −2 | −12 | −14 天 |
| seed340 | 4N O2 | 121 | 123 | 111 | +2 | −12 | −10 天 |

这说明用户原表中的 19/15 天和 14/10 天是相对各条真实实验谱系的总差异；如果只讨论模型内部 `selected_day` 与 `duration_days` 的差异，则应使用 17 天或 12 天。

## 6. 直接导致的结果

### 6.1 相同传代次数对应不同累计生物学时间

模型仍模拟与实验相同次数的 passage，但每次从 `selected_day` 提取状态。因此，只要 `selected_day` 提前：

- 下一次传代也从更早的状态开始；
- 整条谱系的有效时间提前；
- 单位累计时间内发生的传代和稀释瓶颈次数被人为提高。

对 seed340 的 2N O1/O2，23 次低氧传代仍全部存在，但被压缩到 40 天，而实验为 122 天。模型只保留了实验累计时间的 32.8%。

### 6.2 低氧暴露和进化过程被压缩

模型中的细胞分裂、增长抑制、低氧相关死亡、染色体错误分离、WGD、错误分离后存活选择和染色体组成变化都随时间累积。当有效时间从 122 天缩短到 40 天时，模型用于解释终点 ploidy 的过程不再对应实验的完整低氧暴露历史。

因此，passage 编号相同不等于生物学时间相同。seed340 的终点分布是“经过 40 个模型日和 23 次传代瓶颈后的状态”，不能直接称为对 122 天实验终点的时间一致预测。

### 6.3 终点选择和增长拟合被同一个自由时间对齐规则耦合

`selected_day` 是根据预测细胞数与观察终末细胞数的距离选出的，然后同一天又被用于计算预测增长率和提取终点 ploidy。这样会形成参数依赖的时间重对齐：

- 如果模型增长轨迹不能在真实传代日解释终末细胞数，它仍可选择更早的一天；
- 这个更早的时间随后又进入增长率分母；
- 优化过程中可以在“终点组成更接近”与“增长预测恶化”之间进行非生物学的时间补偿。

这也是为什么不能把 `selected_day` 仅解释为输出或绘图层面的技术选择。

### 6.4 “最接近目标”不等于达到目标

以 passage index 34 的 2N、0% O₂ segment 为例：

- 实验/模型时间窗为 5 天；
- 目标终末细胞数约为 12.34 million；
- seed10 选择第 5 天，预测活细胞数约为 7.08 million，预测增长率为 0.371/day；
- seed340 选择第 1 天，预测活细胞数约为 1.05 million，预测增长率为 −0.051/day。

seed340 并没有接近实验终末规模；只是其预测细胞数随后继续下降，因此第 1 天成为整个候选时间窗中“距离最小”的一天。这个例子说明，提前的 `selected_day` 可能标志目标不可达，而不是合理的提前传代。

### 6.5 负增长下的强制传代和群体数量抬升

#### 6.5.1 问题如何产生

设上一 passage 在模型终点可用的活细胞状态为：

\[
\boldsymbol{x}^{\mathrm{available}}_s
=
\boldsymbol{x}_s(selected\_day_s),
\qquad
N^{\mathrm{available}}_s
=
\sum_N x^{\mathrm{available}}_{s,N}.
\]

下一 passage 记录的接种需求为：

\[
N^{\mathrm{required}}_{s+1}
=
N^{\mathrm{obs}}_{\mathrm{initial},s+1}.
\]

当前代码没有先判断：

\[
N^{\mathrm{available}}_s
\ge
N^{\mathrm{required}}_{s+1}.
\]

它只保留上一段的染色体组成比例，然后无条件把总量设置为下一段的 `initial_cells`：

\[
\boldsymbol{x}_{s+1}(0)
=
\frac{
\boldsymbol{x}^{\mathrm{available}}_s
}{
N^{\mathrm{available}}_s
}
\times
N^{\mathrm{required}}_{s+1}.
\]

当 \(N^{\mathrm{available}}_s > N^{\mathrm{required}}_{s+1}\) 时，这个操作可以解释为从足够大的收获群体中抽取一部分细胞进行正常接种。

但是，当 \(N^{\mathrm{available}}_s < N^{\mathrm{required}}_{s+1}\) 时，缩放因子为：

\[
J_s
=
\frac{
N^{\mathrm{required}}_{s+1}
}{
N^{\mathrm{available}}_s
}
>1.
\]

此时程序并不是传代或稀释，而是在 passage 边界把所有染色体状态同时乘以 \(J_s\)，从而凭空增加活细胞总数。

这个问题不是由一条显式的“负增长时增加细胞”代码触发的。程序对所有 passage 都采用相同的无条件缩放规则；负增长会降低 \(N^{\mathrm{available}}\)，因此特别容易把正常的向下稀释变成人为向上扩增。

产生该行为的代码路径为：

- `oxygen/code/O2_supply_demand_MAP/util/o2_supply_demand_map_invitro_lineage_simulation_utils.R` 第 247–258 行：用上一段 `selected_frac` 乘以当前段记录的 `initial_cells`；
- 同文件第 208–223 行：`selected_frac` 只保存相对组成，`reseeded_state` 直接按指定接种量重设总数；
- 当前可行性不进入目标函数；`n_growth_negative_pred` 只被统计，没有形成硬约束，见 `o2_supply_demand_map_invitro_objective_utils.R` 第 630–658 行。

#### 6.5.2 当前结果中的规模

真实实验中：

- 108/108 个相邻 passage 的前次真实终末细胞数都足以提供下一次真实接种量；
- 因此，不存在必须靠“向上缩放”才能完成下一次传代的实验记录。

以每条实验谱系的实际下一次接种量检查预测收获量：

- seed10：108/108 个相邻连接均足够；
- seed340：只有 85/108 个连接足够，23/108 个连接不足。

上述 23/108 是相对于 O1、O2 各自记录的下一次接种量进行的实验记录级检查。程序实际使用 O1/O2 segment 中位数接种量。按程序的共享 segment 结构计算：

- seed340 有 12 个独立模型边界发生 \(N^{\mathrm{available}} < N^{\mathrm{required}}\)；
- 其中 2N O1/O2 共享路径有 8 个边界，且这 8 个边界的预测增长率全部为负；
- 4N O1/O2 共享路径另有 4 个边界，其预测增长率为正，但收获量仍低于下一次中位数接种需求；
- 展开到 O1/O2 两条实验路径后，相当于 24/108 个路径连接会被程序向上缩放；
- 最大向上缩放因子约为 1.95，即 passage 边界处的群体数量可被直接提高约 95%。

因此，更准确的表述是：**无条件接种量重设在负增长的 2N 谱系中产生了系统性的强制传代和群体补充，同时也可在正增长但收获量不足的 passage 发生。**

#### 6.5.3 对科学结果的影响

1. **破坏细胞数量连续性。**  
   passage 前后的活细胞总数不再满足“从已有细胞中取样”的质量约束。模型可以从不足的母群体生成更大的下一代起始群体。

2. **把不可传代的轨迹强制延续。**  
   对 seed340 的 2N、0% O₂ 谱系，负增长本应导致后续接种失败或至少产生明显的拟合惩罚；当前实现却在边界补足细胞，使该谱系继续完成全部 passage。

3. **掩盖灭绝或衰退风险。**  
   一个相对高 ploidy 比例可以由低 ploidy 细胞更快死亡形成，但此时总活细胞数可能正在下降。边界补充细胞后，相对组成被保留、总量被恢复，模型便可持续传播这个高 ploidy 比例，从而把“相对比例上升”误读为“可持续群体扩张”。

4. **改变后续动力学。**  
   即使向上缩放不立即改变染色体组成比例，它会改变绝对群体规模和后续可达到的细胞数。在任何依赖绝对群体规模、拥挤或资源负荷的计算中，这个跳变还会继续影响下一段动力学。

5. **削弱增长和 ploidy 的联合约束。**  
   增长轨迹可以失败，但 ploidy 组成仍被强行送入下一次 passage。这样会人为减弱“一个机制必须同时产生足够细胞并给出正确组成”的科学约束。

6. **与 `selected_day` 提前形成反馈。**  
   负增长时模型倾向选择较早且细胞数相对较高的 `selected_day`；即便该日仍不足，下一段又会被强制补足。时间提前和数量补充共同让原本不可行的轨迹继续存在。

### 6.6 seed10 也不是完全没有时间偏差

逐条比较 114 个完整实验记录：

- seed10：87 次 `selected_day` 与该实验记录持续时间相同，27 次提前，0 次延后；
- seed340：53 次相同，59 次提前，2 次相对该条实验记录延后。

seed340 的两次“延后”来自 4N O2 的原始持续时间短于 O1/O2 中位数；模拟并没有超过其共享的模型 `duration_days`。

seed10 没有 seed340 那样严重的 2N 低氧时间崩溃，但其 2N/4N control 均压缩 20.7%，4N O1/O2 的共享模型时间压缩 13.8%。因此，seed10 可以比 seed340 更可靠，但不能据此认为当前时间处理已经严格正确。

### 6.7 seed340 的增长崩溃与时间错位相互关联

在 114 个有增长观察的传代中：

- seed10 没有负的预测增长率；
- seed340 有 34 个负的预测增长率，均出现在 2N、0% O₂ 条件。

seed340 的高终点 ploidy 不能脱离这一事实评价。它在 2N 低氧谱系中通过低倍体群体衰退、终点时间提前和状态重新缩放得到更高的高染色体比例，但没有同时再现实验的增长和完整暴露时间。因此，seed340 不能作为“改善终点 ploidy 且仍保持实验过程合理”的候选解。

## 7. 科学性判断

### 7.1 问题级别

这是模型结构和观测模型层面的时间对齐问题，不是变量命名、绘图或报告格式问题。

它直接影响：

- passage growth 的预测定义；
- 每次传代终点的染色体组成；
- 下一次传代的初始组成；
- 累计低氧暴露时间；
- CIN、WGD 和选择过程的累计机会；
- 终点 ploidy 的机制解释；
- 不同 seed 之间是否处于同一实验时间尺度的比较。

### 7.2 对当前结果的判断

- **seed340：不能用于主要科学结论。** 其 2N O1/O2 时间压缩达到 67.2%，并伴随负增长和 23 个相邻传代连接的预测接种量不足。其终点 ploidy 改善不是在保持实验时间和增长约束下获得的。
- **seed10：可作为当前结果中相对更合理的拟合基线，但仍需时间固定后的敏感性分析。** 特别是 control 和 4N 谱系的累计时间偏差不能忽略。
- **当前 endpoint ploidy 结果只能称为 passage-index-aligned prediction，不能称为 calendar-time-aligned prediction。**
- 在修正或敏感性分析完成前，不应使用当前链式轨迹对 CIN、WGD、低氧暴露累计效应或“达到某 ploidy 所需时间”作强机制结论。

## 8. 必须进行的修正或敏感性分析

### 8.1 首选方案：固定真实传代时间

对每个实验 passage：

1. 直接模拟到该条记录的 `passage_duration`；
2. 在真实终点提取细胞数量和染色体组成；
3. 用真实 `passage_duration` 计算预测增长率；
4. 把真实终点的染色体组成传给下一次传代；
5. 单独比较预测终末细胞数与观察终末细胞数；
6. 明确检查预测终末细胞数是否足以提供下一次接种量。

在这个方案中，细胞数拟合误差应由细胞数或增长率的观测模型承担，不能通过改变 passage 时间来消除。

### 8.2 O1/O2 应保留各自的实验时间

O1 和 O2 可以共享生物学参数，但不应因为共享参数而强制共享：

- `passage_duration`；
- 初始细胞数；
- 终末细胞数；
- passage 终点状态。

尤其是 4N A15 和 A19，应分别按 O1 和 O2 的真实时间模拟。中位数汇总会同时缩短 O1、延长 O2，并消除实验线之间已经记录的时间差异。

### 8.3 传代可行性门控：细胞不足时只记录 passage，不增加群体数量

建议采用用户提出的“无稀释保持”规则。设终点活细胞数为 \(N^{\mathrm{available}}_s\)，下一次正常接种所需的最低细胞数为 \(N^{\mathrm{required}}_{s+1}\)。

#### 情形一：细胞数足够

若：

\[
N^{\mathrm{available}}_s
\ge
N^{\mathrm{required}}_{s+1},
\]

则执行正常传代，从已有群体中抽取所需数量：

\[
\boldsymbol{x}_{s+1}(0)
=
\frac{
\boldsymbol{x}^{\mathrm{available}}_s
}{
N^{\mathrm{available}}_s
}
\times
N^{\mathrm{required}}_{s+1}.
\]

这一步只允许保持或降低总活细胞数，不允许增加。

#### 情形二：细胞数不足

若：

\[
N^{\mathrm{available}}_s
<
N^{\mathrm{required}}_{s+1},
\]

则：

1. passage index/实验时间段标记继续推进，以保持与观察表的 passage 顺序对齐；
2. 不执行稀释；
3. 不按下一次观察 `initial_cells` 放大群体；
4. 下一段直接继承当前完整活细胞状态：

\[
\boldsymbol{x}_{s+1}(0)
=
\boldsymbol{x}^{\mathrm{available}}_s;
\]

5. 记录该边界为未实际完成的传代：

\[
passage\_performed_s = \mathrm{FALSE};
\]

6. 记录接种缺口：

\[
shortfall_s
=
N^{\mathrm{required}}_{s+1}
-
N^{\mathrm{available}}_s.
\]

这就是“只做传代标记，不实际改变群体数量”的明确状态转移定义。它保证：

\[
\sum_N x_{s+1,N}(0)
\le
\sum_N x^{\mathrm{available}}_{s,N},
\]

因此 passage 操作本身永远不能创造活细胞。

#### 该规则还必须配合可行性标记或惩罚

保持群体数量可以修复细胞凭空增加的问题，但如果模型仍把这些边界当作正常传代用于完整拟合，优化器仍可能利用大量“未实际完成的 passage”解释后续数据。

由于真实实验 108/108 个相邻连接都具有足够接种细胞，合理拟合也应尽量满足这个条件。因此建议至少输出：

- `passage_performed`；
- `available_live_cells`；
- `required_inoculum_cells`；
- `reseed_shortfall_cells`；
- `reseed_scale_factor`；
- 每条谱系的失败传代次数和最大缺口。

对于正式拟合，可选择：

- 把任何 `passage_performed = FALSE` 视为不可行参数组合并给予硬惩罚；或
- 对缺口大小和失败次数加入显式似然/惩罚项。

对当前实验设计，硬约束更容易解释：实验确实完成了所有传代，因此一个无法产生最低接种量的模型轨迹不应被视为可接受的实验解释。

#### 必须通过的实现测试

1. 若 \(N^{\mathrm{available}}\ge N^{\mathrm{required}}\)，下一段总量严格等于 \(N^{\mathrm{required}}\)；
2. 若 \(N^{\mathrm{available}}<N^{\mathrm{required}}\)，下一段状态向量与上一段终点状态向量完全相同；
3. 不足分支的总细胞数不得在 passage 边界增加；
4. 不足分支仍正确推进 passage index 和指定 O₂ 时间段标记；
5. `passage_performed` 和 `shortfall` 必须出现在逐 passage 输出中；
6. 使用 seed340 回放时，应识别出上述 12 个独立不足边界，而不是静默补足；
7. 使用 seed10 回放时，正式实验路径不应触发不足规则。

### 8.4 如果坚持使用密度触发传代

如果模型假设实验是在达到某个细胞密度时才传代，则必须显式建模“达到阈值的时间”，并同时满足：

- 对观察 `passage_duration` 的拟合；
- 对终末细胞数的拟合；
- 目标细胞数必须在允许时间内可达；
- 收获细胞数必须足以提供下一次接种；
- 未达到阈值时应产生明确惩罚或判为不可行，而不是仍选择最近的一天。

在没有这些约束时，`selected_day` 不能替代真实传代时间。

## 9. 最终建议

在稿件继续使用体外模型的 endpoint ploidy 和机制结论之前，至少应完成一次“固定真实传代时间”的重算或敏感性分析，并报告：

- 每条谱系的真实累计时间、模型累计 `duration_days` 和累计 `selected_day`；
- passage 级预测终末细胞数、增长率和时间误差；
- 收获细胞数是否足以支持下一次接种；
- 每个 passage 边界是否实际执行传代，以及不足时是否保持群体数量而没有向上缩放；
- 每条谱系的传代失败次数、接种缺口和最大向上缩放风险；
- 固定时间后 seed10 的增长、ploidy、flow 拟合是否仍成立；
- 固定时间后 O1/O2 终点分布的拟合程度；
- 主要机制结论对时间处理方式是否稳健。

在完成上述分析前，最稳妥的表述是：当前模型结果证明某些参数组合能够在 passage 序列上产生相似的染色体组成变化，但尚未证明这些变化能够在真实累计培养时间和真实传代约束下同时发生。

## 10. 复核所用文件

### 实验元数据

- `/Users/4482173/Documents/GitHub/soft_coupling/oxygen/ploidyOxygen/data/fit_objects/fit_data.Rds`
- `/Users/4482173/Documents/GitHub/soft_coupling/oxygen/ploidyOxygen/data/fit_objects/jobs_2N.Rds`
- `/Users/4482173/Documents/GitHub/soft_coupling/oxygen/ploidyOxygen/data/fit_objects/jobs_4N.Rds`

### 模型代码

- `/Users/4482173/Documents/GitHub/soft_coupling/oxygen/code/O2_supply_demand_MAP/util/o2_supply_demand_map_invitro_lineage_utils.R`
- `/Users/4482173/Documents/GitHub/soft_coupling/oxygen/code/O2_supply_demand_MAP/util/o2_supply_demand_map_invitro_objective_utils.R`
- `/Users/4482173/Documents/GitHub/soft_coupling/oxygen/code/O2_supply_demand_MAP/util/o2_supply_demand_map_invitro_lineage_simulation_utils.R`
- `/Users/4482173/Documents/GitHub/soft_coupling/oxygen/code/O2_supply_demand_MAP/util/o2_supply_demand_map_invitro_summary_utils.R`

### seed 结果

- `/Users/4482173/Documents/GitHub/soft_coupling/oxygen/results/fit_invitro_O2_buffering_500seed/seed10/invitro_lineage_summary.tsv`
- `/Users/4482173/Documents/GitHub/soft_coupling/oxygen/results/fit_invitro_O2_buffering_500seed/seed10/invitro_growth_loglik.tsv`
- `/Users/4482173/Documents/GitHub/soft_coupling/oxygen/results/fit_invitro_O2_buffering_500seed/seed340/invitro_lineage_summary.tsv`
- `/Users/4482173/Documents/GitHub/soft_coupling/oxygen/results/fit_invitro_O2_buffering_500seed/seed340/invitro_growth_loglik.tsv`
