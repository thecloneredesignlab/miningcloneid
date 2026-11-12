# 全域eFAST与Sobol最终比较

Sobol固定为N131072；新eFAST分辨率为[2049, 4097, 8193, 16385]，每档5次独立相位、201个固定氧点。没有预设最大N。
新eFAST与Sobol使用相同参数范围/分布、冻结固定氧矩阵和异常恢复；方法一致性检查见validation/method_alignment.json。
旧N513仅作为历史参考；旧元数据未完整保存best_params与fit_config哈希，不能赋予与新运行同等的身份认证。

## S1与ST差异

- Ploidy / S1：平均绝对差0.00248，95分位绝对差0.01005，最大绝对差0.01558。
- Ploidy / ST：平均绝对差0.08118，95分位绝对差0.16469，最大绝对差0.20238。
- Net live-cell growth rate / S1：平均绝对差0.00219，95分位绝对差0.00914，最大绝对差0.02207。
- Net live-cell growth rate / ST：平均绝对差0.01001，95分位绝对差0.03534，最大绝对差0.04615。

## Noemi假设的方法一致程度

- Ploidy：最高新eFAST分辨率与Sobol在三项ST描述性检验中，2/3项对比方向一致。
  - low_O2_death_largest_descriptive_group_mean：Sobol对比值-0.00277，eFAST对比值-0.06473；各自支持重复数0/5、1/5。
  - missegregation_high_above_low：Sobol对比值-0.02749，eFAST对比值+0.05837；各自支持重复数0/5、5/5。
  - buffering_high_above_low：Sobol对比值+0.09455，eFAST对比值+0.13198；各自支持重复数5/5、5/5。
- Net live-cell growth rate：最高新eFAST分辨率与Sobol在三项ST描述性检验中，3/3项对比方向一致。
  - low_O2_death_largest_descriptive_group_mean：Sobol对比值-0.09348，eFAST对比值-0.10722；各自支持重复数0/5、0/5。
  - missegregation_high_above_low：Sobol对比值-0.06575，eFAST对比值-0.06565；各自支持重复数0/5、0/5。
  - buffering_high_above_low：Sobol对比值-0.00445，eFAST对比值-0.00169；各自支持重复数0/5、2/5。

## 科学解释与使用边界

Noemi假设的Sobol逐项解释见../sobol/interpretation_zh.md；各eFAST分辨率的同口径过程描述和假设检查分别保存在efast/N*/源表中。
S1/ST比较评价方法稳健性；S2仍只来自Sobol。不得把eFAST的ST−S1标为S2。
收敛仅作诊断，所有参数、氧点和重复保留；未通过不能被隐藏为缺失。通过也不等于交互显著或生物学机制被证实。
五次重复范围不是95%置信区间。微弱S2的细微排名不作为机制证据。
过程均值不等于正式组Sobol指数；死亡参数也参与误分离通路。响应面固定背景，不是全域交互方向。
本轮只完成批准的四档；更高N可作为后续独立不可变批次追加，无N=8193硬上限，也未无限自动提交未指定的计算。
邻域保持封存。

全部逐点差值、重复范围及未通过项见本目录源表。
