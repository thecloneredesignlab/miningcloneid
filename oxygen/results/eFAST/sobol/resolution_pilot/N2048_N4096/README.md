# 二阶 Sobol 分辨率扩展试运行

2026-10-13：按确认的下一阶段运行 N=2048、4096，保留 N=1024 作比较。
本目录与上一轮 N256/N512/N1024 试运行分开。

## 分析范围

- 全域：完整 independent in-vivo 参数范围，五次随机化 Sobol 重复。
- 邻域：seed25、464、165、2、214、138、66、395 八个代表端点；
  最佳值 ±10% 全局自然参数跨度，裁剪到全局边界；每个端点五次重复。
- 固定氧条件：0、0.5、2.5、5%。
- 14 个有效参数、全部 91 对；同时计算主特征向量的平均倍性和主特征值的净活细胞增长率。
- 使用原 SIF、SALib 1.5.2、相同模型方程、参数分布和 scramble seeds。

本轮完整设计对应 22,118,400 次算子求值。旧 N1024 的 5,529,600 次求值
在核验模型代码、输入、嵌套设计及输出校验和后复用；新增 16,588,800 次。
N4096 每次完整二阶设计有 30×4096=122,880 个参数向量。
N 是基础采样分辨率，不是染色体数或拟合 seed 数。

## 运行与诊断

代码提交：`42dc4c08b46e09ebb8b204bc09c0b76786238c8c`。
Slurm 配置为 xxlarge、12 小时，不指定节点、不设 array 并发上限。
六段任务链为运行验证 → 1350 个增量求值块 → 九个范围汇总 → 指数报告
→ 54 个关键配对条件响应面 → 最终报告。

逐单元格诊断比较 N2048→N4096 均值差（≤0.02）和五重复范围（≤0.05）。
保留所有八个端点及原始负估计；这些标准是运行诊断，不构成数学收敛证明。
`resolution_comparison.tsv` 将对齐三种 N 的 S1、ST、S2 和重复范围。
`reused_indices_validation.json` 将检查本轮 N1024 指数是否精确重现旧结果。

实际作业号见 `slurm/submission_jobs.tsv`，准备参数及校验和见
`slurm/submission_plan.json`。只有 HPC `status.json` 为 complete 且最终
报告作业 COMPLETED/0:0 时，才可收集并称本轮分析完成。

HPC 输出根：
`/share/lab_crd/taoli/Project/soft_couping_eFAST/oxygen/results/eFAST/sobol/resolution_pilot/N2048_N4096`。
运行与收集说明见 [代码说明](../../../../../code/O2_supply_demand_MAP/analysis/eFAST/sobol/README.md)。

本轮仍是代表性试运行，尚未扩展为 500 个邻域或 201 个氧点的正式分析。
