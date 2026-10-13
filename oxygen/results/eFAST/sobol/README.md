# In-vivo 二阶 Sobol 试运行

完成日期：2026-10-13。最终报告任务 `20836987` 为 `COMPLETED / 0:0`。

## 分析范围

- 全域：in-vivo 完整参数上下界，1 个全局设计，5 次独立随机化重复。
- 邻域：8 个代表性拟合端点，分别使用最佳值 ± 全局自然参数跨度的 10%，截断到全局上下界；每个端点 5 次重复。
- 14 个有效参数、全部 91 个参数对，N=256、512、1024。
- O2=0、0.5、2.5、5%；两个输出为主特征向量对应的平均倍性和主特征值对应的净活细胞增长率。
- 共 1,350 个主采样矩阵块、9 个范围汇总、6 个关键配对 × 9 个范围的 54 个响应面。
- 这是代表性试运行，未开展全部 500 个端点或完整 201 点氧网格的 Sobol 分析。

## 阅读顺序

1. [参数与方程关系](model_structure/model_equations.md)、[逐参数映射](model_structure/parameter_equation_map.tsv)、[91 个配对的机制映射](model_structure/parameter_pair_mechanism_map.tsv)。
2. [全域 S2 热图](global/pilot/figures/s2_pair_by_oxygen.pdf)与[邻域 S2 热图](neighborhood10pct/pilot/figures/s2_pair_by_oxygen.pdf)。
3. [关键配对清单](key_pairs.tsv)及每个范围 `figures/responses/` 下的响应面。
4. [试运行解读与限制](analysis_notes.md)、各范围的 `interpretation.md`。
5. `indices_repeats.tsv.gz` 保留逐重复 S1/ST/S2；`indices_summary.tsv` 保留逐范围、逐采样量的诊断。邻域汇总是先平均每个端点的五次重复，再取八个端点的中位数。

## 质量与来源

7 项数值测试及 9 个配置的参考输出检查通过。主采样 5,529,600 次算子求值中，1 个点采用原矩阵的 50/100 位精度复核；最大特征对残差约 3.31e-14。敏感性收敛失败不会删除端点或配对。

响应面准备曾因最佳值与对数转换后的边界相差约 1e-16 而停止。原始采样和指标表保留，修正仅对响应网格采用 64 个机器 epsilon 的边界容差；原始计划与后处理修订见 `slurm/`。

[收集记录](artifact_collection_receipt.json)记录本地复制文件的 SHA-256 验证及最终 Slurm 状态。原始矩阵块、采样数组、恢复矩阵、日志和完整 Figure 4 外部参考表留在 HPC；本地保留源表、诊断、图与来源记录。

复现入口：[代码与运行说明](../../../code/O2_supply_demand_MAP/analysis/eFAST/sobol/README.md)。
