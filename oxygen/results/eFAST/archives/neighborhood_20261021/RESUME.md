# 将来恢复邻域分析

**当前不执行本页的计算或投递步骤。** 只有用户明确要求恢复 eFAST 或 Sobol 邻域后，
才能解除对应的暂停约束。先确认恢复哪一套及哪一阶段；全域分析不需要随邻域重跑。

## 1. 验证封存与依赖

从当前仓库根目录，可只做读取校验：

```bash
python3 oxygen/code/O2_supply_demand_MAP/analysis/eFAST/archive_neighborhood.py verify \
  --repo . \
  --manifest oxygen/results/eFAST/archives/neighborhood_20261021/local_results.jsonl.gz
```

HPC 上改为同目录下的 `hpc_results.jsonl.gz`，并使用 HPC 仓库作为 `--repo`。
该检查分别报告 SHA256 与大小/mtime 验证数。大小/mtime 相同不足以认证内容。
即使文件清单检查通过，也必须沿用原始每个轨迹/NPZ 完成回执和设计身份核验。

- 按 `archive.json` 核验快照与清单本身的 SHA256。
- 两端原 Git HEAD 不同，具体见封存 README 与 `code_snapshot_comparison.json`。
  HPC 有未跟踪的原始科学输出；以后更新分支前先核对并保全文件，不能用 `git clean`
  或强制覆盖消除冲突。优先按原阶段记录的冻结代码恢复计算。
- 按 `local_payload.json` 核验本地 tar。若需还原，先在独立临时目录解包并核对成员哈希；
  不直接覆盖现在的工作目录。HPC 的完整大数据并不包含在此 tar 中。
- 原拟合目录如已丢失，`hpc_fit_inputs.tar.gz` 可还原全部 500 个端点的四类输入；
  原始路径与参数/O2/代码等哈希仍须核对。绝对路径变更需要新的路径/身份审计，不能直接改旧计划。
- 校验 `hpc_runtime.json` 指定的原 SIF，使用原 SALib 1.5.2 环境。
- 保留冻结代码、Rcpp 模板及已有绑定方式；不能重编译后默认为与原数值实现等价。
- 维持邻域定义：最佳值 ±0.10×完整自然参数跨度，裁剪到全局界限，继承 uniform/log-uniform 变换。
  保留原五次重复、14 参数及氧网格；收敛是诊断，不是删除 seed 的条件。

## 2. eFAST 恢复入口

eFAST 500 个邻域的主分析已经完成，不需要重新投递 35,000 个轨迹。
优先复用 `neighborhood10pct/summaries/`、主图和原始输出。

- 主工作流：`analysis/eFAST/neighborhood.py`、`neighborhood_slurm.py`。
- 分辨率：`resolution_pilot.py`、`resolution_4097.py`、`resolution_nm.py` 及各自提交回执。
- 诊断：`phase_spectrum.py`、`high_precision_review/`、`boundary_review/`。
- 原配置、N/M/seed/phase、输入身份、任务粒度和冻结代码都以结果目录中的计划与回执为准。
- 需要新的科学计算时创建独立的新阶段，保留本次封存的旧结果和诊断，不静默覆盖。
  若扩大 N_MAX，完整 S1/ST 需要在同一状态域重算相应完整设计；不能只用 166 个点修补索引。

## 3. Sobol 恢复入口

完整 N4096 结果保留原样。N8192 是对 497 个 seed、55,089 个 seed/O2 组合、六对参数
和五次重复的**定向细化**；不是 500 个 seed 在全部 O2 下统一达到 N8192。
`selection.json`、`targets.tsv.gz`、`plan.json` 不能为恢复操作重新选择或静默改写。

恢复前逐阶段处理：

1. 阅读 `sobol_refinement_resume_state.json`、`sobol_failed_model_task_ids.txt`，
   结合 `hpc_job_accounting.tsv.gz` 和原日志定位 32 个 FAILED/127:0 的任务。
   此次确认它们没有最终输出/完成回执；检查点是否存在逐任务登记。
   其余 74,518 个任务的 Slurm COMPLETED 不替代输出内容及 identity 认证。
2. 验证每个已有尾部及原 N4096 前缀的完成回执、哈希、样本顺序、氧索引和五个 scramble。
   未通过时保留证据并定位原因，不能假定只缺文件或悄悄重算整个 array。
3. 在用户授权的恢复阶段，仅对确认缺失/损坏、确实需要模型复算的最小任务集合做 smoke/恢复。
   复用已经认证的检查点。原 `plan.json` 的内容哈希参与 identity，不能为了改路径随意编辑。
   新的投递/恢复记录另存；已有完整块保持原字节。任何新目录映射需要显式身份审计。
4. 认证全部所需尾部后，才恢复原 2,485 个统计任务以及后续邻域报告。
   入口是冻结的 `refinement_20261018/slurm/code_snapshot/run_refinement.sbatch`，
   由 `SOBOL_REFINEMENT_ROLE`、`SOBOL_REFINEMENT_SCOPE=neighborhood10pct`、
   `SOBOL_REFINEMENT_N=8192` 和准确 task ID 定位。当前被取消的旧 job ID 不可直接“恢复”。
5. **不要直接运行整套 `submit_refinement.py`**：它包含全域作业，且会重复完整生产计划。
   根据已认证的缺口单独准备恢复提交清单，先确认资源/队列要求，新的延时监控也只对应新作业。
6. 任何混合分辨率汇总必须保存 `N_used`，注明哪些细胞/参数对用 N4096、哪些用 N8192；
   不把定向细化称作完整均匀 N8192，保留 undefined 和未通过收敛诊断的估计。

原有 pilot/conditional response surfaces 保留，不因恢复敏感性统计而重复调用模型。

## 4. 恢复后的记录

新增恢复说明、精确的任务/输入映射、认证结果、代码和环境标识、逐任务终态、结果收集记录。
如需要继续原位写入尚缺块，先明确记录解除封存的范围，并另存新版本清单，不能覆盖本目录的封存清单。
只有用户授权的对应范围才能解除 `NEIGHBORHOOD_ARCHIVED.json` 与 `AGENTS.md` 的暂停标记。
