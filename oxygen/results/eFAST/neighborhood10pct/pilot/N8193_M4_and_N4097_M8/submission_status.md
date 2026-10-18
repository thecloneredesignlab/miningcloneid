# 本轮投递记录

计算代码提交：`9cea6b22feee20cae9dea2cc64467c8b2199bc05`，已 push，并同步至 HPC。

| 配置 | seed | smoke | 主 array（140） | seed 汇总（2） | 最终比较 |
|---|---|---:|---:|---:|---:|
| N8193_M4 | seed25, seed138 | 21138912 | 21138913 | 21138914 | 21138915 |
| N4097_M8 | seed101, seed174 | 21138916 | 21138917 | 21138918 | 21138919 |

## 核验与当前状态

- 指定 SIF 内 9 项测试全部通过，无跳过；Python 3.9.25、SALib 1.5.2。
- 调度器实际设置已核验：全部 xxlarge、12 小时，无指定节点、无 array 并发限制。
- 共 280 个主轨迹任务；五次重复、14 参数、201 个 O2 点，保持原有邻域和变换。
- 两个 smoke 在调度器快照中为 RUNNING；主 array 等待对应 smoke 的 afterok 依赖。
- 汇总和最终比较任务已提交；缺失或失败输入会明确报错，不会被标记为分析完成。
- 本地已核对 19 个源文件的 SHA-256，并将 450,240 个基线数值逐项与历史结果比较，一致。

## 比较边界

- seed25/138：N4097/M4 与 N8193/M4，检查 N 分辨率收敛。
- seed101/174：N4097/M4 与 N4097/M8，检查 M 设置稳健性；不能据此宣称 N 收敛。
- 五次相位范围不是置信区间；记录完整 O2 网格和低/高氧窗口的 p90、最大偏差及位置。
- 重点检查 ploidy ST，保留全部参数、S1 和 growth 的配套结果。

当前记录是成功投递与输入核验，不是数值收敛或分析完成报告。
调度器状态时间及实际字段见各配置的 `slurm/submission_validation.json` 和 `slurm/squeue_snapshot.tsv`。

本地结果根目录：`oxygen/results/eFAST/neighborhood10pct/pilot/N8193_M4_and_N4097_M8/`。
HPC 对应根目录：`/share/lab_crd/taoli/Project/soft_couping_eFAST/oxygen/results/eFAST/neighborhood10pct/pilot/N8193_M4_and_N4097_M8/`。
