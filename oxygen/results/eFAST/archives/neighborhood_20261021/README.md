# eFAST 与 Sobol 邻域分析阶段封存 — 2026-10-21

用户已明确暂停 **eFAST 和 Sobol 的全部邻域分析**。恢复前需要用户重新明确提出要求。
封存入口为本目录；恢复步骤见 [RESUME.md](RESUME.md)。全域分析不属于此次暂停范围。

## 保存了什么

- eFAST：500 个拟合端点的邻域主分析，分辨率试运行，相位/频谱诊断，
  50/100 位特征对复核，以及 166 点的染色体状态边界复核。
- Sobol：初始邻域 pilot、分辨率 pilot、完整 N4096 邻域结果，以及 N8192
  定向细化的已完成模型输出、未完成检查点、失败记录和调度依赖。
- 两端代码快照：`local_code.tar.gz`、`hpc_code.tar.gz`，包含分析代码、
  模型、工具、固定氧实现及 Git 中的 container 构建文件；每个成员均有 SHA256。
  两端 Git HEAD 与快照内容分别登记，不假定两端代码完全相同。
  差异清单为 `code_snapshot_comparison.json`；共同的计算代码没有内容差异，
  共同路径的差异仅为边界复核 README，另有本地新增的解释/重绘和封存文件。
- `hpc_fit_inputs.tar.gz`：500 个 **独立 in-vivo** 拟合端点的 `best_params.tsv`、
  `fit_config.rds`、`parameter_table.csv`、`fit_summary.tsv`。前三者逐个与既有
  `seed_manifest.tsv` 的 SHA256 核对；不是 in-vitro 或 joint fitting 输入。
- 参数上下界、分布变换、氧浓度网格、采样设计、数值恢复记录、缓存模板、
  原始输出及结果表/图的路径清单。共享 Sobol 输入、冻结代码和缓存记录作为恢复依赖保留；
  这些共享文件的登记不表示暂停全域分析。
- Slurm 逐任务终态与退出码、取消回执、失败任务到 seed/重复/采样块的映射。
- 固定 SALib 1.5.2 SIF 的实际完整 SHA256 校验记录 `hpc_runtime.json`。

## 封存时的真实阶段

| 分析阶段 | 保存的状态 |
|---|---|
| eFAST 500 邻域主分析 | 原 `status.json` 为 complete：35,000 个轨迹任务、2,500 个设计；收敛诊断独立保留 |
| eFAST 后续试运行与数值/边界诊断 | 保存现有结果及各自核验，不继续扩大 N、M 或 N_MAX |
| Sobol 完整 N4096 | 原 `status.json` 为 complete，500 个邻域均已运行；具体 undefined/数值/收敛状态以原表为准 |
| Sobol N8192 定向邻域细化 | **未完成**：74,518 个模型任务 COMPLETED/0:0，32 个 FAILED/127:0；失败的 32 个块均无最终输出、完成回执或检查点 |
| Sobol N8192 后处理 | 2,485 个分析任务及报告原处于等待，现已按用户要求取消；不生成缺失细化结果的完整汇总 |

取消的作业父 ID：`21140548`、`21140549`，及其延时监控 `21140562`、`21140563`。
取消后已检查队列及 accounting；保留的其他用户作业未被取消。未重新投递任何分析任务。
完整证据见 `evidence/hpc_stop_receipt.json`、`hpc_job_accounting.json` 和
`hpc_job_accounting.tsv.gz`。失败块清单见 `sobol_refinement_resume_state.json`。

原有完成状态、索引、收敛判定和解释未被重写。封存不代表所有估计已收敛；
166 点的 308→440 边界稳定也不等于重新计算了完整轨迹的 S1/ST。
此前 N_MAX=154 的正式索引仍受该状态域限制。

## 数据保存方式与校验范围

1. **本地独立备份**：`payload/local_results.tar`，34,054 个文件，
   9,644,707,840 字节。复制及读回时逐文件核对 SHA256；整包校验见 `local_payload.json`。
   原始本地文件也保留在原处。压缩的逐文件清单为 `local_results.jsonl.gz`；
   `evidence/local_verification.json` 记录原文件的再次完整内容核验。
2. **HPC 大型原始数据原位封存**：原结果、检查点和缓存继续保存在原绝对路径，
   由 `hpc_results.jsonl.gz` 登记。没有复制整套 HPC 原始数据或 SIF，也没有删除/改名原目录。
   共享盘接近满额，仅 eFAST 邻域目录即约 258 GB；本地包不是 HPC 全部原始数据的副本。
3. **明确区分内容校验和文件登记**：HPC 对小于等于 256 KiB 的文件、JSON、Markdown、
   冻结代码、汇总表和图生成新 SHA256；其余大型轨迹/NPZ/缓存包登记大小和 mtime。
   每行 `integrity` 为 `fresh_sha256` 或 `size_mtime_only`，不把历史校验声称为本次完整重算。
   后者在恢复复用前必须通过原完成回执中的内容哈希和设计身份核验。
4. **代码与输入真正复制**：两端代码和 500 端点输入均另存独立 tar.gz，逐成员读回核验。
   其 tar.gz 校验及成员清单在 `*_manifest.json`。SIF 保留原路径，仅核验整文件哈希。
5. 这是带清单和恢复约定的阶段封存，不是文件系统不可写快照。
   `NEIGHBORHOOD_ARCHIVED.json` 及分析代码目录的 `AGENTS.md` 记录暂停约束。

Git 保存本目录的说明、清单、代码/输入快照和封存回执；9.64 GB 的 `payload/`
由 `.gitignore` 排除，**不上传 GitHub**。从 Git 重新 clone 不会自动取得这个大包或 HPC 原始输出。
迁移或清理原目录前，应先把这些数据复制到有容量的独立存储并核验，不能仅保留 Git。

## 原路径

- 本地仓库：`/Users/4482173/Documents/GitHub/soft_couping_eFAST`
- HPC 仓库：`/share/lab_crd/taoli/Project/soft_couping_eFAST`
- eFAST：`oxygen/results/eFAST/neighborhood10pct`
- Sobol 完整邻域：`oxygen/results/eFAST/sobol/full_N4096/neighborhood10pct`
- Sobol 新尾部：`oxygen/results/eFAST/sobol/neighborhood_refinement`
- Sobol 细化控制/冻结代码：`oxygen/results/eFAST/sobol/refinement_20261018`
- 原拟合：`/share/lab_crd/taoli/Project/soft_couping_org/oxygen/results/fit_invivo_unified_np256_500seed_all_xxlarge_r442_exact_20260828_145253`
- SIF：`/share/lab_crd/taoli/Docker/o2_supply_demand_map_r442_hpc_exact_salib152_20260923.sif`
- SIF SHA256：`0b60f6cdab8a91f6bbeea1ab6cf02dd79660a295f24bc6e8ad6c5fed04d982ad`

完整机器可读入口为 `archive.json`。不要因原控制文件仍写着 submitted/running 就自动续跑；
本次封存的终态以取消回执与 accounting 为依据。

此次只同步新增的封存文件/约定/核验工具，没有覆盖 HPC 原科学文件或改变其 Git HEAD。
HPC 原 HEAD 为 `9166a2ddeccbd1d0e04405c4e4fb3f5172a12b9e`，
本地封存前 HEAD 为 `b1b7fc69e6137738a06c4c4a49ec754adaa6e94c`。
封存说明与元数据提交在当前分支，并用 `archive/neighborhood-20261021` 标签定位。
