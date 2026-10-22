# Figure7/8 统一入口

本地与 HPC 使用相同当前代码和专用模型布局。
完整入口登记表：`Code/config/manifests/figure78_entrypoints.tsv`。

## 路径与模型规则

- 专用模型固定为工作区内 `Code/vendor/figure78_model`；12个源文件逐一验证锁定 SHA-256。
  不接受共享/外部模型覆盖，不修改 HPC 共享模型。
- 原始拟合输入按工作区相对路径定位到相邻 `soft_couping_org/oxygen/results` 的已批准
  joint-fit 目录；`--fit-root` 或 `FIGURE_JOINT_RESULT_ROOT` 可显式重定位原始拟合输入。
- 现有结果只从 `data/Figures/Figure7`、`Figure8` 读取，经验数据在
  `data/share_bundle_figures_3_5`。旧版本目录及归档不是运行输入或回退路径。
- 当前数据内的 provenance 冻结快照保留为来源证据，不执行为旧数值实现；不修改原始
  receipt、RDS 指纹。初始化函数身份认证改用专用锁定模型内的相同源码。

R 入口自行解析路径。Shell/HPC 入口自动读取 `figure78_hpc_paths.env`，不修改登录配置。
`figure78_Rscript.sh` 在本地使用 PATH 中的 Rscript，或显式指定的
`FIGURE78_NATIVE_RSCRIPT`；在 HPC 使用已批准的 R 4.4.2 Apptainer 镜像。
其他运行环境必须由 `FIGURE78_CONTAINER` 显式指定并另获批准。

## 入口对应

| 用途 | 当前入口 | 边界 |
|---|---|---|
| Figure7 完整生成 | `run_Figure7.R`、`data_Figure7.R` | 当前统一 prepare/task/aggregate/validate/plot |
| Figure8 完整生成 | `run_Figure8.R`、`data_Figure8.R` | 当前统一流程 |
| 保留数据绘图 | `draw_Figure7.R`、`draw_Figure8.R` 及单图包装 | 共26个产品，不模拟、不读归档 |
| 结果来源核验 | `validate_Figure78.R` | 只认证数值记录和身份；可能写 audit |
| 静态检查 | `figure78_entry_tools.py static` | AST、R parse-only、Shell语法、路径及哈希 |
| 本地/HPC 分派 | `Code/hpc/run_Figure78.sh` | 显式 figure、stage，不提交 |
| 清单中的一项任务 | `run_Figure78_task.sh` | 显式 prepared 根目录及真实 Slurm 索引 |
| 未来获批提交 | `submit_Figure78.sh` | 显式 `--submit=yes`、xxlarge、完整清单、无并发上限 |

现有保留数据不是新 prepared campaign；没有新任务清单时，worker/提交入口会拒绝。
完整生成需要另行批准的独立输出目录，或已有显式 rebuild 保护，不能隐式覆盖当前表。

旧 data/finalize/audit/HPC 文件名只保留兼容保护，不再调用旧 worker。生成别名必须显式
使用当前 `--stage=...`。旧 `--action`、任务编号、配置、位置参数和版本目录不能安全映射，
因此明确报错并提示当前接口，不猜测对应关系、不默认执行全流程。旧校验别名需要
`--mode=retained`，Python 别名也支持 `--mode=static`；旧 Slurm spool 包装拒绝旧数组索引。
旧实现与配置移入本次 audit 备份，不是当前运行依赖。

## 本次允许的静态命令

在 iteration5 下，使用现有 Python 3.6+ 与已批准的 R 环境：

```sh
python3 Code/Figures/figure78_entry_tools.py static --workspace=. --report=audit/figure78_static/checks.json
```

只写 audit 报告。R 文件仅 parse，不 source/eval；不执行 plan、prepare、模拟、绘图、
数值验收或 Slurm 命令。静态通过不代表数值/绘图运行验收通过。

## 需另行授权的运行

```sh
bash Code/hpc/run_Figure78.sh --figure=8 --stage=retained-validate
bash Code/hpc/run_Figure78.sh --figure=7 --stage=plot
```

前者认证现有记录并可能写 audit，后者重绘 Figure7 全部分析；本次静态验收不运行它们。
cold-start、独立数值验收和图形回归工具也会运行模型/出图，必须另行授权。
计算口径、科学配置、数值内核、随机流、绘图实现、样式和正式图片名保持不变。
