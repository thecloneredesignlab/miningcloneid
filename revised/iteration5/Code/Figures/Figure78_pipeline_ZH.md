# Figure7、Figure8：当前完整计算流程

本文是从原始完整拟合参数开始的执行说明，不是仅利用保留数据重绘。
英文对应：[English guide](Figure78_pipeline.md)。各图详细公式仍见中英文 calculation methods 文档。

## 入口与组织

- `run_Figure7.R`：即时/历史两种口径的主图及补图1–7。
- `run_Figure8.R`：主图及补图1–9。
- `data_Figure7.R`、`data_Figure8.R`：调用同一套真实数据生成流程，不再是仅核验入口。
- `draw_Figure7.R`、`draw_Figure8.R`：只绘图，不隐式启动模拟。

各图只有一个当前配置（`Code/config/figure7.json`、`figure8.json`）、一个数据目录
（`data/Figures/Figure7`、`Figure8`）。公共模型和传代模块共用，不复制为多套版本。
已有 `_ver2` 图片名保留，表示历史来源这一科学分析口径，并非代码/数据版本目录。

## 原始输入与参数

`Code/config/figure78_inputs.json` 明确 joint-fitting 根目录、两家族目录、Figure7 已批准的
C01 seed94/C02 seed80，以及当前工作区的专用锁定模型 `Code/vendor/figure78_model`。
联合拟合本身作为上游输入，不由本流程重新拟合；缺少原始输入直接报错，不回退到旧预测。

in vivo 读取各 seed 的26个 `shared_invivo` 参数及其配置；in vitro 读取33个
`invitro_effective` 参数和该 seed 的 `fit_result.rds$ctx$invitro_cfg`。
逐参数核对传递过程，只允许情景明确指定的氧浓度和原生 `p_misseg` 覆盖。
不固定有效错误分离概率，不改写 `p_mis_base`、`k_o_mis`，不按重复拟合终点去重。
Figure8 正式规模为每家族500 seeds。

模型文件由 `Code/config/manifests/figure78_model_sources.tsv` 锁定 SHA256。
准备时将核对过的当前模型闭包复制到 `prepared/model_source` 后编译，不改上游仓库，
也不调用旧 Figure7/8 worker。原始拟合位置可通过 `--fit-root` 或相应环境变量指定。
模型必须是当前工作区内的专用副本；`--model-root` 或 `FIGURE_MODEL_CODE_ROOT` 如提供，
必须解析到该副本，不再接受共享/外部模型。默认路径在本地与 HPC 按工作区相对位置解析。
完整入口说明见 [统一入口](Figure78_entrypoints_ZH.md)。

经验数据直接从 `data/share_bundle_figures_3_5` 的两份 Figure3 样本表及19行
`organ_lut.csv` 重新计算。传代规则直接从 C01 seed1 原始 fit result 的实验 adapters
恢复六条 lineage，再做 lineage 等权几何汇总：接种569799，阈值6180783。
不读取旧 Figure6/7/8 的规则表或模拟结果。

## 分阶段运行

在 iteration5 目录中：

```sh
Rscript Code/Figures/run_Figure7.R --stage=plan
Rscript Code/Figures/run_Figure8.R --stage=plan
```

`plan` 只检查输入与模型，报告任务清单规模，不模拟、不提交任务。
正式规模为 Figure7 的320项任务（各任务内部包含三种 β），以及 Figure8 的9160项
seed-block 任务；实际清单为准。

在独立目录生成完整数据示例：

```sh
Rscript Code/Figures/run_Figure7.R --stage=prepare --data-root=/absolute/new/data
Rscript Code/Figures/run_Figure7.R --stage=simulate --data-root=/absolute/new/data
Rscript Code/Figures/run_Figure7.R --stage=aggregate --data-root=/absolute/new/data
Rscript Code/Figures/run_Figure7.R --stage=validate --data-root=/absolute/new/data
Rscript Code/Figures/run_Figure7.R --stage=plot --data-root=/absolute/new/data --output-root=/absolute/new/rendered
```

Figure8 替换入口名称即可。`--stage=all` 真实执行完整流程，可能需要大量计算。
`simulate --task-index=N` 执行清单中的一项；`--task-ids=N,M` 指定多项；不指定则顺序执行全部。
兼容的已完成任务可以复用；Figure8 在任务内部逐 seed 保存检查点；缺任务禁止完整汇总。

已有正式/保留数据不会被准备步骤隐式覆盖。确实需要替换时，显式使用 `prepare --rebuild=yes`，
先将原目录移动到可恢复 audit 备份，再在同一个无版本号目录生成；不是删除，也不是新增 Vxx。
单纯检查或重绘时不要使用此选项。

## 保持不变的科学口径

Figure7 当前代码自行生成总群体轨迹，再在同一轨迹和相同传代样本上计算两种被动来源标签。
即时 A 为 `100 F_HL/(F_HL+F_LL)`；即时 B 为 `100 sum(J_LL D_L)/sum(J_LL X_L)`。
D 只包括当前存活且尚未再次分裂的直接高→低子细胞，标签不传给后代。
历史 A 为当前低群体中的永久来源谱系数量占比，历史 B 为该谱系当前的低细胞增殖贡献。
各 growth 占比均不扣死亡。Figure7 在两种环境中都保留纯2N/4N起始这一诊断设计。

Figure8 是有限时间预测，不用主特征向量代替终点。in vivo 使用稳定的有限时间矩阵指数。
in vitro 使用每个 seed 的拟合高斯起始、真实 DT 和密度配置，以及阈值触发随机传代：
整数天更新后检查真实群体量、系统随机整数化、超几何无放回抽样。
不是固定每5天稀释，也不改算成连续增长。

培养结果先在 seed 内对传代重复求均值，再在500或合并1000个 seed 之间取 type-7 中值、Q25、Q75。
混合起始按细胞状态/数量组合，不平均终点 ploidy。十个时间快照来自同一连续推进的1000天轨迹。
只有保守上界证明不可能触发传代时，才将200条相同的确定性重复折算为一条 Euler 轨迹，
同时记录执行重复数和名义重复数。

## 统一数据目录

```text
Figure7/ 或 Figure8/
  prepared/    完整 seed/context 参数、任务清单、bundle、当前模型闭包
  tasks/       当前任务结果和可继续的检查点
  seed_data/   Figure8 逐 seed 终点及时间快照
  tables/      主图和全部补图的绘图表
  provenance/  原始输入哈希、参数审计、当前代码哈希、原始传代推导、运行环境和数值校验
  .cache/      编译缓存，不是旧科学结果
```

新生成数据建立自己的 provenance，不伪装成旧 worker 的产物。
数值代码、配置或输入指纹改变后，禁止混用旧检查点。新数据绘图只核对当前生成记录，
不依赖历史整理清单或归档代码。旧保留数据的独立兼容绘图路径仍可使用，图形构建器和26个
主图/补图入口的风格保持不变。

## 代表性从零验收

```sh
Rscript Code/Figures/validate_Figure78_cold_start.R --data-root=/absolute/empty/test/data
Rscript Code/Figures/validate_Figure78_pipeline.R --data-root=/absolute/empty/test/data
```

第一步使用同一生产 worker、完整1000天，但只用2次传代重复；Figure8 每家族取2个 seeds。
R 层 readRDS/fread/source/sys.source 保护与日志禁止读取归档、旧冻结 worker 和保留的
Figure7/8 科学结果。计算保护启动前，只对正式文件做字节哈希以检查是否被改动，不将参考数值传给模拟。
入口自行生成保护基线，两条验收命令不依赖人工准备的 audit 文件。
可以复用编译缓存，但目录中不能已有科学 bundle、任务结果或汇总表。
第二步独立检查指标、传代、混合和分位数，再与保留结果做只读回归比较。
代表性测试数据不会被允许冒充500+500 seeds 的正式图片输入。

## HPC 边界

`Code/hpc/run_Figure78_task.sh` 通过已准备的清单及 `SLURM_ARRAY_TASK_ID` 执行同一 worker，
自身不提交任务。以后批准全量提交时使用 `xxlarge`，array 索引覆盖实际清单，不添加 `%N` 并发上限。
本机代表性验收不等于已在集群执行或完成全量任务。
