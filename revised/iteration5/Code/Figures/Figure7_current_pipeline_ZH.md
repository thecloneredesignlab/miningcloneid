# Figure7 当前正式流程

统一入口是 `run_Figure7.R`；`data_Figure7.R` 调用同一入口。专用调度为
`util/analysis/figure7_current_pipeline.R`，培养数值实现为已独立验收的
`figure7_passaged.R/.cpp`。Figure8与历史 `_ver2` 不走新计算或写入。

## 参数、协议和统计量

使用C01 seed94、C02 seed80完整joint-fitting winners：in vivo26个参数，in vitro33个。
培养配置仅来自 `fit_result.rds$ctx$invitro_cfg`。固定原生p_misseg=.1/.3、场景氧、
纯44/88条根起始及设计阻断；不把有效错误分离概率固定为同一个常数。
不修改专用12文件锁定模型或HPC共享模型。

control/O1/O2从原实验DAG继承组成，各自循环末5个完整实测元组。完整窗口选取
最接近目标的正观察日，再按下一元组实际接种量重标定。不使用统一接种量、
阈值触发随机抽样、接种不足停止、N<1停止或未达目标停止。
累计1000天包含记录期；最后一段先完整选日再截断。重接种可条件性放大期望数量，
不应解释为真实库存守恒。完整科学定义见 manuscript/figure_legends/Figure7_calculation_methods.md。

A/B分别是同一日存活低子细胞的高母/低母来源，分母均为HL+LL；有产出时互补，
day0或零产出均NA。逐路线计算后取三路线中值和type-7 IQR，表示协议差异而非CI。
supp_fig7_8保留非遗传第一代D和永久历史谱系两种诊断，速率采用当前低倍性状态。

## 完整运行

必须提供独立输出目录，不能用prepare/rebuild隐式覆盖正式数据：

```sh
Rscript Code/Figures/run_Figure7.R --stage=plan --profile=full
Rscript Code/Figures/run_Figure7.R --stage=prepare --profile=full --data-root=/absolute/campaign/data
Rscript Code/Figures/run_Figure7.R --stage=simulate --profile=full --data-root=/absolute/campaign/data --task-index=48
Rscript Code/Figures/run_Figure7.R --stage=aggregate --profile=full --data-root=/absolute/campaign/data
Rscript Code/Figures/run_Figure7.R --stage=validate --profile=full --data-root=/absolute/campaign/data
Rscript Code/Figures/run_Figure7.R --stage=plot --profile=full --data-root=/absolute/campaign/data --output-root=/absolute/campaign/rendered
```

320个条件扩展为660个worker：170个培养条件各137/139/141N，加150个连续条件。
各worker内部计算三个β，培养每β三路线；每task有原子RDS和SHA收据，每路线也有
SHA及输入/源码指纹。缓存不符立即拒绝，成功结果可复用。中断留下无收据文件时，
先移动到可恢复audit位置，不覆盖也不当作成功。缓存与模型副本按worker物理隔离。

完整网格门槛为两次相邻比较都通过：份额≤.5个百分点，平均ploidy≤.01N，
log数量≤.01，CDF≤.01，上界丢失产出比例≤.001，传代日程完全一致。
网格失败不发布、不自动扩大上界。数值完成不等于图片已发布。

## HPC

```sh
bash Code/hpc/submit_Figure7_current.sh --data-root=/absolute/prepared/campaign/data --submit=yes
```

red partition、xxlarge QoS；完整array范围不附加 `%` 并发限制。依赖后处理仅执行
aggregate/validate，不自动绘图或覆盖正式数据。失败/缺失项可用显式 `--task-ids=...`
补跑；重提交前检查并保留已完成的收据。

## 验收与发布

`test_Figure7_current.R` 做正式清单、原生协议、独立Euler及选定任务检查；
`figure78_entry_tools.py static --source-only=yes` 仅做源码/入口/锁定模型/原始输入检查，
明确不认证旧正式数据与新代码的数值一致性。最终普通static和render preflight仍要求
真实完成收据、数值门槛、表和源文件SHA匹配，不能用source-only代替。

只生成当前Figure7及supp_fig7_1–8；保持版式和样式。先独立重绘并检查PDF/PNG，
全部通过后按白名单更新当前表、来源收据和图片。历史表保留原始认证链，不使用
新campaign指纹为旧数据背书。Figure8和所有 `_ver2` 数据、来源文件、图片保持原字节。

`validate_Figure7_current_figures.py` 仅检查独立重绘的当前版成品，不模拟、不绘图。
检查9个当前产品及16个矩阵分图：每套镜像50个PDF/PNG文件、共86个PDF页面，
并导出逐页栅格供人工检查。新数值故意改变，不能要求它们与旧图片像素一致；
应保持原尺寸/页数/版式，再验证数据绑定、文字和可读性。

`figure7_passaged.json` 的小规模入口授权标记仅约束独立验收程序；正式调度复用其
科学协议与网格门槛，不借该小规模入口运行全量。正式计算需显式独立campaign，
且必须有本次获批执行与备份记录。
