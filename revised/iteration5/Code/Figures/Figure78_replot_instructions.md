# Figure7 / Figure8：统一代码与数据入口

当前范围说明（2026-10-17）：本文只说明重绘。完整模拟和数据生成现由 `run_Figure7.R` / `run_Figure8.R` 提供；`data_Figure7.R` / `data_Figure8.R` 已接入真实生成流程，不再是仅核验入口。见[完整流程说明](Figure78_pipeline_ZH.md)。

当前每个 Figure 只有一套作图实现和一个数据根目录，不使用日期或 Vxx 选择结果。
所有命令只读取已经验证的结果；不拟合、不模拟、不积分、不重新汇总 seed，也不读取 arxives。

## 运行

在 revised/iteration5 下运行；本机 Rscript 的完整路径为 /usr/local/bin/Rscript。

```sh
Rscript --vanilla Code/Figures/draw_Figure7.R
Rscript --vanilla Code/Figures/draw_Figure7.R --analysis=history
Rscript --vanilla Code/Figures/draw_Figure7.R --analysis=all
Rscript --vanilla Code/Figures/draw_Figure8.R
Rscript --vanilla Code/Figures/draw_Figure78_all.R
```

Figure7 默认绘制即时分析8张；history 为历史谱系8张；all 为两种分析共16张。
Figure8 包含主图及补图1–9，共10张。全套26张，4份矩阵PDF各8页，合计54页。
历史分析的正式输出名称仍为 figure7_ver2 / supp_fig7_*_ver2，以免影响论文引用；
它们是科学分析类型而非数据/代码版本。新入口使用 --analysis；--version 仅保留为旧命令兼容别名。

独立补图和旧分组脚本只作薄包装，全部调用相同实现，不加载旧计算模块。
仅绘主图可加 --products=assembled_fig7 或 --products=assembled_fig8；多个图名用逗号分隔。

```sh
Rscript --vanilla Code/Figures/draw_Figure78_all.R --output-root=/absolute/path/to/trial
Rscript --vanilla Code/Figures/draw_Supp_Figure8_8.R --formats=pdf
Rscript --vanilla Code/Figures/draw_Figure7.R --preflight-only=yes
Rscript --vanilla Code/Figures/draw_Figure8.R --preflight-only=yes
```

只核验而不模拟时，使用绘图入口的 --preflight-only=yes。
data_Figure7.R 和 data_Figure8.R 现已接入完整当前计算流程，不能再当作仅核验命令使用。
完整阶段命令见 Figure78_pipeline_ZH.md。原始冻结记录只核验保留的历史发布结果；
新生成数据使用自己的当前 provenance，不要求历史迁移或旧 worker 清单。

## 数据和代码结构

```text
data/Figures/Figure7/
  tables/       即时、历史谱系及来源诊断表
  provenance/   原始数值结果、完成/验证记录、冻结源码与配置
data/Figures/Figure8/
  tables/       主图及全部补图汇总表
  seed_data/    in vivo / in vitro 逐seed时间快照
  provenance/   参数传递、传代规则、完成/验证记录、冻结源码与配置
Code/Figures/util/analysis/
  figure7_plotting.R       Figure7唯一绘图实现
  figure8_plotting.R       Figure8唯一绘图实现
  figure78_plot_common.R   共同样式及读取辅助
  figure78_render_only.R   产品清单、统一配置与命令行
  figure78_render_plots.R  产品组装
  figure78_render_io.R     暂存、保护、发布与矩阵打包
  figure78_provenance.R    不可变来源核验
Code/config/figure7.json
Code/config/figure8.json
```

参数、模拟方式、传代过程、分子/分母、中值/IQR、筛选及图形风格均未改变。
即时flux矩阵和即时growth使用同一份当前即时表；旧flux参考表仅用于原始来源链核验。
旧固定p表按 pooled_fixed_low/high 命名，避免旧补图号与当前补图号混淆。

## 路径迁移与来源校验

Code/config/manifests 下：

- figure78_plot_contract.tsv：26个当前图的分析含义、呈现方式、区间定义和页数。
- figure78_current_plot_inputs.tsv：当前全部输入的语义组、逻辑名称、相对路径、大小、SHA-256。
- figure78_data_migration.tsv：原56个文件到新路径的对应关系，所有文件内容不变。
- figure78_frozen_sources.tsv：原始源码/配置路径与冻结副本的对应关系。
- figure78_bundle_bindings.tsv：原始RDS中代码/模型指纹的有序绑定。
- figure78_layout.json：仅用于解释原始记录中的路径，不是旧路径回退。

冻结副本只用于身份核验，不加载为绘图实现。原始source/code清单、RDS参数和指纹不改写，
新代码通过迁移映射验证原始记录；模型与Gaussian初始化函数仍与原始RDS指纹比较。
绘图无需原版本目录、归档、原拟合文件或外部模型源码目录。
当前还会验证 `Code/vendor/figure78_model` 的12个锁定源文件；初始化函数身份核验使用
该专用副本，不再加载数据内的冻结初始化源码，不执行模型模拟。见[统一入口说明](Figure78_entrypoints_ZH.md)。
每次预检验证所选Figure全部当前文件，任何缺失/损坏均失败，不自动更新指纹或重算。

## 输出与验收

默认更新 iteration5/Figures 和 manuscript/Figures 两份正式输出；加 --output-root 可隔离试绘。
先暂存，通过输入/源码及非目标图片保护后才发布，并备份被更新的目标图。
每次运行的审计、缓存和备份写入输出根目录 audit/figure78_render_only，不写入科学数据目录。

其他选项：--preflight-only=yes、--list-products=yes、--formats=pdf,png、--data-root、--audit-root。
未知选项/图名、冲突的分析选择，以及指向 arxives 的路径或软链接均拒绝。
Python需有pypdf；FIGURE78_PYTHON可以指定已有解释器，不自动安装或提交HPC。

验收入口为 validate_Figure78_render_only.py；覆盖全套/单图、PDF页数、PNG一致性、镜像、
输入缺失/损坏、非目标图保护、统一目录和不重跑模拟。迁移审计见
audit/figure78_canonical_consolidation_20261016。
