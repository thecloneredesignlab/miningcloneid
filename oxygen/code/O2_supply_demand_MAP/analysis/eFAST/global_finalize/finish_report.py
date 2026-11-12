#!/usr/bin/env python3
"""Write the final source-backed bilingual synthesis from verified exports."""
import argparse
from pathlib import Path
from common import read,load,save,now,table,sha


def finish(root):
    root=Path(root);p=load(root/'plan.json');latest=max(p['resolutions'])
    checks=load(root/'validation/independent_exports.json')
    if checks['status']!='passed':raise ValueError('Independent export verification required')
    diagnostics=read(root/'comparison/diagnostic_summary_disjoint_bands.tsv')
    scores=read(root/'comparison/method_difference_summary.tsv')
    sb=read(root/'sobol/oxygen_band_indices.tsv.gz')
    eb=read(root/'efast'/('N%d'%latest)/'oxygen_band_indices.tsv.gz')
    null=read(root/'comparison/structural_null_O2_zero.tsv')
    last_diagnostic={(r['output'],r['index']):r for r in diagnostics if r['method']=='eFAST' and int(r['N'])==latest and r['oxygen_band']=='all'}
    last_score={(r['output'],r['index']):r for r in scores if int(r['efast_N'])==latest}
    lines=['# 全域敏感性分析收尾：Sobol与eFAST','',
        '本轮计算、方法比较和数据核验已完成。Sobol N131072作为主要科学依据；eFAST首批四档完整保留，'
        '但最高N的倍性ST仍未稳定，不能宣称两种方法已经相互验证一致。','',
        '## 1. 完成范围','',
        '- 独立in-vivo固定氧模型；14个活跃参数，91个不重复二阶参数对。涉及p_mis_base、mu_hp、gamma_mu的焦点参数对共36对。',
        '- Sobol使用既有N131072、5次独立scramble；本轮没有重算Sobol模型或估计器。',
        '- eFAST使用SALib 1.5.2、M=4，N=%s；每档5次独立相位、完整201个氧点（0–5%%，每0.025%%）。'%(', '.join(map(str,p['resolutions']))),
        '- 两种方法的上下界、独立分布、变换、模型矩阵及数值恢复规则一致；603个缓存模型点与新求解路径逐值一致。',
        '- seed25仅提供既定配置和参考；全域14参数全部重新独立采样，不是seed25邻域，也不是500个拟合端点的后验抽样。',
        '- 生产计算共432,286,680个参数—氧点组合，两个输出同时计算；另有603个方法一致性验证点。四档2,380个模型分片和280个估计器分片均完成。',
        '- 邻域保持封存。本批没有N硬上限，完成的是已批准的四档；更高N可另建不可变批次追加。','',
        '## 2. 收敛诊断与方法比较','',
        '逐个输出×参数/参数对×氧点检查：五次重复最大值−最小值≤0.05，并且相邻N均值绝对差≤0.02。'
        '这些是操作性诊断，不是置信水平或交互显著性。所有单元保留，不根据是否通过筛选结果。',
        'Sobol的47,838个单元全部通过。eFAST首个新分辨率没有同协议的相邻较低N，标记为未评估相邻N诊断，'
        '不把零差值误记为通过；旧N513只作为历史参考。','',
        '| 方法/N | 输出 | 指标 | 通过/总数 | 通过率 | 最大重复范围 | 最大相邻N均值差 |',
        '|---|---|---|---|---|---|---|']
    for r in diagnostics:
        if r['oxygen_band']!='all':continue
        if r['method']=='Sobol' and r['index']=='S2':continue
        name='倍性' if r['output']=='dominant_mean_ploidy' else '增长率'
        rate='未评估' if not r['pass_fraction'] else '%.2f%%'%(100*float(r['pass_fraction']))
        count='—' if not r['n_passed'] else r['n_passed']+'/'+r['n_total']
        delta='—' if not r['resolution_delta_max'] else '%.5f'%float(r['resolution_delta_max'])
        lines.append('| %s %s | %s | %s | %s | %s | %.5f | %s |'%(r['method'],r['N'],name,r['index'],count,rate,float(r['repeat_range_max']),delta))
    lines += ['', '最高eFAST分辨率与Sobol的逐点差异：','',
        '| 输出 | 指标 | 平均绝对差 | 95分位绝对差 | 最大绝对差 |','|---|---|---|---|---|']
    for r in scores:
        if int(r['efast_N'])==latest:
            lines.append('| %s | %s | %.5f | %.5f | %.5f |'%('倍性' if r['output']=='dominant_mean_ploidy' else '增长率',
                r['index'],float(r['mean_absolute_difference']),float(r['p95_absolute_difference']),float(r['max_absolute_difference'])))
    lines += ['', 'S1的平均差较小，但仍存在未通过单元；不能仅凭相似颜色宣称收敛。'
        '倍性ST的重复波动和方法差异仍较大，增加N没有带来单调改善。增长率ST的方法差较小，但也没有全体通过。','',
        '### 结构性零效应检查','',
        '在O2=0、C=O2_crit>0、n=n_O>0时，h=C^n/(C^n+0^n)=1。当前模型中这两个参数只通过h进入方程，'
        '因此其理论S1/ST均为零。Sobol得到精确零；最高eFAST分辨率的ST如下：','',
        '| 输出 | 参数 | eFAST平均ST | 五次重复最小–最大 | 理论值 |','|---|---|---|---|---|']
    for out in ['dominant_mean_ploidy','dominant_growth_rate']:
        for param in ['O2_crit','n_O']:
            v=[float(r['estimate']) for r in null if r['method']=='eFAST' and int(r['N'])==latest and r['output']==out and r['parameter']==param and r['index']=='ST']
            lines.append('| %s | %s | %.5f | %.5f–%.5f | 0 |'%('倍性' if out=='dominant_mean_ploidy' else '增长率',param,sum(v)/len(v),min(v),max(v)))
    lines += ['', '这说明有限频谱估计仍会给结构上无作用的输入分配非零ST；不能解释为零氧下真实存在参数作用。'
        '该结果与频谱泄漏/混叠相容，但本轮没有改变M，未单独确定其具体来源或证明增大M能够消除。'
        '即使某个单元通过当前重复/分辨率阈值，也不代表不存在估计偏差。','',
        '## 3. 对Noemi假设的回答','',
        '对倍性，已收敛的Sobol支持“死亡/应激相关作用在低氧更重要，基线误分离与缓冲在高氧相对增强”的部分趋势。'
        '不支持把它简化为“低氧由死亡独占主导”或“所有误分离参数都在高氧增强”。lambda_max与死亡耦合的p_misseg同样重要。',
        '注意p_misseg是死亡/应激耦合的误分离放大参数，不能与基线p_mis_base混为一谈。'
        'mu_hp和gamma_mu也同时影响死亡和误分离，过程并非彼此独立。','',
        '| 倍性/ST参数 | Sobol低氧均值 | Sobol高氧均值 | eFAST低氧均值 | eFAST高氧均值 |',
        '|---|---|---|---|---|']
    parameter_rows=[]
    for param in ['p_mis_base','p_misseg','mu_hp','gamma_mu','buffer_smax','buffer_beta','buffer_n_exp','lam_max']:
        a={r['oxygen_band']:float(r['band_mean']) for r in sb if r['output']=='dominant_mean_ploidy' and r['index']=='ST' and r['parameter_i']==param}
        b={r['oxygen_band']:float(r['band_mean']) for r in eb if r['output']=='dominant_mean_ploidy' and r['index']=='ST' and r['parameter_i']==param}
        lines.append('| %s | %.5f | %.5f | %.5f | %.5f |'%(param,a['low'],a['high'],b['low'],b['high']))
        parameter_rows.append(dict(parameter=param,sobol_low=a['low'],sobol_high=a['high'],efast_low=b['low'],efast_high=b['high']))
    table(root/'comparison/ploidy_ST_band_comparison.tsv',parameter_rows)
    lines += ['', '这里的高低氧变化是**敏感性大小的氧偏好**，不是增加参数后倍性变化的正负方向。'
        '例如gamma_mu的ST氧偏好在两种方法间相反，且eFAST倍性ST尚未稳定，因此不能说三个焦点参数已经获得跨方法一致验证。',
        '按预先声明的过程内参数索引平均，Sobol低氧死亡ST均值0.36327，误分离0.36604，不能据此称死亡单独最大；'
        '缓冲均值从0.09053升至0.18508。这些平均值只是描述性统计，不是过程方差份额，ST不能相加。',
        '增长率是独立的补充输出：Sobol中lam_max的ST从0.47699升至0.70979，p_mis_base从0.16676升至0.27543；'
        '缓冲组均值反而从0.01726降至0.01281。因此增长率没有显示与倍性相同的缓冲增强，不能把两类输出混写。','',
        '## 4. 两两交互与响应面','',
        'S2来自Sobol，不来自eFAST。全部91对均已整理，36对焦点参数另列。按氧点峰值，'
        'lam_max×p_misseg在两种输出中最大：倍性0.10495（O2=0%），增长率0.23646（O2=0.075%）。',
        '低氧倍性的重要交互还包括p_misseg×mu_hp（峰值0.05665）、lam_max×mu_hp（0.03442）；'
        '高氧buffer_smax×buffer_n_exp增强（峰值0.02415，5%）。'
        '增长率中lam_max×p_mis_base在高氧突出（峰值0.13176，5%）。',
        '峰值排名与物理氧区积分排名分别提供，不相互替代。弱S2的细微名次不作为机制证据。'
        'ST−S1包含涉及该参数的所有交互阶数，不能把整个差额归给一个S2搭档。',
        '复用六对原有全域条件响应面及40,344个原始网格点，没有新增模型计算。其余参数固定在seed25背景，'
        '响应面只显示0、0.5、2.5和5%四个代表氧点；敏感性分析仍覆盖201点。'
        '等值线可帮助描述该背景下的补偿组合，参考点非加性对比的正负不等于Sobol S2的正负。',
        '没有对应响应面的强S2参数对，不能仅凭S2补写协同/拮抗方向；TII和高阶方向分析未在本轮开展。','',
        '## 5. 图怎么读','',
        '- 最终图组在publication/。Sobol汇总图左列倍性、右列增长率；各行是S1、ST和36对焦点S2。eFAST汇总图只有S1/ST。',
        '- 横轴为O2：0–0.025%采用线性段保留零点，其后为对数段；1%和3%虚线分隔低、中、高氧。原Figure4B参考图保留原始坐标。',
        '- 每个输出分别标色条。跨方法S1/ST图在同一输出、同一指标内使用相同色条上限及Sobol倍性行序。相同颜色跨不同色条不一定等值。',
        '- P为原Figure4B展示用过程分组；C为原相关强度的氧偏好，不能替代正负rho；S为当前方法的倍性敏感性氧偏好。'
        '增长率图仍沿用倍性S注释。P展示分组与描述性过程统计的primary_group对应表在process_group_membership.tsv。',
        '- S1/ST有P、C、S三列；S2只有P、C，分裂色块分别代表参数对两成员。排序按高氧≥低氧时高氧降序、其余低氧升序，不使用第三列决定顺序。',
        '- 曲线阴影是五次重复最小–最大，不是95%置信区间。负S1/S2估计值保留，不解释为负方差或拮抗。',
        '- S1/ST/S2不提供参数效应的正负方向；原Figure4B Spearman图提供拟合解之间的方向信息。','',
        '## 6. 使用边界与复现','',
        '本分析针对既定输入范围、独立分布和固定染色体状态域22–154。采样收敛不证明状态截断无偏差，也不证明参数可辨识。'
        '全域包含净增长率为负的组合，其倍性是主导特征向量组成，不能等同于正在稳定扩增的肿瘤。',
        'p_wgd范围为1e-5–1e-3，不含精确零；拟合接近下界或敏感性较低均不能证明WGD生物学缺失。',
        '数据、脚本、原始HPC导出和出版图分开保存。原HPC导出保持哈希不变；两类曲线图在publication/中仅调整标题/图例布局。',
        '最终诊断计数使用互斥氧区：低[0,1]、中(1,3)、高[3,5]；带平均值另按[0,1]、[1,3]、[3,5]在物理O2上积分。'
        '原HPC efast_diagnostics表的窗口端点有重叠，不应相加；最终计数以diagnostic_summary_disjoint_bands.tsv为准。',
        '复现入口见代码目录README；完整性记录见validation/independent_exports.json、method_alignment.tsv及PDF验证记录。'
        '第一份Sobol报告遇到模块名冲突，已用报告专用快照恢复；模型和估计器数组未重新投递，原失败与恢复证据均保留。']
    (root/'interpretation_zh.md').write_text('\n'.join(lines)+'\n')
    original=(root/'methods_results_en.md').read_text()
    extra=('\n## Additional robustness findings\n\n'
        'All approved FAST resolutions completed without filtering entries by convergence. At N=%d, %s/2814 ploidy ST entries '
        'and %s/2814 growth ST entries met the operational diagnostics. Mean absolute S1 differences from Sobol were '
        '%.5f for ploidy and %.5f for growth, whereas mean absolute ploidy ST difference was %.5f. '
        'At zero oxygen, O2_crit and n_O are structurally inactive: Sobol gave zero, while FAST mean ploidy ST estimates '
        'were %.5f and %.5f. Thus repetition/resolution checks alone do not exclude spectral estimation error. '
        'Cross-method agreement for ploidy total effects has not been established; gamma_mu also had opposite oxygen-preference '
        'patterns. Primary mechanistic interpretation therefore relies on the completed Sobol analysis, with FAST retained '
        'as a diagnostic comparison. The evidence supports relative shifts in selected death/stress, baseline missegregation '
        'and buffering parameters, but not exclusive death dominance or a uniform high-oxygen increase of all missegregation '
        'parameters. Sensitivity magnitudes do not identify effect direction or parameter identifiability.\n'%(
            latest,last_diagnostic['dominant_mean_ploidy','ST']['n_passed'],last_diagnostic['dominant_growth_rate','ST']['n_passed'],
            float(last_score['dominant_mean_ploidy','S1']['mean_absolute_difference']),float(last_score['dominant_growth_rate','S1']['mean_absolute_difference']),
            float(last_score['dominant_mean_ploidy','ST']['mean_absolute_difference']),
            *[sum(float(r['estimate']) for r in null if r['method']=='eFAST' and int(r['N'])==latest and r['output']=='dominant_mean_ploidy' and r['index']=='ST' and r['parameter']==param)/5 for param in ['O2_crit','n_O']]))
    (root/'publication/methods_results_en.md').write_text(original+extra)
    (root/'README.md').write_text(
        '# Completed global Sobol / eFAST comparison\n\n'
        '**Primary result:** existing Sobol N131072. **FAST batch:** N2049,4097,8193,16385, M4, five phases, 201 oxygen points. '
        'All computations completed; FAST ploidy ST remains unstable. Convergence is diagnostic-only. Neighborhoods remain archived.\n\n'
        '- [Chinese interpretation](interpretation_zh.md)\n'
        '- [English Methods / Results with robustness findings](publication/methods_results_en.md)\n'
        '- [Final figure package](publication/README.md)\n'
        '- [Sobol S1/ST/focal-S2 summary](publication/sobol/sensitivity_summary.pdf)\n'
        '- [Highest-N eFAST S1/ST summary](publication/efast/N16385/sensitivity_summary.pdf)\n'
        '- [Matched ST heatmaps](publication/comparison/method_comparison_st.pdf)\n'
        '- [Structural zero-effect check](publication/comparison/structural_null_O2_zero.pdf)\n'
        '- [Conditional response surfaces](response_surfaces/README_zh.md)\n'
        '- [Method alignment](method_alignment.tsv)\n'
        '- [Independent numerical verification](validation/independent_exports.json)\n'
        '- [Final nonoverlapping-band diagnostic counts](comparison/diagnostic_summary_disjoint_bands.tsv)\n'
        '- [Original transfer manifest](collection_manifest.json)\n\n'
        'Source tables and five-repeat indices are in sobol/, efast/N*/, comparison/ and inputs/. '
        'Original HPC exports remain unchanged. Publication curves correct title/oxygen-label and legend spacing without '
        'rerunning models or estimators. Runtime, source commits, scheduling and recovery receipts are retained.\n')
    save(root/'validation/report_synthesis.json',dict(status='complete',latest_efast_N=latest,source_verification_status=checks['status'],
         convergence_only_diagnostic=True,neighborhood_resumed=False,model_calls=0,estimator_calls=0,completed_utc=now()))
    manifest=load(root/'publication/manifest.json')
    manifest['files']=[dict(path=str(path.relative_to(root)),sha256=sha(path),size=path.stat().st_size)
        for path in sorted((root/'publication').rglob('*')) if path.is_file() and path.name!='manifest.json']
    save(root/'publication/manifest.json',manifest)
    print('Bilingual synthesis complete')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--root',type=Path,required=True)
    finish(parser.parse_args().root)
