#!/usr/bin/env python3
"""Write a bilingual synthesis from verified results, without fixed FAST conclusions."""
import argparse
from pathlib import Path
from common import read,load,save,now,table,sha,resolutions,efast_directory


def oxygen_preference_comparison(sobol,fast):
    rows=[]
    for out in ['dominant_mean_ploidy','dominant_growth_rate']:
        for param in sorted({r['parameter_i'] for r in fast if r['index']=='ST'}):
            def bands(data):
                return {r['oxygen_band']:float(r['band_mean']) for r in data
                        if r['output']==out and r['index']=='ST' and r['parameter_i']==param}
            a,b=bands(sobol),bands(fast)
            ad,bd=a['high']-a['low'],b['high']-b['low']
            rows.append(dict(output=out,parameter=param,sobol_low=a['low'],sobol_high=a['high'],
                efast_low=b['low'],efast_high=b['high'],sobol_high_minus_low=ad,efast_high_minus_low=bd,
                same_trend=(ad>0)==(bd>0) and (ad<0)==(bd<0)))
    return rows


def finish(root):
    root=Path(root);p=load(root/'plan.json');latest=max(p['resolutions'])
    checks=load(root/'validation/independent_exports.json')
    if checks['status']!='passed':raise ValueError('Independent export verification required')
    diagnostics=read(root/'comparison/diagnostic_summary_disjoint_bands.tsv')
    scores=read(root/'comparison/method_difference_summary.tsv')
    sb=read(root/'sobol/oxygen_band_indices.tsv.gz')
    eb=read(efast_directory(root,p,latest)/'oxygen_band_indices.tsv.gz')
    null=read(root/'comparison/structural_null_O2_zero.tsv')
    latest_d=[r for r in diagnostics if r['method']=='eFAST' and int(r['N'])==latest and r['oxygen_band']=='all']
    ploidy_st=next(r for r in latest_d if r['output']=='dominant_mean_ploidy' and r['index']=='ST')
    passed=int(ploidy_st['n_passed']);total=int(ploidy_st['n_total'])
    state='全部通过操作性诊断' if passed==total else '仍有未通过操作性诊断的单元'
    alignment=load(root/'validation/method_alignment.json')
    prep=load(root/'preparation.json')
    band_rows=oxygen_preference_comparison(sb,eb)
    table(root/'comparison/ST_band_comparison.tsv',band_rows)
    table(root/'comparison/ploidy_ST_band_comparison.tsv',[r for r in band_rows if r['output']=='dominant_mean_ploidy'])
    focus=[r for r in band_rows if r['output']=='dominant_mean_ploidy' and r['parameter'] in ['p_mis_base','mu_hp','gamma_mu']]
    opposite=[r['parameter'] for r in focus if not r['same_trend']]
    lines=['# 全域敏感性分析：Sobol与eFAST追加结果','',
        'Sobol继续使用既有N131072。eFAST最高N=%d，倍性ST%s（%d/%d通过）；所有结果均保留。'%(latest,state,passed,total),
        '这比较的是各方法当前完成的设计，不是相同计算预算下的效率试验；N的定义和模型调用次数另列。','',
        '## 完成范围','',
        '- 独立in-vivo固定氧模型；14个参数，主导平均倍性和渐近净活细胞增长率两个输出。',
        '- 新eFAST档位：%s；已认证只读参考档位：%s。'%(', '.join(map(str,p['resolutions'])),', '.join(map(str,p.get('reference_resolutions',[]))) or '无'),
        '- SALib 1.5.2，M=4，各5次独立相位，201个氧点（0–5%，步长0.025%）。',
        '- 本轮新设计共%d个参数—氧点位置，%d个模型分片和%d个估计器分片；两个输出同时计算。'%(p['model_points'],p['model_tasks'],70*len(p['resolutions'])),
        '- 实际额外一致性验证点%d；复用既有验证点%d。'%(prep.get('model_alignment_validation_points',0),prep.get('reused_model_alignment_validation_points',0)),
        '- 参数范围、独立分布、模型和求解规则一致，校验见method_alignment.tsv和validation/。',
        '- Sobol模型及估计器不重算；邻域保持封存。','',
        '## 收敛诊断','',
        '逐单元要求五次重复最大值减最小值≤0.05，并且相邻N均值绝对差≤0.02。'
        '通过率是数值诊断，不是置信水平或生物学显著性；未通过项不排除。'
        '追加档位承接已认证的前一档，不把N32769误记为没有相邻参考。','',
        '| 方法/N | 输出 | 指标 | 通过/总数 | 通过率 | 最大重复范围 | 最大相邻差 |',
        '|---|---|---|---|---|---|---|']
    for r in diagnostics:
        if r['oxygen_band']!='all' or (r['method']=='Sobol' and r['index']=='S2'):continue
        assessed=r['diagnostic_status']=='assessed'
        count=r['n_passed']+'/'+r['n_total'] if assessed else '未评估'
        rate='%.2f%%'%(100*float(r['pass_fraction'])) if assessed else '—'
        delta='%.5f'%float(r['resolution_delta_max']) if assessed else '—'
        lines.append('| %s %s | %s | %s | %s | %s | %.5f | %s |'%(r['method'],r['N'],r['output'],r['index'],count,rate,float(r['repeat_range_max']),delta))
    lines+=['','## 与既有Sobol结果的差异','',
        '| 输出 | 指标 | 平均绝对差 | 95分位绝对差 | 最大绝对差 |','|---|---|---|---|---|']
    for r in scores:
        if int(r['efast_N'])==latest:
            lines.append('| %s | %s | %.5f | %.5f | %.5f |'%(r['output'],r['index'],float(r['mean_absolute_difference']),float(r['p95_absolute_difference']),float(r['max_absolute_difference'])))
    lines+=['','重复范围、相邻差和跨方法差值分别报告。即使通过操作性诊断，也不能仅凭相似颜色宣称无偏或一致。','',
        '## 结构性零效应','',
        '零氧时h=C^n/(C^n+0^n)=1，因此O2_crit和n_O均无作用，理论S1/ST为零。',
        '| 输出 | 参数 | 最高N的eFAST平均ST | 五次重复最小–最大 |','|---|---|---|---|']
    null_means={}
    for out in ['dominant_mean_ploidy','dominant_growth_rate']:
        for param in ['O2_crit','n_O']:
            v=[float(r['estimate']) for r in null if r['method']=='eFAST' and int(r['N'])==latest and r['output']==out and r['parameter']==param and r['index']=='ST']
            if len(v)!=5:raise ValueError('Missing five-repeat structural null check')
            null_means[out,param]=sum(v)/5
            lines.append('| %s | %s | %.6f | %.6f–%.6f |'%(out,param,sum(v)/5,min(v),max(v)))
    lines+=['','非零估计是相对理论零值的数值误差。本轮固定M=4；不能仅靠这些结果认定其具体来源，或保证改变M能消除误差。','',
        '## Noemi假设与焦点参数','',
        '主要科学依据仍为已完成的全域Sobol，详细过程检验和91对S2排名见sobol/。'
        '它支持部分低氧死亡/应激与高氧基线误分离、缓冲增强的趋势；死亡并非独占主导，整个误分离组也并非一致增强。'
        '增长率单独解释，不套用倍性的缓冲结论。',
        '以下是倍性ST的高低氧偏好，不能解释为增大参数后倍性变化的正负方向：','',
        '| 参数 | Sobol低氧 | Sobol高氧 | eFAST低氧 | eFAST高氧 | 趋势相同 |','|---|---|---|---|---|---|']
    for r in band_rows:
        if r['output']=='dominant_mean_ploidy':
            lines.append('| %s | %.5f | %.5f | %.5f | %.5f | %s |'%(r['parameter'],r['sobol_low'],r['sobol_high'],r['efast_low'],r['efast_high'],'是' if r['same_trend'] else '否'))
    lines+=['','三个焦点参数中，ST高低氧均值趋势不一致的参数为：%s。'%(', '.join(opposite) or '无'),
        '趋势相同只代表带均值方向一致，仍须结合相应单元诊断与效应差值；不等于参数效应符号一致。',
        '过程参数索引平均是描述性统计，不是组方差份额，ST不能相加。S2只来自既有Sobol；ST−S1包含所有阶数的交互。','',
        '## 图与使用范围','',
        'P/C/S注释、倍性行序、高氧降序与低氧升序规则沿用既定图形。S2只有P/C。'
        'S注释在增长率图中仍表示倍性氧偏好。O2保留0–0.025%的线性段，其后对数显示，1%和3%虚线为氧区分界。',
        '同一输出、同一指标的跨方法面板共用色条；倍性和增长率各自定标。五次重复阴影为最小–最大，不是95%置信区间。',
        '条件响应面固定其余参数背景，不能代表全域交互方向。本轮不新增响应面模型计算；既有全域响应面继续有效。',
        '输入独立分布、固定染色体状态域与负增长组合等使用边界保持不变。采样诊断不证明参数可辨识或WGD生物学缺失。']
    (root/'interpretation_zh.md').write_text('\n'.join(lines)+'\n')
    english=('\n## Refinement diagnostics\n\nAt eFAST N=%d, %d/%d ploidy ST entries met the operational diagnostics. '
        'All entries were retained. The three focal inputs with differing high-minus-low ploidy ST band-mean signs were: %s. '
        'This concerns oxygen preference of sensitivity magnitude, not the sign of a parameter effect. '
        'The mean structural-null ST errors at zero oxygen were %.6f for O2_crit and %.6f for n_O (ploidy). '
        'Matched N values are not matched computational budgets. Numerical stability and cross-method agreement must be assessed separately. '
        'Sobol models and estimators were not rerun; endpoint neighborhoods remained archived.\n'%(latest,passed,total,
         ', '.join(opposite) or 'none',null_means['dominant_mean_ploidy','O2_crit'],null_means['dominant_mean_ploidy','n_O']))
    (root/'publication').mkdir(exist_ok=True)
    (root/'publication/methods_results_en.md').write_text((root/'methods_results_en.md').read_text()+english)
    (root/'summary_for_noemi_en.md').write_text('# Global sensitivity refinement\n\n'+english)
    (root/'README.md').write_text('# Completed global eFAST refinement\n\n'
        'New FAST resolutions: %s. Five phases; 201 oxygen levels; M4. Existing Sobol N131072 is read-only. '
        'Latest ploidy ST diagnostics: %d/%d passed. Neighborhoods remain archived.\n\n'%(', '.join(map(str,p['resolutions'])),passed,total)+
        '- [Chinese interpretation](interpretation_zh.md)\n- [English methods and results](publication/methods_results_en.md)\n'
        '- [Sobol summary](publication/sobol/sensitivity_summary.pdf)\n'
        '- [Latest FAST summary](publication/efast/N%d/sensitivity_summary.pdf)\n'%latest+
        '- [Matched ST heatmaps](publication/comparison/method_comparison_st.pdf)\n'
        '- [Independent verification](validation/independent_exports.json)\n'
        '- [Diagnostic counts](comparison/diagnostic_summary_disjoint_bands.tsv)\n')
    save(root/'validation/report_synthesis.json',dict(status='complete',latest_efast_N=latest,
        source_verification_status=checks['status'],convergence_only_diagnostic=True,
        neighborhood_resumed=False,model_calls=0,estimator_calls=0,completed_utc=now()))
    manifest=load(root/'publication/manifest.json')
    manifest['files']=[dict(path=str(path.relative_to(root)),sha256=sha(path),size=path.stat().st_size)
        for path in sorted((root/'publication').rglob('*')) if path.is_file() and path.name!='manifest.json']
    save(root/'publication/manifest.json',manifest)
    print('Bilingual synthesis complete from actual resolution-specific source tables')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--root',type=Path,required=True)
    finish(parser.parse_args().root)
