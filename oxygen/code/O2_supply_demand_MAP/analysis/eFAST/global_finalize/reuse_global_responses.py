#!/usr/bin/env python3
"""Verify and copy the six existing global conditional response surfaces."""
import argparse
import math
from pathlib import Path
import shutil

from common import read,load,save,sha,now,table


def reuse(source,root):
    source,root=Path(source),Path(root)
    original=load(root/'inputs/source_sobol_submission_plan.json')
    case=original['cases'][0]
    rows=[r for r in read(source/'response_surface_inventory.tsv') if r['scope']=='global']
    if len(rows)!=6:raise ValueError('Expected the six existing global surfaces')
    dest=root/'response_surfaces';dest.mkdir(exist_ok=True)
    files=[];summary=[]
    for r in rows:
        directory=source/r['part'];meta=load(directory/'metadata.json')
        if sha(directory/'metadata.json')!=r['metadata_sha256']:raise ValueError('Response metadata differs')
        if meta['scope']!='global' or meta['fit_seed']!='seed25' or meta['method']!='Sobol_response_surface':
            raise ValueError('Unexpected response scope')
        if meta['best_params_sha256']!=case['best_params_sha256'] or meta['fit_config_sha256']!=case['fit_config_sha256']:
            raise ValueError('Conditional response background differs from the certified model')
        data=read(directory/'response_surface_source.tsv.gz')
        expected={(round(float(x),14),round(float(y),14),round(float(o),12))
                  for x in meta['grid_i'] for y in meta['grid_j'] for o in meta['oxygen_pct']}
        actual={(round(float(x['value_i']),14),round(float(x['value_j']),14),round(float(x['O2_pct']),12)) for x in data}
        if len(data)!=int(r['n_points']) or len(data)!=len(actual) or actual!=expected:
            raise ValueError('Missing or duplicate response grid points')
        columns=['dominant_mean_ploidy','dominant_growth_rate','ploidy_nonadditive_contrast','growth_nonadditive_contrast']
        if not all(math.isfinite(float(x[c])) for x in data for c in columns):raise ValueError('Nonfinite response data')
        name=directory.name
        pairs=[(directory/'metadata.json',dest/name/'metadata.json'),
               (directory/'response_surface_source.tsv.gz',dest/name/'response_surface_source.tsv.gz')]
        for suffix in ['pdf','png']:
            pairs.append((source/'global/full/global/figures/responses'/(name+'.'+suffix),dest/(name+'.'+suffix)))
        for src,target in pairs:
            target.parent.mkdir(parents=True,exist_ok=True)
            if target.exists() and sha(target)!=sha(src):raise ValueError('Refusing to replace a differing surface export')
            if not target.exists():shutil.copyfile(str(src),str(target))
            if sha(target)!=sha(src):raise ValueError('Response readback differs')
            files.append(dict(source=str(src.resolve()),path=str(target.relative_to(root)),sha256=sha(target),size=target.stat().st_size))
        for o2 in meta['oxygen_pct']:
            selected=[x for x in data if abs(float(x['O2_pct'])-o2)<1e-12]
            for c in columns:
                values=[float(x[c]) for x in selected]
                summary.append(dict(parameter_i=r['parameter_i'],parameter_j=r['parameter_j'],O2_pct=o2,quantity=c,
                    minimum=min(values),maximum=max(values),fixed_background='seed25; global bounds',
                    interpretation='conditional response or reference contrast; not marginalized S2'))
    table(dest/'conditional_surface_ranges.tsv',summary)
    save(dest/'reuse_receipt.json',dict(status='complete',source_stage='full_N4096 conditional surfaces',
        fixed_background='seed25',n_pairs=6,n_surfaces=24,n_response_points=sum(int(r['n_points']) for r in rows),
        model_calls=0,neighborhood_files_read=0,files=files,completed_utc=now()))
    (dest/'README_zh.md').write_text(
        '# 复用的全域条件响应面\n\n'
        '这六对响应面直接复用既有计算，没有重算模型。响应面不依赖Sobol的N，因此可用于解释当前N131072分析，'
        '但保留原始full_N4096来源，不改称N131072模型采样结果。\n\n'
        '每张图展示固定氧条件下两参数变化时的倍性、增长率及参考点非加性对比。其余参数固定为seed25参考值；'
        '显示氧点为0、0.5、2.5和5%。完整敏感性分析仍使用201个氧点。\n\n'
        '先读坐标上的参数及尺度，再读色条输出或非加性对比。相同颜色的等值线表示不同参数组合能产生相似输出；'
        '若要讨论补偿，必须沿等值线说明两参数如何共同变化，以及固定的背景和氧点。'
        '非加性对比的正负相对于所标参考点定义，不是Sobol S2的正负，也不是普遍的生物学协同或拮抗。\n\n'
        '本次覆盖lam_max与p_misseg、p_mis_base、O2_crit、mu_hp的四对，以及buffer_smax与buffer_beta、'
        'p_mis_base与buffer_smax两对。没有响应面的强S2参数对仍只能讨论交互大小，不能据此补写方向。'
        '更高阶交互和TII不在本轮范围。\n')
    print('Verified and reused six global response pairs;',len(files),'source/figure files')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,required=True);parser.add_argument('--root',type=Path,required=True)
    a=parser.parse_args();reuse(a.source,a.root)
