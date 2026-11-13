#!/usr/bin/env python3
"""Freeze an authorized new study from clean numerical imports; no model calls."""
from pathlib import Path
import hashlib
import shutil
import sys
sys.dont_write_bytecode=True
import ledger as L
import workflow as W

project=Path(sys.argv[1]).resolve();study=Path(sys.argv[2]).resolve();code=Path(__file__).resolve().parent
for sub in ['implementation','manifests','state','tasks','logs','submissions','tables','tests','fits','populations','reserve','diagnostics','reports']:
    (study/sub).mkdir(exist_ok=True)
if (study/'manifests/control_manifest.tsv').exists():raise ValueError('Existing frozen study; do not reinitialize')
baseline=project/'oxygen/results/structural_identifiability/20261111_culture_v2_baseline_recovery_nojson'
settings=L.read_kv(baseline/'settings.tsv');settings.update(core_hour_cap=4000,created_utc=L.now(),authorized=True,
  baseline_relative_path=str(baseline.relative_to(project)),global_attempt_cap=36,local_attempt_cap=70,reserve_round_cap=2)
L.write_kv(study/'settings.tsv',settings)
shutil.copy2(project/'docs/culture_seven_point_profile_nojson_plan.md',study/'manifests/approved_plan.md')
L.write_kv(study/'manifests/authorization.tsv',dict(authorized=True,budget=4000,nonfit_cap=120,qos='xxlarge',tumor_fitting=False,utc=L.now(),scope='Seven fixed culture points, guards, bidirectional continuation, at most two reserve rounds'))
ids=['guard_1','guard_2','diagnostic_base']
cases=L.read_table(study/'cases.tsv')
for c in cases:
    p=c['point'];ids.extend([p+'_historical',p+'_repeat',p+'_neighbor_forward',p+'_neighbor_reverse','diagnostic_'+p])
    for fam in ['A','B']:
        for k in [1,2]:ids.extend([p+'_'+fam+str(k)+'_global',p+'_'+fam+str(k)+'_local'])
    for n in [1,2]:
        for fam in ['A','B']:ids.extend([p+'_R'+str(n)+'_'+fam+'_global',p+'_R'+str(n)+'_'+fam+'_local'])
for n in [1,2]:
    for fam in ['A','B']:ids.extend(['base_R'+str(n)+'_'+fam+'_global','base_R'+str(n)+'_'+fam+'_local'])
    ids.extend('reserve_'+str(n)+'_repair_'+str(k) for k in range(1,9))
for p in ['base','pm_00','bb_00','pm_05','bb_07']:ids.extend(['fixture_'+p+'_global','fixture_'+p+'_local'])
assert len(ids)==len(set(ids))
seeds=[dict(task_id=tid,seed=10000000+i*1009) for i,tid in enumerate(sorted(ids))]
L.write_table(study/'seeds.tsv',seeds,['task_id','seed'])
source_rows={r['path']:r for r in L.read_table(baseline/'manifests/source_manifest.tsv')}
def pin(path,role):
    path=Path(path);rel=str(path.relative_to(project))
    source_rows[rel]=dict(path=rel,bytes=path.stat().st_size,sha256=L.sha(path),role=role)
for r in L.read_table(study/'imports.tsv'):
    pin(project/r['source'],'historical_numerical_artifact');pin(study/r['relative_path'],'clean_import')
pin(baseline/'tables/recovery_gate.tsv','completed_baseline_gate')
pin(baseline/'tables/delivery_receipt.tsv','completed_baseline_delivery')
pin(project/settings['review_relative_path']/'tables/review_point_quality.tsv','historical_point_registry')
L.write_table(study/'manifests/source_manifest.tsv',list(source_rows.values()),['path','bytes','sha256','role'])
for p in code.iterdir():
    if p.is_file() and p.suffix in {'.R','.py','.sh'}:shutil.copy2(p,study/'implementation'/p.name)
L.write_table(study/'manifests/code_manifest.tsv',[dict(path=p.name,sha256=L.sha(p)) for p in sorted((study/'implementation').iterdir()) if p.is_file()],['path','sha256'])
control=['settings.tsv','cases.tsv','imports.tsv','seeds.tsv','manifests/approved_plan.md','manifests/authorization.tsv','manifests/source_manifest.tsv']
L.write_table(study/'manifests/control_manifest.tsv',[dict(path=p,sha256=L.sha(study/p)) for p in control],['path','sha256'])
state=L.State(study)
with state.locked():
    data=state.load()
    imported={r['import_id']:r for r in L.read_table(study/'imports.tsv')}
    for point in ['base','pm_00','bb_00','pm_05','bb_07']:
        imp=imported['base_reference' if point=='base' else 'historical_'+point]
        src=str((study/imp['relative_path']).relative_to(project))
        gt='fixture_'+point+'_global';t=W.task(study,gt,'global','fixture',point,'fixture_'+point,source=src,fixture=True);W.add(state,data,t)
        lt='fixture_'+point+'_local';t=W.task(study,lt,'local','fixture',point,'fixture_'+point,source=str((study/'tests/fixtures'/gt/'accepted_fit.rds').relative_to(project)),fixture=True);W.add(state,data,t)
    state.event(data,'study_frozen','No scientific optimization before validation');state.save(data)
L.audit(study);L.verify_code(study)
print('Frozen study:',study)
print('Code manifest SHA256:',L.sha(study/'manifests/code_manifest.tsv'))
