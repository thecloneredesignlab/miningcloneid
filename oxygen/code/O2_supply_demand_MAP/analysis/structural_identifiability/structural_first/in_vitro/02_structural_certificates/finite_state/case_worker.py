"""One registered case/witness; finite precision and candidate scope, no timer."""
from pathlib import Path
import sys
from runtime import boot
from protocol import read,table,sha,kv,MECHANISMS

study=Path(sys.argv[1]);task=sys.argv[2];case_id,wid=task.split('_')
assert case_id in ('B01','B02','B03','B04','B05','B06') and wid in ('W1','W2','W3')
sp,flint,np,out=boot(study,task)
from algebra import point,exact,rates
from model import operator,columns_for,output_rows,GRID
from reference_operator import entry
from certificates import nominate,submatrix,validate_and_save
accepted_phase=kv(study/'accounting/operator_acceptance.tsv')['phase']
assert read(study/'verification'/accepted_phase/'operator_gate.tsv')[0]['status']=='PASS'
assert sha(Path(__file__).parent/'model.py')==sha(study/'releases/r03/model.py')
assert sha(Path(__file__).parent/'algebra.py')==sha(study/'releases/r02/algebra.py')
case=next(r for r in read(study/'inputs/cases.tsv') if r['case']==case_id)
if case_id=='B05':
    prerequisite=read(study/'verification'/('B03_'+wid)/'case_result.tsv')[0]
    assert prerequisite['status']!='CERTIFIED_FULL_LOCAL'
if case_id=='B06':
    prerequisite=read(study/'verification'/('B02_'+wid)/'case_result.tsv')[0]
    assert prerequisite['status']!='CERTIFIED_FULL_LOCAL'
flint.ctx.cap=5
flow=read(study/'inputs/flow_grid.tsv');step=exact(next(r['value'] for r in read(study/'inputs/flow_settings.tsv') if r['key']=='grid_step'))
columns=columns_for(case);target=len(columns)-(2 if case_id=='B01' else 0)
attempts=[];selection=None;row_scales=None;column_scales=None;original_masks=None
best=0;det_calls=0;certified=False

if case['output']=='full' and case['unknown_initial']=='True':
    a,mu,s,z=sp.symbols('adjacent_center mean width log_normalizer',real=True)
    lp=[-(a+r-mu)**2/(2*s**2)-z for r in (-1,0,1)]
    assert sp.simplify(lp[0]+lp[2]-2*lp[1]+1/s**2)==0
    assert sp.simplify(a+s**2*(lp[2]-lp[0])/2-mu)==0
    table(out/'initial_family_global_proof.tsv',[dict(cohort=c,status='PROVED_FOR_IDEAL_FULL_INITIAL_VECTOR',
     inverse='s^2=-1/log(x_minus*x_plus/x_center^2); mu=center+s^2*log(x_plus/x_minus)/2',
     domain='positive width and Gaussian weights; three adjacent chromosome states; full ideal x(0); finite normalization cancels') for c in ('2N','4N')])

def make_reference(p,O,op):
    n=len(GRID);G=flint.arb_mat(n,n);dG=[flint.arb_mat(n,n) for h in MECHANISMS]
    mus=[];dmus=[[] for h in MECHANISMS]
    for j,mother in enumerate(GRID):
        value,grad=rates(p,flint.arb(mother)/44,O);mus.append(value[1])
        for h,name in enumerate(MECHANISMS):dmus[h].append(grad[1][name])
        for i,destination in enumerate(GRID):
            v,dv=entry(p,O,destination,mother,op['masks'][j]['mask']);G[i,j]=v
            for h in range(14):dG[h][i,j]=dv[h]
    return dict(G=G,dG=dG,mu=mus,dmu=dmus)

for precision in (128,256,512,1024):
    flint.ctx.prec=precision;p=point(study,wid);rows=[];labels=[];operators=[];masks=[]
    for Ot in case['oxygens'].split(','):
        O=exact(Ot);op=operator(p,O);operators.append((Ot,O,op));masks.append((Ot,op['mask_sha256']))
        for cohort in ('2N','4N'):
            cols,ll,rr=output_rows(p,case,cohort,O,op,flow,step)
            assert cols==columns
            for l in ll:l['case']=case_id;l['oxygen']=Ot
            rows.extend(rr);labels.extend(ll)
    assert len(rows)<=6000 and len(rows[0])<=19
    if original_masks is None:original_masks=masks
    else:assert original_masks==masks,'Witness changed branch during refinement'
    if selection is None:
        selection,row_scales,column_scales,sv,rn=nominate(rows,labels,columns,np)
        table(out/'candidate_registry.tsv',[dict(candidate=r['candidate'],rows=','.join(map(str,r['rows'])),columns=','.join(columns[j] for j in r['columns']),scope='frozen midpoint nomination only; not a structural rank') for r in selection] or [dict(candidate='none',rows='',columns='',scope='unresolved nomination')])
        table(out/'nominal_singular_values.tsv',[dict(index=i,value=v,nominal_threshold='1e-10 * leading singular value; nomination only') for i,v in enumerate(sv)])
        table(out/'row_registry.tsv',[dict(index=i,scale_power_two=row_scales[i],**l) for i,l in enumerate(labels)])
        table(out/'column_registry.tsv',[dict(index=j,parameter=name,scale_power_two=column_scales[j]) for j,name in enumerate(columns)])
        table(out/'branch_registry.tsv',[dict(oxygen=o,mask_sha256=h) for o,h in masks])
    if not selection:break
    refrows=None
    for candidate in selection:
        M=submatrix(rows,candidate,row_scales,column_scales)
        # The determinant probe itself is made exactly once by validate_and_save.
        if refrows is None:
            refrows=[]
            for Ot,O,op in operators:
                refop=make_reference(p,O,op)
                for cohort in ('2N','4N'):
                    _,_,rr=output_rows(p,case,cohort,O,refop,flow,step);refrows.extend(rr)
        reference=submatrix(refrows,candidate,row_scales,column_scales)
        det_calls+=1;assert det_calls<=8
        table(out/'determinant_call_ledger.tsv',[dict(ordinal=i+1,status='returned' if i<len(attempts) else 'started') for i in range(det_calls)])
        result=validate_and_save(out,candidate,M,reference,precision);attempts.append(result)
        table(out/'determinant_call_ledger.tsv',[dict(ordinal=i+1,status='returned') for i in range(det_calls)])
        table(out/'certificate_attempts.tsv',attempts)
        if result['status']=='CERTIFIED':
            best=max(best,result['size'])
            if best==target:certified=True;break
        print('MINOR_ATTEMPT',case_id,wid,precision,result['candidate'],result['size'],result['status'],flush=True)
    if certified or best and best==max(len(s['columns']) for s in selection) and best<target:break

status='CERTIFIED_RANK_16_WITH_EXACT_SYMMETRY' if certified and case_id=='B01' else 'CERTIFIED_FULL_LOCAL' if certified else 'CERTIFIED_RANK_LOWER_BOUND' if best else 'UNRESOLVED_REGISTERED_SCOPE'
table(out/'case_result.tsv',[dict(case=case_id,witness=wid,status=status,rank_lower_bound=best,
 original_mean_columns=len(columns),registered_target=target,determinant_evaluations=det_calls,max_order=case['max_order'],
 maximum_precision_attempted=precision,model='conditional continuous finite-support fixed-mask model',
 local_domain='neighborhood of declared interior witness and strict branch masks; not entire production parameter box',
 global_identifiability='only the separately proved ideal full-initial Gaussian inverse, where applicable',
 native_tool_initial_support='not_verified',Step_5=False,fits=0)])
table(out/'application_result.tsv',[dict(status='PASS',task=task,scientific_result=status,scope_exhaustion_is_nonidentifiability=False)])
print('REGISTERED_STRUCTURAL_CASE_COMPLETED',case_id,wid,status,best,flush=True)
