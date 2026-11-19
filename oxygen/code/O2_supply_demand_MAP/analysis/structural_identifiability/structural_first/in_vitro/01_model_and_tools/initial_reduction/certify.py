"""Exact finite-jet certificates for the registered rational model only."""
from pathlib import Path
import csv,sys,time
from math import isqrt
import exact_model as m
def tab(p,rows,fields):
 p.parent.mkdir(parents=True,exist_ok=True)
 with p.open('w',newline='') as f:
  w=csv.DictWriter(f,fields,delimiter='\t');w.writeheader();w.writerows(rows)
def prime(n):return n>1 and all(n%d for d in range(2,isqrt(n)+1))
def run(study,task):
 out=study/'results'/task;out.mkdir(exist_ok=True);summary=[];individual=[];start=time.monotonic()
 assert all(prime(p) for p in m.PRIMES)
 for case in m.cases():
  results=[]
  for p in m.PRIMES:
   print('EXACT_CASE',case['case'],p,flush=True);r=m.evaluate_case(case,p);results.append(r)
   target=out/'certificates'/(case['case']+'_'+str(p));target.mkdir(parents=True)
   tab(target/'jacobian.tsv',[dict(row=i,label=label,**dict(zip(r['columns'],row))) for i,(label,row) in enumerate(zip(r['row_labels'],r['jacobian']))],['row','label']+r['columns'])
   tab(target/'minor.tsv',[dict(row_order=i,original_row=row,**{'c'+str(j):r['jacobian'][row][col] for j,col in enumerate(r['pivot_cols'])}) for i,row in enumerate(r['pivot_rows'])],['row_order','original_row']+['c'+str(j) for j in range(r['rank'])])
   for i,H in enumerate(r['hankels']):
    assert m.rank_pivots(H,p)[0]==5,'Linear growth output is not of generic minimal order five'
    tab(target/('hankel_'+str(i)+'.tsv'),[dict(zip(['c'+str(j) for j in range(5)],row)) for row in H],['c'+str(j) for j in range(5)])
   tab(target/'receipt.tsv',[dict(key=k,value=v) for k,v in dict(prime=p,rank=r['rank'],mathematical_upper_bound=r['upper'],nuisance_rank=r['nuisance_rank'],parameter_rank=r['parameter_rank'],nonzero_pivot_product=r['nonzero_pivot_product'],pivot_columns=','.join(str(v) for v in r['pivot_cols']),finite_jet_order=16,floating_point_derivative_operations=0).items()],['key','value'])
  a,b=results;assert a['rank']==b['rank'] and a['nuisance_rank']==b['nuisance_rank']
  complete=a['rank']==a['upper'];kind=case['kind'];single=len(case['oxygens'])==1
  overall='exact_generic_local_certificate' if complete else 'exact_rank_lower_bound_only_unresolved'
  summary.append(dict(case=case['case'],model=case['model'],crowding=case['crowding'],outputs=kind,oxygen_levels=';'.join(case['oxygens']),initial_unknown=case['unknown'],initial_type=case['initial'],cohorts=case['cohorts'],variant=case['variant'],joint_rank=a['rank'],joint_upper_bound=a['upper'],nuisance_rank=a['nuisance_rank'],mechanism_rank=a['parameter_rank'],status=overall,global_status='Unresolved except separately proved explicit symmetries and invariants'))
  for name in m.active(case['variant']):
   if not complete:status='unresolved';evidence='A nonzero minor is a lower bound; no maximal-rank claim'
   elif single and name in ('alpha','dH','K'):status='structurally_nonidentifiable';evidence='Explicit fixed-oxygen response symmetry'
   elif kind=='L' and case['unknown'] and single and not case['crowding']:
    status='structurally_nonidentifiable' if a['augmented'][name]>a['rank'] and b['augmented'][name]>b['rank'] else 'unresolved'
    evidence='Characteristic-polynomial fiber with initial-state remapping; exact coordinate-augmentation witness'
   else:status='locally_identifiable_global_unresolved';evidence='Exact maximal finite-jet minor and exhaustive known continuous-symmetry dimension'
   individual.append(dict(case=case['case'],model=case['model'],parameter=name,status=status,evidence=evidence,scope='Five-class rational model; not original fourteen-parameter simulator'))
  tab(out/'case_results.tsv',summary,list(summary[0]));tab(out/'parameter_results.tsv',individual,list(individual[0]))
  print('CERTIFIED',summary[-1],flush=True)
 tab(out/'result.tsv',[dict(key=k,value=v) for k,v in dict(pass_gate=True,cases=len(summary),exact_maximal_cases=sum(r['status']=='exact_generic_local_certificate' for r in summary),seconds=time.monotonic()-start,synthetic_datasets=0,fit_calls=0,profiles=0,tumor_fits=0,step05_release=False).items()],['key','value'])
 (out/'COMPLETE').write_text('Exact certificates produced; independent minor verification is required.\n')
if __name__=='__main__':run(Path(sys.argv[1]),sys.argv[2])
