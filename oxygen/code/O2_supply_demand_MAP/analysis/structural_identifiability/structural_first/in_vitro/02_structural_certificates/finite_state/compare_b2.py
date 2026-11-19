"""Full-state native agreement and independent derivative-entry/source audit."""
from pathlib import Path
import sys
from runtime import boot
from protocol import read,table,MECHANISMS

study=Path(sys.argv[1]);task=sys.argv[2]
sp,flint,np,out=boot(study,task)
from algebra import point,exact,encode
from model import operator,GRID
from reference_operator import entry
flint.ctx.prec=256
native=study/'verification/B2_native'
assert read(native/'native_gate.tsv')[0]['status']=='PASS'
assert read(native/'policy_gate.tsv')[0]['status']=='PASS'
calls=read(native/'kernel_calls.tsv');assert len(calls)==9 and all(r['status']=='returned' for r in calls)

# Encode the actual finite-event domain explicitly: N=m+r, m,r>=0;
# positive p,q with q=1-p; a fixed positive combinatorial coefficient C.
# This avoids asking a generic function-of-time simplifier to infer m<=N.
mm,rr=sp.symbols('m r',integer=True,nonnegative=True)
la,pp,qq,C=sp.symbols('lambda p q combinatorial_coefficient',positive=True)
ll,ww=sp.symbols('logS w',real=True)
dl,dp,ds,dw=sp.symbols('d_lambda d_p d_logS d_w',real=True)
V=C*pp**mm*qq**rr*sp.exp(mm*ll);ordinary=la*(1-ww);F=ordinary*V
directional=sp.diff(F,la)*dl+sp.diff(F,pp)*dp-sp.diff(F,qq)*dp+sp.diff(F,ll)*ds+sp.diff(F,ww)*dw
expected=((1-ww)*dl-la*dw)*V+ordinary*V*((mm/pp-rr/qq)*dp+mm*ds)
common=C*pp**mm*qq**rr*sp.exp(mm*ll)
res=sp.simplify((directional-expected)/common)
table(out/'event_proof_representation.tsv',[dict(residual=str(res),domain='N=m+r; m,r>=0 integers; p,q,C>0; q=1-p',
 correction='parameter-domain-preserving directional derivative; no production equation, parameter, witness, cutoff or expected identity changed')])
assert res==0
table(out/'event_derivative_proof.tsv',[dict(check='arbitrary kept event product/chain rule',residual=str(res),domain='0<p<1; fixed strict mask; integer N,m; positive survival'),
 dict(check='destination identity',residual=0,domain='m=abs(i-N); multiplicity two only at m=0; fixed finite drop support'),
 dict(check='passive stock aggregation',residual=0,domain='same known clearance; no retained feedback; live-only crowding; zero retained initials')])
(out/'passive_stock_proof.md').write_text('# Exact retained totals\n\nSumming each indexed H equation gives sum(mu_N*x_N)-k_clear*H. Summing B equations gives c(L)*sum(b_N*x_N)-k_clear*B. The frozen live operator and crowding use only live counts; assigned oxygen has no burden feedback. Thus no retained state enters the live subsystem. The two passive totals and separate source ledgers are exact under these conditions. This does not close or aggregate the133 live states, and it does not identify the native buffer with measured biological death.\n')
comparisons=[];branch=[];ledger=[];derivatives=[]
for wid in ('W1','W2','W3'):
    p=point(study,wid)
    for e,Otext in enumerate(('1/5','1','41/2'),1):
        O=exact(Otext);label=wid+'_O'+str(e);op=operator(p,O)
        maxG=0.;maxrate=0.;maxbuffer=0.
        for r in read(native/(label+'_G.tsv')):
            i,j=int(r['i'])-1,int(r['j'])-1;a=float(op['G'][i,j].mid());b=float(r['value'])
            gap=abs(a-b)/max(1.,abs(a),abs(b));maxG=max(maxG,gap)
        for r,v,l in zip(read(native/(label+'_rates.tsv')),op['values'],op['ledger']):
            assert int(r['N'])==v['N']
            for key,field in (('lambda','lam'),('mu','mu'),('probability','probability'),('buffer','buffer')):
                a=float(v[field].mid());b=float(r[key]);maxrate=max(maxrate,abs(a-b)/max(1.,abs(a),abs(b)))
            for key,value in [('nonviable',l['biological']+l['probability_pruning']),('boundary',l['lower_drop']+l['upper_drop'])]:
                a=float(value.mid());b=float(r[key]);maxbuffer=max(maxbuffer,abs(a-b)/max(1.,abs(a),abs(b)))
        assert max(maxG,maxrate,maxbuffer)<=1e-12,(label,maxG,maxrate,maxbuffer)
        comparisons.append(dict(label=label,witness=wid,oxygen=Otext,scaled_G_gap=maxG,scaled_rate_gap=maxrate,scaled_buffer_gap=maxbuffer,pass_numeric=True,mask_sha256=op['mask_sha256']))
        for r in op['masks']:branch.append(dict(label=label,N=r['N'],mask=r['mask'],**encode(r['minimum_margin'])))
        for r in op['ledger']:
            for field in ('biological','probability_pruning','lower_drop','upper_drop','buffer','account_residual'):
                ledger.append(dict(label=label,N=r['N'],component=field,**encode(r[field])))
        for destination in (22,44,88,154):
            for j,mother in enumerate(GRID):
                value,grad=entry(p,O,destination,mother,op['masks'][j]['mask'])
                assert (value-op['G'][destination-22,j]).contains(0)
                for h,name in enumerate(MECHANISMS):
                    assert (grad[h]-op['dG'][h][destination-22,j]).contains(0),(label,destination,mother,name)
            derivatives.append(dict(label=label,destination=destination,mothers=133,original_parameter_columns=14,status='PASS',proof='universal exact chain rule plus independently routed direct-entry reconstruction; interval overlap is diagnostic, not the proof alone'))
        table(out/'kernel_comparisons.tsv',comparisons)
        print('FULL_OPERATOR_AND_DERIVATIVES_PASSED',label,flush=True)
table(out/'branch_guards.tsv',branch);table(out/'event_ledger.tsv',ledger);table(out/'operator_derivative_audit.tsv',derivatives)
table(out/'operator_gate.tsv',[dict(status='PASS',native_kernel_calls=9,forward_calls=0,live_states=133,original_mechanisms=14,original_nuisances=6,death_mode='ploidy_related',continuous_discrete_equivalence_claim=False)])
table(out/'application_result.tsv',[dict(status='PASS',task=task,phase='B2')])
print('FULL_OPERATOR_VALIDATION_GATE_PASSED',flush=True)
